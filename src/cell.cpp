/*© 2023. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos
National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S.
Department of Energy/National Nuclear Security Administration. All rights in the program are.
reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear
Security Administration. The Government is granted for itself and others acting on its behalf a
nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare.
derivative works, distribute copies to the public, perform publicly and display publicly, and to permit.
others to do so.*/

#include "cell.h"

#include "utilities.h"
#include "constants.h"
#include "configuration.h"
#include "network.h"
#include "cellobserver.h"
#include "elements.h"
#include "sput_params.h"
#include "sputter.h"

#include <iomanip>
#include <limits>
#include <boost/math/quadrature/gauss.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <boost/numeric/odeint.hpp>
#include <boost/format.hpp>
#include <plog/Log.h>
#include <algorithm>
#include <boost/format.hpp>
#include <boost/math/interpolators/makima.hpp>
#include <boost/range/algorithm.hpp>
#include <boost/range/numeric.hpp>
#include <boost/range/adaptors.hpp>
#include <boost/range/algorithm_ext.hpp>
#include <chrono>
#include <cmath>
#include <fstream>
#include <numeric>
#include <plog/Log.h>
#include <string>
#include <valarray>
#include <vector>
#include <sstream>

using boost::math::interpolators::makima;
using namespace boost::numeric::odeint;
using namespace std::chrono;

// defined parameters needed later and call funcitons to initialize data with input data
cell::cell ( network* n, params* sputARR, configuration* con, uint32_t id, const spec_v &init_s,
             const cell_input &input_data )
    : net (n), sputARR (sputARR), config (con), cid (id), elm("data/elements.json")
{
  using constants::N_MOMENTS;
  cell_st.numReact = net->n_nucleation_reactions;
  cell_st.numBins = config->bin_number;
  cell_st.numGas = init_s.size();
  set_init_data(init_s, input_data);
  set_env_data(input_data);
  cell_st.parts.resize(cell_st.numReact);
  reaction_switch.resize(cell_st.numReact);
  std::fill(reaction_switch.begin(), reaction_switch.end(), true);
}

// resize and set initial data 
void
cell::set_init_data(const spec_v& init_s, const cell_input& init_data)
{
  using constants::N_MOMENTS;
  
  // solution vector and abundances
  cell_st.abund_moments_sizebins.assign(init_data.inp_solution_vector.begin(),init_data.inp_solution_vector.end());
  cell_st.init_abund.resize(cell_st.numGas);
  for (size_t i = 0; i < cell_st.numGas; ++i) {
    auto idx = net->get_species_index(init_s[i]);
    if (idx != -1) {
      cell_st.init_abund[idx] = init_data.inp_init_abund[i];
    }
  }
  // vectors for tracking multiple grain values. makes printout easier
  cell_st.cbars.resize(net->n_nucleation_reactions);
  cell_st.S.resize(net->n_nucleation_reactions);
  cell_st.Js.resize(net->n_nucleation_reactions);
  cell_st.dadt.resize(net->n_nucleation_reactions);
  cell_st.ncrit.resize(net->n_nucleation_reactions);

  cell_st.start_time = init_data.sim_start_time;
  // vectors for binning and destruction/growth
  cell_st.grn_sizes.assign(init_data.inp_binSizes.begin(),init_data.inp_binSizes.end());
  cell_st.edges.assign(init_data.inp_binEdges.begin(),init_data.inp_binEdges.end());
  cell_st.vd.assign(init_data.inp_vd.begin(),init_data.inp_vd.end());
  cell_st.rebin_chng.resize(cell_st.numBins * cell_st.numReact);
  cell_st.runningTot_size_change.assign(init_data.inp_delSZ.begin(),init_data.inp_delSZ.end());

}

// load data and setup interpolators
void
cell::set_env_data(const cell_input& input_data)
{
  if(config->environment_file.empty())
  {
    cell_st.temperature = input_data.inp_shock_temp;
    PLOGD << "No environment file loaded.";
    return;
  }
  
  env_times.assign(input_data.inp_times.begin(), input_data.inp_times.end());
  env_temp.assign(input_data.inp_temp.begin(),input_data.inp_temp.end());
  env_volumes.assign(input_data.inp_volumes.begin(),input_data.inp_volumes.end());
  env_rho.assign(input_data.inp_rho.begin(), input_data.inp_rho.end());
  env_pressure.assign(input_data.inp_pressure.begin(),input_data.inp_pressure.end());

  cell_st.temperature = input_data.inp_shock_temp;
  cell_st.rho      = env_rho[0];
    
  cell_st.volume_0      = env_volumes[0];
  cell_st.volume        = env_volumes[0];

  // return if there is only one initial time, otherwise continue and setup splines
  if(input_data.inp_times.size()==1)
  {
    return;
  }

  std::vector<double> times = env_times;
  env_temp_interp = makima(std::move(times), std::move(env_temp));
  std::vector<double> times1 = env_times;
  env_volume_interp = makima(std::move(times1), std::move(env_volumes));
  std::vector<double> times2 = env_times; 
  env_rho_interp     = makima(std::move(times2), std::move(env_rho));
  std::vector<double> times5 = env_times;
  env_pressure_interp = makima(std::move(times5), std::move(env_pressure));
}

// check there's no nans or negatives in the solution
bool
cell::check_solution(const std::vector<double>& x)
{
  for (const auto& val: x) {
    if (val < 0.0 || std::isnan(val))
      return false;
  }
  return true;
}

// setup and start integrator, initialize cellObserver which deals with data output. There're a few checks to stop integration. 
void
cell::solve()
{
  
  using constants::k_B;
  using constants::kB_eV;
  using numbers::one;
  using constants::yrs2s;
  
  auto abs_err = config->ode_abs_err, rel_err = config->ode_rel_err;
  double max_dt = config->ode_dt_max;// min_dt = config->ode_dt_min;
  double time_start, time_end;
  if(!config->environment_file.empty() && env_times.size()!=1)
  {
    time_start = env_times[0];
    time_end = env_times.back();
    PLOGI << "Times taken from environment file.";
    PLOGI << "Start time: " << time_start << ", end time: " << time_end;
  }
  else
  {
    time_start = cell_st.start_time;
    time_end = cell_st.start_time + 3.14e7;
    PLOGI << "End Time is not specified, running simulation out for another year. End Time: " << time_end;
  }
  auto dt0             = config->ode_dt_0;
  auto rkd             = runge_kutta_dopri5<std::vector<double>>{};
  auto stepper         = make_dense_output(abs_err, rel_err, max_dt, rkd);
  auto dumpN           = config->io_dump_n_steps;
  auto RSN             = config->io_restart_n_steps;
  size_t n_solve_steps   = 0;
  size_t n_stepper_reset = 0;
  cell_st.kT = k_B * cell_st.temperature; // ergs
  cell_st.kTeV = kB_eV * cell_st.temperature;
  cell_st.invkT = 1.0 / cell_st.kT;
  stepper.initialize(cell_st.abund_moments_sizebins, time_start, dt0);
  calc_state_vars(cell_st.abund_moments_sizebins, time_start);
  CellObserver observer(cid,net,config);
  observer.init_dump(cell_st);
  
  while ((stepper.current_time() < time_end)) {
    auto t0               = stepper.current_time();
    auto dt               = stepper.current_time_step();
    auto x0               = stepper.current_state();
    cell_st.time             = t0;
    cell_st.dt               = dt;
    integration_abandoned = false;
    if (t0 + dt > time_end) {
      PLOGI << "finished integration cell: " << cid << ", t_current: " << t0;
      break;
    }
    stepper.do_step(std::ref(*(this)));
    check_reactions(stepper.current_state());
    if (integration_abandoned) {
      stepper.initialize(x0, t0, dt * 0.5);
      ++n_stepper_reset;
    } 
    else 
    {
      observer(cell_st);
      n_stepper_reset = 0;
    }
    if (n_solve_steps > CELL_MAX_STEPS) {
      PLOGI << "too many solve steps, exiting cell " << cid << " at t: " << stepper.current_time();
      break;
    }
    if (n_stepper_reset > CELL_MAXIMUM_STEPPER_RESETS) {
      PLOGI << "TOO MANY RESTARTS, exiting cell " << cid << " at t = " << stepper.current_time();
      break;
    }
    // the user can specify restart and data dumps here or in the cellObserver operator
    /*
    if(n_solve_steps%dumpN==0.0)
    {
      observer.dump_data(cell_st);
    }
    if(n_solve_steps%RSN==0.0)
    {
      observer.restart_dump(cell_st);
    }
    */
    ++n_solve_steps;
  }
  PLOGI << "done cell: " << cid;
  observer.finalSave(cell_st);
}

// check the abundances are greater than the min abundance
void
cell::check_reactions(const std::vector<double>& x)
{
  for (size_t i = 0; i < cell_st.numReact; ++i) {
    reaction_switch[i] = true;
    for (const auto& r_idx: net->reactants_idx[i]) {
      if (x[r_idx] < CELL_MINIMUM_ABUNDANCE) {
        reaction_switch[i] = false;
        break;
      }
    }
  }
}

// update state variables from interpolator or if no spline was created, update cell temperature
void
cell::calc_state_vars(const std::vector<double>& x, const double time)
{
  using constants::k_B;
  using constants::kB_eV;

  if(!config->environment_file.empty() && env_times.size()!=1)
  {
    cell_st.temperature = env_temp_interp(time);
    cell_st.volume   = env_volume_interp(time);
    cell_st.rho      = env_rho_interp(time);
    cell_st.drho     = env_rho_interp.prime(time);
    cell_st.pressure = env_pressure_interp(time);
    cell_st.dP       = env_pressure_interp.prime(time);
  }
  
  cell_st.kT = k_B * cell_st.temperature; // ergs
  cell_st.kTeV = kB_eV * cell_st.temperature;
  cell_st.invkT = 1.0 / cell_st.kT;

}

// solve ODEs for nucleation, grain growth, key species depeletion, etc.
void cell::nucleate(const std::vector<double>& x)
{
  using constants::amu2g;
  using constants::pi;
  using constants::istdP;
  using constants::stdP;
  using constants::N_MOMENTS;

  int sd_start = cell_st.numGas + cell_st.numReact * N_MOMENTS;
  for (size_t gidx = 0; gidx < cell_st.numReact; ++gidx) 
  {
    // finding key specie for reaction and identifying the reactants and their index
    auto reaction_idx = net->nucleation_reactions_idx[gidx];
    auto num_ks       = net->ks_lists_idx[gidx].size();
    auto key_spec_idx = net->ks_lists_idx[gidx][0];
    if (num_ks > 1) 
    {
      for (size_t kidx = 0; kidx < num_ks; ++kidx) 
      {
        auto r_idx = net->ks_lists_idx[gidx][kidx];
        if (x[key_spec_idx] > x[r_idx])
        {
          key_spec_idx = r_idx;
        }
      }
    }
    cell_st.parts[gidx].ks_idx = key_spec_idx;
    // initializing and zeroing nuclation arrays
    if (x[key_spec_idx] < CELL_MINIMUM_ABUNDANCE) 
    {
      cell_st.parts[gidx].is_nucleating = 0;
      cell_st.S[gidx]                   = 0.0;
      cell_st.Js[gidx]                  = 0.0;
      cell_st.dadt[gidx]                = 0.0;
      cell_st.ncrit[gidx]               = 0.0;
      continue;
    }
    cell_st.parts[gidx].ks_react_mass =
      elm.elements.at(net->species[cell_st.parts[gidx].ks_idx]).mass * amu2g;
    std::vector<double> react_nu;
    std::vector<int> react_idx;
    // stoichiometry calculations for each reaction
    auto stoich_ks = 0.0;
    for (const auto& kv: net->nucleation_species_count[gidx]) 
    {
      if (kv.first == key_spec_idx) 
      {
        stoich_ks = kv.second;
      }
      react_idx.emplace_back(kv.first);
    }
    for (const auto& kv: net->nucleation_species_count[gidx])
    {
      react_nu.emplace_back(kv.second / stoich_ks);
    }
    cell_st.parts[gidx].react_idx = react_idx;
    cell_st.parts[gidx].react_nu  = react_nu;
    // nozawa et al. 2003 equ. 4, 2nd term r.h.s.
    // term for saturation
    double psum = 0.0;
    for (size_t ridx = 0; ridx < react_idx.size(); ridx++) 
    {
      if (x[react_idx[ridx]] != 0) 
      {
        if (react_idx[ridx] != cell_st.parts[gidx].ks_idx)
        {
          psum = psum + log(x[react_idx[ridx]] * cell_st.kT * istdP) * react_nu[ridx];
        }
      }
    }
    // updating concentrations
    double c1 = x[key_spec_idx];
    cell_st.cbars[gidx] = cell_st.init_abund[key_spec_idx] * cell_st.volume_0 / cell_st.volume;
    // change in Gibbs free energy 
    auto delg_reduced = (net->reactions[reaction_idx].alpha / cell_st.temperature -
                          net->reactions[reaction_idx].beta) + psum;
    // saturation
    // nozawa et al. 2003 equ 4
    cell_st.parts[gidx].lnS  = log(c1 * cell_st.kT * istdP) + delg_reduced;
    // weights from reaction
    double w = 1.0;
    for (size_t ridx = 0; ridx < react_idx.size(); ++ridx) 
    {
      if (react_idx[ridx] != cell_st.parts[gidx].ks_idx)
      {
        w = w + react_nu[ridx];
      }
    }
    // function of partial gas presures
    // yamamoto et al 2001 equ 16
    double Pii = 1.0;
    for (size_t ridx = 0; ridx < react_idx.size(); ++ridx) 
    {
      if (react_idx[ridx] != cell_st.parts[gidx].ks_idx)
      {
        Pii = Pii * std::pow(x[react_idx[ridx]] / c1, react_nu[ridx]);
      }
    }
    if (cell_st.parts[gidx].lnS > 0.0) 
    {
      double iw = 1. / w;
      Pii       = std::pow(Pii, iw);
      // nozawa et al. 2003 energy barrier for nucleation
      double mu = 4.0 * pi * std::pow(net->reactions[reaction_idx].a_rad, 2.) *
                  net->reactions[reaction_idx].sigma / cell_st.kT;
      // nozawa et al. 2003 equ 3 term in exponential
      double expJ = -4.0 / 27.0 * std::pow(mu, 3.) / std::pow(cell_st.parts[gidx].lnS, 2.);
      // nozawa et al. 2003 equ 3 term in 1st square root r.h.s.
      double Jkin = std::pow(2.0 * net->reactions[reaction_idx].sigma /
                          (pi * cell_st.parts[gidx].ks_react_mass),0.5);
      // saturation nozawa et all 2003 exponential of equ 4 
      cell_st.S[gidx]    = exp(cell_st.parts[gidx].lnS);
      // steady state nucleation rate nozawa et al. 2003 equ 3
      cell_st.Js[gidx]   = net->reactions[reaction_idx].omega0 * Jkin * c1 * c1 * Pii * exp(expJ);
      // growth rate, nozawa et al. 2003 equ 8
      cell_st.dadt[gidx] = net->reactions[reaction_idx].omega0 *
                      std::pow(0.5 * cell_st.kT / (pi * cell_st.parts[gidx].ks_react_mass), 0.5) *
                      c1 * (1. - 1. / cell_st.S[gidx]);
      // critical radius nozawa et al. 2003
      cell_st.ncrit[gidx] = std::pow(2.0 / 3.0 * (mu / cell_st.parts[gidx].lnS), 3.0) + iw;
      // finding growth from the dadt and storing it to determine if rebinning is needed
      double growth = cell_st.dadt[gidx] * cell_st.dt;
      for (int bidx = 0; bidx < cell_st.numBins; ++bidx)
      {
        auto idx = (gidx*cell_st.numBins)+bidx;
        if(cell_st.abund_moments_sizebins[idx+sd_start]==0.0) continue;
        cell_st.runningTot_size_change[idx] += growth;
      }
    } 
    else 
    { // we want to force these to zero in case a value is unchanged for
      // the next timestep
      cell_st.S[gidx]     = 0.0;
      cell_st.Js[gidx]    = 0.0;
      cell_st.dadt[gidx]  = 0.0;
      cell_st.ncrit[gidx] = 0.0;
    }
    if (cell_st.ncrit[gidx] > 0.0) 
    {
      cell_st.parts[gidx].grains_nucleating = cell_st.Js[gidx] * cell_st.ncrit[gidx];
      cell_st.parts[gidx].is_nucleating     = 1;
    } 
    else 
    {
      cell_st.parts[gidx].grains_nucleating = 0.0;
      cell_st.parts[gidx].is_nucleating     = 0;
    }
  }
}

// checking if a grain nucleates, finds the size and add it to the solution vector
void cell::add_new_grn(const std::vector<double>& x)
{
  using constants::N_MOMENTS;
  int sd_start = cell_st.numGas + cell_st.numReact * N_MOMENTS;
  for (int gidx = 0; gidx < cell_st.numReact; ++gidx) 
  {
    if (cell_st.parts[gidx].lnS > 0.0) 
    {
      auto momIDX = cell_st.numGas + constants::N_MOMENTS * gidx;
      if ((x[momIDX + 3] > 0.0) && (x[momIDX + 0] > 0.0)) 
      {
        auto reaction_idx = net->nucleation_reactions_idx[gidx];
        // new grain size nozawa et al. 2003 equation 12
        double new_grn_size =
          (net->reactions[reaction_idx].a_rad) *
          std::pow(x[momIDX + 3] / x[momIDX + 0], 1. / 3.);
        auto addToBin = 0;
        double dr    = 0;
        for (size_t bidx = 0; bidx < cell_st.grn_sizes.size(); ++bidx) 
        {
          if ((new_grn_size < cell_st.edges[bidx + 1]) && (new_grn_size > cell_st.edges[bidx])) 
          {
            addToBin = bidx;
            dr       = cell_st.edges[bidx + 1] - cell_st.edges[bidx];
          }
        }
        cell_st.abund_moments_sizebins[sd_start + gidx*cell_st.numBins+addToBin] += cell_st.Js[gidx] * cell_st.dt / dr;
      }
    }
  }
}

// calculate the sputtering yield for an impactor
double cell::Y(const double& E, const int grnid, const int gsID)
{
    using constants::pi;
    using constants::echarge_sq;
    using numbers::onehalf;
    using numbers::twothird;
    using utilities::square;
    using numbers::sp_3_441;
    using numbers::sp_2_718;
    using numbers::one;
    using numbers::sp_6_35;
    using numbers::sp_6_882;   
    using numbers::sp_1_708;
    using constants::y_pref;

    if(E<sputARR->Eth[grnid][gsID]){return 0.0;}
    // biscaro & cherchneff 2016 eq 7
    auto eps = sputARR->eiCoeff[grnid][gsID] * E;
    // sqrt(eps)
    auto sqrt_eps = std::sqrt(eps);
    // biscaro & cherchneff 2016 si(epsi) eq 6, matsunami et al 1980
    auto sieps = sp_3_441 * sqrt_eps * std::log(eps + sp_2_718)
                / (one + sp_6_35 * sqrt_eps + eps * (sp_6_882 * sqrt_eps - sp_1_708));
    // biscaro & cherchneff 2016 Eth/E
    auto eth_ratio = sputARR->Eth[grnid][gsID]/E;
    // biscaro & cherchneff 2016 equ 4
    auto Si = sputARR->SiCoeff[grnid][gsID] * sieps;
    // biscaro & cherchneff 2016 equ 2
    auto suffix = (one - std::pow(eth_ratio,twothird)) * square(one - eth_ratio);
    // biscaro & cherchneff 2016 equ 2
    auto preret = Si * sputARR->alpha[grnid][gsID] * suffix / 
                (sputARR->u0[grnid] * (sputARR->K[grnid] * sputARR->mu[grnid][gsID] + one));
    // sputtering yield biscaro & cherchneff 2016 equ 2
    return y_pref * preret;
}

// thermal sputtering rate caclulations
double cell::Therm(const int gidx, const int gsID)
{
    // nozawa et al 2006 equ 22
    double pref = sputARR->msp_2rhod[gidx]  * 
                    sqrt(sputARR->y8_piMi[gsID]*cell_st.kT);
    // nozawa et al 2006 equ 22
    auto f = [&](const double& x) { return x * exp(-x)*Y(x*cell_st.kTeV,gidx,gsID); };
    double lowLim = sputARR->Eth[gidx][gsID]/cell_st.kTeV;
    double Q = boost::math::quadrature::gauss_kronrod<double, 7>::integrate(f, 
                lowLim, std::numeric_limits<double>::infinity());
    // thermal sputtering rate. nozawa et al 2006 equ 22, dwek et al 1996
    return pref * Q * cell_st.abund_moments_sizebins[gsID]; 
}

// calculate the non-thermal sputtering rate
double cell::NonTherm(const double& sidx, const int gidx, const int gsID)
{
    using constants::JtoEV;
    using constants::kB_eV;
    using constants::cm2m;
    using constants::amu2g;
    using constants::amu2Kg;
    using numbers::onehalf;
    using utilities::square;
    int iid = gidx*cell_st.numBins+sidx;
    double x = onehalf * sputARR->miKG[gsID] * 
                square(cell_st.vd[iid]*cm2m)*JtoEV; 
    double pref = sputARR->msp_2rhod[gidx] * cell_st.vd[iid] * cell_st.abund_moments_sizebins[gsID];
    // nonthermal sputtering rate. nozawa et al 2006 equ 23
    return  pref * Y(x,gidx,gsID); 
}

// determin which sputtering occurs, clalculate it, store the erosion amount to determine if rebinning is needed.
void cell::destroy()
{
  using constants::pi;
  using constants::k_B;
  using constants::amu2g;
  using constants::JtoEV;
  using constants::kB_eV;
  using constants::amu2Kg;
  using constants::kBeta;
  using constants::cm2m;
  using numbers::onehalf;
  using utilities::square;
  using numbers::ten;
  using constants::N_MOMENTS;
  int sd_start = cell_st.numGas + cell_st.numReact * N_MOMENTS;

  for ( auto gidx = 0; gidx < cell_st.numReact; ++gidx )
  {
    for( int sidx =0; sidx < cell_st.numBins; ++sidx)
    {
        double dadt = 0.0;
        int idx = (gidx*cell_st.numBins)+sidx;
        // remember grain sizes are in cm
        if(cell_st.abund_moments_sizebins[sd_start + idx] != 0.0)
        {
            for(auto gsID=0; gsID < cell_st.numGas; ++gsID)
            {
                if(cell_st.abund_moments_sizebins[gsID] == 0.0) continue;
                // s_i2 is unitless, invkT is in cgs, vd is in cm/s
                double s_i2 = sputARR->miGRAMS[gsID] * onehalf * cell_st.invkT * square(cell_st.vd[idx]);
                // non-thermal sputtering occurrs
                if( s_i2 > ten) 
                {
                    dadt +=  NonTherm(sidx+sd_start,gidx,gsID);
                }
                // thermal sputtering occurs
                else 
                {
                    dadt += Therm(gidx,gsID);
                }
            }
        }
        cell_st.runningTot_size_change[idx] -= dadt*cell_st.dt;
        // calculate the slow down of the shock and update
        double temp_velo = calc_dvdt(cell_st.grn_sizes[sidx],cell_st.vd[idx],gidx) * cell_st.dt; // in cm/s
        if(cell_st.vd[idx] - std::abs(temp_velo) >= 0.0)
        {
            cell_st.vd[idx] = cell_st.vd[idx] - std::abs(temp_velo);
        }
        else{
            cell_st.vd[idx] = 0.0;
        }
    }
  }
}

// calculate the slowing of the shocked grains
double cell::calc_dvdt(const double& cross_sec, const double& vd, const int grnid) // cross sec is in cm
{
    using constants::pi;
    using constants::k_B;
    using constants::amu2g;
    using constants::amu2Kg;
    using numbers::eight_threeRootPi;
    using numbers::ninePi_sixyfour;
    using utilities::square;
    using numbers::one;
    using constants::N_MOMENTS;
    int sd_start = cell_st.numGas + cell_st.numReact * N_MOMENTS;

    // nozawa et al 2006 equ 19
    double G_tot = 0;
    for(int gsID=0; gsID < cell_st.numGas; ++gsID)
    {
        // units of # of particles * mass in grams. might just need the mass not the * # of particles
        double m = sputARR->miGRAMS[gsID]; // should be in grams now
        double s2 = m * square(vd) / (2.*cell_st.kT); // assumes cgs units
        G_tot += cell_st.abund_moments_sizebins[sd_start+gsID] * std::sqrt(s2) * eight_threeRootPi * 
                std::sqrt(1.+s2*ninePi_sixyfour);
    }
    auto yield = sputARR->three_2Rhod[grnid] * cell_st.kT/(cross_sec)*G_tot;
    // deceleration nozawa et al 2006 equ 19
    return  -std::abs(yield); // should be in cm/s, cross sec is in cm
}

// rebin grains based on the growth and erosion totals
void cell::rebin(const std::vector<double>& x, std::vector<double>& dxdt)
{
  using constants::kBeta;
  using constants::equPres;
  using constants::pi;
  using constants::istdP;
  using constants::k_B;
  using constants::stdP;
  using constants::amu2g;
  using constants::kB_eV;
  using constants::N_MOMENTS;
  int sd_start = cell_st.numGas + cell_st.numReact * N_MOMENTS;

  std::fill(cell_st.rebin_chng.begin(),cell_st.rebin_chng.end(),0.0);
  // now we find which grains move up, which move down, and which stay the same
  for ( auto gidx = 0; gidx < cell_st.numReact; ++gidx )
  {
    std::vector<double> binsMoveUp(cell_st.numBins,0.0);
    std::vector<double> binsMoveDown(cell_st.numBins,0.0);
    for (size_t bidx = 0; bidx < cell_st.numBins; ++bidx)
    {
      auto idx = (gidx*cell_st.numBins)+bidx;
      if (cell_st.abund_moments_sizebins[sd_start + idx]==0.0) continue;
      // move down a bin
      if(cell_st.grn_sizes[bidx]+cell_st.runningTot_size_change[idx] < cell_st.edges[bidx])
      {
        if(bidx==0)continue;
        else
        {
          double binWidth = cell_st.edges[bidx+1]-cell_st.edges[bidx];
          double lowerBinWidth = cell_st.edges[bidx]-cell_st.edges[bidx-1];
          binsMoveDown[bidx-1] = cell_st.abund_moments_sizebins[sd_start + idx]*(binWidth/lowerBinWidth);
          cell_st.runningTot_size_change[idx] = 0.0;
        }
        cell_st.rebin_chng[bidx] -= 1.0;
        cell_st.rebin_chng[bidx-1] += 1.0;
      }
      // move up a bin
      if(cell_st.grn_sizes[bidx]+cell_st.runningTot_size_change[idx] > cell_st.edges[bidx+1])
      {
        if(bidx==cell_st.numBins-1) continue;
        else
        {
          double binWidth = cell_st.edges[bidx+1]-cell_st.edges[bidx];
          double higherBinWidth = cell_st.edges[bidx+2]-cell_st.edges[bidx+1];
          binsMoveUp[bidx+1] = cell_st.abund_moments_sizebins[sd_start + idx]*(binWidth/higherBinWidth);
          cell_st.runningTot_size_change[idx] = 0.0;
        }
        cell_st.rebin_chng[bidx] -= 1.0;
        cell_st.rebin_chng[bidx+1] += 1.0;
      }
    }
    // now update the size bins
    for (size_t bidx = 0; bidx < cell_st.numBins; ++bidx)
    {
      auto idx = (gidx*cell_st.numBins)+bidx;
      dxdt[sd_start + idx] += binsMoveUp[bidx]+binsMoveDown[bidx];
    }
  }
}

// called by integrator, updates x, dxdt, calls the relevant calculations
void
cell::operator()(const std::vector<double>& x, std::vector<double>& dxdt, const double t)
{
  using constants::N_MOMENTS;
  double fi;
  std::fill(std ::begin(dxdt), std ::end(dxdt), 0.0);
  if (!check_solution(x)) {
    integration_abandoned = true;
    return;
  }
  calc_state_vars(x, t);
  if(config->do_nucleation==1)
  {
    nucleate(x);
  }
  if(config->do_destruction==1)
  {
    destroy();
  }
  rebin(x, dxdt);
  add_new_grn(x);
  for (size_t i = 0; i < cell_st.numReact; ++i) {
    if ((cell_st.parts[i].is_nucleating) && (cell_st.ncrit[i] > 2.0)) {
      auto gidx   = cell_st.numGas + N_MOMENTS * i;
      dxdt[gidx] = cell_st.Js[i] / cell_st.cbars[i];
      for (int j = 1; j < N_MOMENTS; ++j) {
        dxdt[gidx + j] =
          dxdt[gidx] * std::pow(cell_st.ncrit[i], (j / 3.0)) +
          (j / net->reactions[i].a_rad) * cell_st.dadt[i] * x[gidx + j - 1];
      }
      for (size_t idx = 0; idx < cell_st.parts[i].react_idx.size(); ++idx) {
        auto r_idx = cell_st.parts[i].react_idx[idx];
        auto r_nu  = cell_st.parts[i].react_nu[idx];
        dxdt[r_idx] -= cell_st.cbars[i] * dxdt[gidx + 3] * r_nu;
      }
    }
  }

  for (size_t i = 0; i < cell_st.numGas; ++i)
    dxdt[i] += cell_st.drho / cell_st.rho * x[i];
  for (size_t i = 0; i < cell_st.numReact; ++i) {
    auto reaction_idx = net->nucleation_reactions_idx[i];
    if (!reaction_switch[reaction_idx])
      continue;
    for (const auto& r: net->reactants_idx[reaction_idx])
      dxdt[r] -= cell_st.parts[i].grains_nucleating;
    for (const auto& p: net->products_idx[reaction_idx])
      dxdt[p] += cell_st.parts[i].grains_nucleating;
  }

  for (size_t i = 0; i < net->n_chemical_reactions; ++i) {
    auto reaction_idx = net->chemical_reactions_idx[i];
    if (!reaction_switch[reaction_idx])
      continue;
    fi = 1.0;
    for (const auto& r: net->reactants_idx[reaction_idx])
      fi *= x[r];
    fi *= net->reactions[reaction_idx].rate(cell_st.temperature);
    for (const auto& r: net->reactants_idx[reaction_idx])
      dxdt[r] -= fi;
    for (const auto& p: net->products_idx[reaction_idx])
      dxdt[p] += fi;
  }
}
