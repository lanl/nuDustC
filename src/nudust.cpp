/*© 2023. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos
National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S.
Department of Energy/National Nuclear Security Administration. All rights in the program are.
reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear
Security Administration. The Government is granted for itself and others acting on its behalf a
nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare.
derivative works, distribute copies to the public, perform publicly and display publicly, and to permit.
others to do so.*/

#include "nudust.h"

#include "constants.h"
#include "utilities.h"
#include "sputter.h"
#include "sput_params.h"

#include <vector>
#include <string>
#include <chrono>
#include <fstream>
#include <filesystem>
#include <boost/filesystem.hpp>

#include <boost/foreach.hpp>

#include <plog/Log.h>
#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>

namespace options = boost::program_options;

nuDust::nuDust ( const std::string &config_file, int sz, int rk) : par_size(sz), par_rank(rk), sputter("data/sputterDict.json")
{
    nu_config.read_config ( config_file );
    // these are always called
    load_network();
    load_initial_abundances();

    ////////////////////////////////////////////////
    // these are called depending on the config file
    // either make a new size bin or read in the size bin from a file
    int (nu_config.sizeDist_file.empty()) ? gen_size_dist() : load_sizeDist();
    // nucleation and destrcution + nucleation path

    if(!nu_config.environment_file.empty())
    {
        load_environment_data();
    }
    // destruction
    if(nu_config.do_destruction==1)
    {
        load_sputter_params();
        if(!nu_config.shock_file.empty())
        {
            // read in shock file
            load_shock_params();
        }
        if(!isnan(nu_config.shock_velo))
        {
            // create shock velo and temp arrays from user specifies shock temp and velocity
            gen_shock_array_frm_val();
        }
    }
    ////////////////////////////////////////////

    // always run but needed at the end
    generate_sol_vector();
    load_outputFL_names();
    create_simulation_cells();
}

// make sure the network is loaded
void
nuDust::load_network()
{
    net.read_network ( nu_config.network_file );
    net.post_process();  
    PLOGI << "loaded network file";
}

// return the index of the element in the abundance vector
int
nuDust::get_element_index(const std::string& elem) const
{
  auto it =
    std::find(std ::begin(initial_elements), std ::end(initial_elements), elem);
  assert(it != std ::end(initial_elements));
  auto idx = -1;
  idx      = std::distance(std ::begin(initial_elements), it);
  return idx;
}

// create a size distribution from the user specified size parameters
void
nuDust::gen_size_dist()
{
    // generating an empty size distribution from user specified boundaries and number of bins
    if(isnan(nu_config.low_sd_exp) || isnan(nu_config.high_sd_exp) || isnan(nu_config.bin_number))
    {
        std ::cout << "! Missing parameters needed to generate a binned size distribution && no size distribution file was specified. Try again.\n"; 
        exit(1);
    }
    double low = nu_config.low_sd_exp;
    double high = nu_config.high_sd_exp;
    numBins = nu_config.bin_number;
    auto expDel = (high - low)/static_cast<double>(numBins);
    init_bin_edges.resize(numBins+1);
    std::generate(init_bin_edges.begin(), init_bin_edges.end(), [&, n = 0] () mutable { return std::pow(10.,low+static_cast<double>(n++)*expDel); } );

    init_size_bins.resize(numBins);
    // this is compiling, double check it generates the correct bin sizes
    std::generate(init_size_bins.begin(), init_size_bins.end(), [&, n = 0,m=1] () mutable { return (init_bin_edges[static_cast<int>(n++)]+init_bin_edges[static_cast<int>(m++)])/2.0; } );    
    for ( const auto &ic : cell_inputs)
    {
        auto cell_id = ic.first;
        cell_inputs[cell_id].inp_binEdges.resize(numBins+1);
        cell_inputs[cell_id].inp_binSizes.resize(numBins);
        cell_inputs[cell_id].inp_size_dist.resize(net.n_reactions*numBins);
        cell_inputs[cell_id].inp_delSZ.resize(net.n_reactions*numBins);
        std::copy(init_bin_edges.begin(),init_bin_edges.end(),cell_inputs[cell_id].inp_binEdges.begin());
        std::copy(init_size_bins.begin(),init_size_bins.end(),cell_inputs[cell_id].inp_binSizes.begin());
        std::fill(cell_inputs[cell_id].inp_size_dist.begin(),cell_inputs[cell_id].inp_size_dist.end(),0.0);
        std::fill(cell_inputs[cell_id].inp_delSZ.begin(),cell_inputs[cell_id].inp_delSZ.end(),0.0);
    }
    PLOGI << "generated dust size distribution";
}

// load the shock data from the user's specified shock parameters: velocity, temperature, pile up factor
void
nuDust::gen_shock_array_frm_val()
{
    if(isnan(nu_config.shock_velo) || isnan(nu_config.shock_temp) || isnan(nu_config.sim_start_time))
    {
        std ::cout << "! Missing parameters needed to generate the shock arrays for each cell. \n";
        std::cout << "! shock_velo, shock_temp, or shock_time are missing from the config.\n"; 
        exit(1);
    }
    // generate shock velo arrays and a shock temp and time of shock from config
    double shock_velo = nu_config.shock_velo;
    double shock_temp = nu_config.shock_temp;
    double shock_time = nu_config.sim_start_time;
    for ( const auto &ic : cell_inputs)
    {
        auto cell_id = ic.first;
        cell_inputs[cell_id].inp_vd.resize(numBins*net.n_reactions);
        std::fill(cell_inputs[cell_id].inp_vd.begin(), cell_inputs[cell_id].inp_vd.end(), shock_velo);
        cell_inputs[cell_id].inp_shock_temp = shock_temp;
        cell_inputs[cell_id].inp_shock_time = shock_time;
    }
    account_for_pileUp();
    PLOGI << "Set up shock arrays from config params";
}

// load the size ditribution from file
void
nuDust::load_sizeDist()
{
    // This loads just the input size distribution file. It get grain sizes and species included.
    std::ifstream sd_file ( nu_config.sizeDist_file );
    std::string line_buffer;
    std::vector<std::string> line_tokens;

    if ( sd_file.is_open() )
    {   
        std::vector<std::string>  SD_grn_names;
        std::getline ( sd_file, line_buffer );
        boost::split ( line_tokens, line_buffer, boost::is_any_of ( " \t" ), boost::token_compress_on );
        SD_grn_names.assign ( line_tokens.begin(), line_tokens.end() );

        std::getline ( sd_file, line_buffer );
        boost::split ( line_tokens, line_buffer, boost::is_any_of ( " \t" ), boost::token_compress_on );

        for ( auto it = line_tokens.begin() ; it != line_tokens.end(); ++it )
        {
            size_bins_init.push_back ( boost::lexical_cast<double> ( *it ));
        }
 
        numBins = size_bins_init.size();
        double low = std::floor(std::log10(size_bins_init[0]));
        double high = std::floor(std::log10(size_bins_init[numBins-1]))+1;
        auto expDel = (high - low)/static_cast<double>(numBins);
        init_bin_edges.resize(numBins+1);
        std::generate(init_bin_edges.begin(), init_bin_edges.end(), [&, n = 0] () mutable { return std::pow(10.,low+static_cast<double>(n++)*expDel); } );
        
        std::vector<int> grn_idx(net.n_reactions);
        // match grians: network ID with input file ID
        for( auto fn_id = 0; fn_id < SD_grn_names.size(); fn_id++){
            for(auto gn_id =0; gn_id < net.n_reactions; gn_id++){
                if(SD_grn_names[fn_id]==net.reactions[gn_id].prods[0])
                {
                    grn_idx[gn_id]=fn_id;
                }
            }
        }
        while ( std::getline ( sd_file, line_buffer ) )
        {
            boost::split ( line_tokens, line_buffer, boost::is_any_of ( " \t" ), boost::token_compress_on );
            auto cell_id = boost::lexical_cast<int> ( line_tokens[0] );
            auto cell_time = boost::lexical_cast<double> ( line_tokens[1] );
            cell_inputs[cell_id].inp_cell_time = cell_time;
            std::vector<double> input_SD;
            cell_inputs[cell_id].inp_binSizes.resize(numBins);
            std::copy(size_bins_init.begin(),size_bins_init.end(),cell_inputs[cell_id].inp_binSizes.begin());
            cell_inputs[cell_id].inp_binEdges.resize(numBins+1);
            std::copy(init_bin_edges.begin(),init_bin_edges.end(),cell_inputs[cell_id].inp_binEdges.begin());

            for ( auto it = line_tokens.begin() + 2; it != line_tokens.end(); ++it )
            {
                try
                {   
                    input_SD.push_back ( boost::lexical_cast<double> ( *it ));
                }
                catch(const std::exception& e)
                {
                    PLOGE << "bad value in sd file line "<< cell_id;
                    input_SD.push_back ( boost::lexical_cast<double> ( 0.0 ) );
                    exit(1);
                }
            }
            cell_inputs[cell_id].inp_size_dist.resize(net.n_reactions*numBins);
            cell_inputs[cell_id].inp_delSZ.resize(net.n_reactions*numBins);
            std::fill(cell_inputs[cell_id].inp_delSZ.begin(),cell_inputs[cell_id].inp_delSZ.end(),0.0);
            for( auto gid=0; gid<net.n_reactions; gid++)
            {
                for( auto sid=0; sid<numBins; sid++)
                {
                    cell_inputs[cell_id].inp_size_dist[gid*numBins+sid]=input_SD[grn_idx[gid]*numBins+sid];
                }
            }
        }
    }
    else
    {
        PLOGE << "Cannot open size dist file " << nu_config.sizeDist_file;
        exit(1);
    }
}

// load the abundance file
void
nuDust::load_initial_abundances()
{
    std::ifstream abundance_file ( nu_config.abundance_file );
    std::string line_buffer;
    std::vector<std::string> line_tokens;

    if ( abundance_file.is_open() )
    {
        bool missingCO = true;
        bool missingSiO = true;
        std::getline ( abundance_file, line_buffer );
        boost::split ( line_tokens, line_buffer, boost::is_any_of ( " \t" ), boost::token_compress_on );
        initial_elements.assign ( line_tokens.begin() + 1, line_tokens.end() );
        // checking to see if CO is in the input file
        if(std::find(initial_elements.begin(), initial_elements.end(), "CO") != initial_elements.end())
        {
            missingCO = false;
        }
        else
        {
            // only add an extra spot for CO if it isn't in the input abundance file
            initial_elements.emplace_back("CO");
        }
        // checking to see if SiO is in the input file
        if(std::find(initial_elements.begin(), initial_elements.end(), "SiO") != initial_elements.end())
        {
            missingSiO = false;
        }
        else
        {
            // only add an extra spot for SiO if it isn't in the input abundance file
            initial_elements.emplace_back("SiO");
        }

        while ( std::getline ( abundance_file, line_buffer ) )
        {
            boost::split ( line_tokens, line_buffer, boost::is_any_of ( " \t" ), boost::token_compress_on );

            auto cell_id = boost::lexical_cast<int> ( line_tokens[0] );
            // getting indices for premaking CO and SiO
            auto CO_idx  = get_element_index("CO");
            auto C_idx   = get_element_index("C");
            auto O_idx   = get_element_index("O");
            auto SiO_idx = get_element_index("SiO");
            auto Si_idx  = get_element_index("Si");

            for (auto it = line_tokens.begin() + 1; it != line_tokens.end(); ++it) 
            {
                cell_inputs[cell_id].inp_init_abund.push_back(boost::lexical_cast<double>(*it));
            }
            // only add an extra spot for CO if it isn't in the input abundance file
            if(missingCO)
            {
                cell_inputs[cell_id].inp_init_abund.push_back(0.0);
            }
            // only add an extra spot for CO if it isn't in the input abundance file
            if(missingSiO)
            {
                cell_inputs[cell_id].inp_init_abund.push_back(0.0);
            }
            // premake CO and SiO
            premake(C_idx, O_idx, CO_idx, cell_id);
            premake(Si_idx, O_idx, SiO_idx, cell_id);
        }
    }
    else
    {
        PLOGE << "Cannot open abundance file " << nu_config.abundance_file;
        exit(1);
    }
    PLOGI << "loaded abundance file. loaded " << initial_elements.size() << " abundances.";
}

// check if SiO2 or CO are specified in the abundance data, if not add them. Convert all Si, O, and C to SiO2 or CO
void
nuDust::premake(const int s1, const int s2, const int sp, const int cell_id)
{
  auto x1 = cell_inputs[cell_id].inp_init_abund[s1];
  auto x2 = cell_inputs[cell_id].inp_init_abund[s2];

  if (x2 > x1) {
    cell_inputs[cell_id].inp_init_abund[sp] = x1;
    cell_inputs[cell_id].inp_init_abund[s2] = x2 - x1;
    cell_inputs[cell_id].inp_init_abund[s1] = 0.0;
  } else {
    cell_inputs[cell_id].inp_init_abund[sp] = x2;
    cell_inputs[cell_id].inp_init_abund[s1] = x1 - x2;
    cell_inputs[cell_id].inp_init_abund[s2] = 0.0;
  }
}

// adjust the number density of species to account for pile up
void
nuDust::account_for_pileUp()
{
    double pileUpFactor = nu_config.pile_up_factor;
    for ( const auto &ic : cell_inputs)
    {
        auto cell_id = ic.first;
        for (size_t idx = 0; idx < initial_elements.size(); idx++) 
        {
            cell_inputs[cell_id].inp_init_abund[idx] *= pileUpFactor * std::pow(cell_inputs[cell_id].inp_shock_time/cell_inputs[cell_id].inp_cell_time,-3.0);
        }
    }
    PLOGI << cell_inputs[1].inp_init_abund[0] << " after adjust for time";
}

// load the sputter parameters and calcuate additional constants needed for thermal and nonthermal sputtering
void
nuDust::load_sputter_params()
{
    using constants::amu2g;
    using constants::amu2Kg;
    using constants::echarge_sq;
    using constants::bohrANG;

    using numbers::onehalf;
    using numbers::onethird;
    using numbers::twothird;
    using numbers::onefourth;
    using numbers::four;
    using numbers::eight;
    using numbers::one;

    using utilities::square;

    size_t numReact = net.n_nucleation_reactions;
    size_t numGas = initial_elements.size();
    sputARR.alloc_vecs(numReact,numGas);
    for(size_t gsID=0; gsID < numGas; ++gsID)
    {   
        auto mi = sputter.ions.at(initial_elements[gsID]).mi;
        sputARR.mi[gsID] = mi;
        sputARR.zi[gsID] = sputter.ions.at(initial_elements[gsID]).zi;
        sputARR.miGRAMS[gsID] = mi*amu2g;
        sputARR.miKG[gsID] = mi*amu2Kg;
        sputARR.y8_piMi[gsID] = eight / (M_PI * mi*amu2g);
    }
    for ( size_t gidx = 0; gidx < numReact; ++gidx )
    {
        std::string grain_name = net.reactions[gidx].prods[0];

        auto md = sputter.grains.at(grain_name).md;
        auto zd = sputter.grains.at(grain_name).zd;
        auto rhod = sputter.grains.at(grain_name).rhod;

        sputARR.u0[gidx] = sputter.grains.at(grain_name).u0;
        sputARR.md[gidx] = md;
        sputARR.mdGRAMS[gidx] = md*amu2g;
        sputARR.zd[gidx] = zd;
        sputARR.K[gidx] = sputter.grains.at(grain_name).K;
        sputARR.rhod[gidx] = rhod;

        sputARR.msp_2rhod[gidx] = md*amu2g * onehalf / rhod;
        sputARR.three_2Rhod[gidx] = 3.0 / (2.0*rhod);

        for(size_t gsID=0; gsID < numGas; ++gsID)
        {
            auto mi = sputARR.mi[gsID];
            auto mu = md/mi;
            auto invMU = mi/md;
            auto zi = sputARR.zi[gsID];
            sputARR.mu[gidx][gsID] = mu;
            // Biscaro 2016 eq. 8 
            if(mu>one)
                sputARR.alpha[gidx][gsID] = 0.3 * std::pow((mu-0.6),twothird);
            else if (mu <= onehalf)
                sputARR.alpha[gidx][gsID] = 0.2;
            else
                sputARR.alpha[gidx][gsID] = 0.1 * invMU + onefourth * square(mu-onehalf);

            // Biscaro 2016 eq. 3
            if(invMU > 0.3)
                sputARR.Eth[gidx][gsID] = eight * sputARR.u0[gidx] * std::pow(invMU,onethird);
            else
            {
                auto g = four * mi * md / square(mi+md);
                sputARR.Eth[gidx][gsID] = sputARR.u0[gidx] / (g * (1 - g));
            }

            // Biscaro 2016 eq 5
            sputARR.asc[gidx][gsID] = 0.885 * bohrANG * std::pow( (std::pow(zi,twothird) + std::pow(zd,twothird)), - onehalf);

            // up to here matches python

            // Biscaro 2016 eq 4 coefficient
            sputARR.SiCoeff[gidx][gsID] = four * M_PI * sputARR.asc[gidx][gsID] * zi * zd * echarge_sq * mi / (mi + md);
            
            // Biscaro 2016 eq 7 coefficient
            sputARR.eiCoeff[gidx][gsID] = md / (mi + md) * sputARR.asc[gidx][gsID] / (zi * zd * echarge_sq);            
        }
    }
    PLOGI << "calculated sputtering terms";
}

// load the trajectory (environment) data
void
nuDust::load_environment_data()
{
    std::ifstream env_file(nu_config.environment_file);
    std::string line_buffer;
    std::vector<std::string> line_tokens;
    double time;
    if (env_file.is_open()) {
        while (std::getline(env_file, line_buffer)) 
        {
            boost::split(line_tokens, line_buffer, boost::is_any_of(" \t "));
            if (line_tokens.size() == 1) 
            {
                time = boost::lexical_cast<double>(line_tokens[0]);
            } 
            else 
            {
                auto cid   = boost::lexical_cast<int>(line_tokens[0]);
                auto temp  = boost::lexical_cast<double>(line_tokens[1]);
                auto vol   = boost::lexical_cast<double>(line_tokens[2]);
                auto rho   = boost::lexical_cast<double>(line_tokens[3]);
                auto press = boost::lexical_cast<double>(line_tokens[4]);
                auto velo = boost::lexical_cast<double>(line_tokens[5]);
                auto x_cm = boost::lexical_cast<double>(line_tokens[6]);

                cell_inputs[cid].inp_times.push_back(time);
                cell_inputs[cid].inp_temp.push_back(temp);
                cell_inputs[cid].inp_volumes.push_back(vol);
                cell_inputs[cid].inp_rho.push_back(rho);
                cell_inputs[cid].inp_pressure.push_back(press);
                cell_inputs[cid].inp_velo.push_back(velo);
                cell_inputs[cid].inp_x_cm.push_back(x_cm);
            }
        }
    }
    else {
        PLOGE << "Cannont open environment file " << nu_config.environment_file;
    }
    PLOGI << "loaded environment file";
}

// load the shock time, temperature, and velocity from the user specified shock file
void
nuDust::load_shock_params()
{
    std::ifstream sd_file ( nu_config.shock_file );
    std::string line_buffer;
    std::vector<std::string> line_tokens;

    int cell_id  = 0;
    double shock_velo = 0.0;

    if ( sd_file.is_open() )
    {           
        while ( std::getline ( sd_file, line_buffer ) )
        {
            boost::split ( line_tokens, line_buffer, boost::is_any_of ( " \t" ), boost::token_compress_on );
            cell_id = boost::lexical_cast<int> ( line_tokens[0] );
            cell_inputs[cell_id].inp_shock_time = boost::lexical_cast<double> ( line_tokens[1] );
            cell_inputs[cell_id].inp_shock_temp = boost::lexical_cast<double> ( line_tokens[2] );
            shock_velo = boost::lexical_cast<double> ( line_tokens[3] );
            cell_inputs[cell_id].inp_vd.resize(size_bins_init.size()*net.n_reactions);
            std::fill(cell_inputs[cell_id].inp_vd.begin(), cell_inputs[cell_id].inp_vd.end(), shock_velo);
        }
        //account_for_pileUp();
    }
    else
    {
        PLOGE << "Cannot open shock param file " << nu_config.shock_file;
        exit(1);
    }
    PLOGI << "loaded shock parameters file";
}

// not currently used. this calculates where (time, cell ID) a shock is identified in the input data.
void 
nuDust::find_shock()
{
    auto key = cell_inputs.begin()->first;
    std::vector<double> times;
    times.assign(cell_inputs[key].inp_times.begin(),cell_inputs[key].inp_times.end());

    std::vector<size_t> cell_ids;
    for(const auto &ic: cell_inputs)
    {
        cell_ids.push_back(ic.first);
    }
    for(size_t tidx = 1; tidx < times.size()-2; ++tidx)
    {
        for (size_t cid = 2; cid < cell_ids.size()-2; ++cid)
        {
            auto shock_measure = 0.0;
            bool shockDetected = false;
            auto delp1 = cell_inputs[cell_ids[cid+1]].inp_pressure[tidx] - cell_inputs[cell_ids[cid-1]].inp_pressure[tidx];
            shock_measure = std::abs(delp1) / std::min(cell_inputs[cell_ids[cid+1]].inp_pressure[tidx], cell_inputs[cell_ids[cid-1]].inp_pressure[tidx]) - 0.33;
            shock_measure = std::max(0.0, shock_measure);
            if (shock_measure > 0.0)
                shockDetected =  true;
            if (cell_inputs[cell_ids[cid-1]].inp_velo[tidx] < cell_inputs[cell_ids[cid+1]].inp_velo[tidx])
            {
                shockDetected = false;
            }    
            if(shockDetected)
            {
                cell_inputs[cid].inp_shock_times_arr.push_back(times[tidx]);
                cell_inputs[cid].inp_shock_velo_arr.push_back( (cell_inputs[cell_ids[cid]].inp_x_cm[tidx-1]-cell_inputs[cell_ids[cid+1]].inp_x_cm[tidx-1])/(times[tidx]-times[tidx-1]));
                cell_inputs[cid].inp_shock_bool_arr.push_back(true);
            }
        }
    }
}

// generates the solution vector based on the number of grains, size bins, and gas species
void
nuDust::generate_sol_vector()
{
    using constants::N_MOMENTS;
    size_t numReact = net.n_nucleation_reactions;
    size_t numGas = initial_elements.size();
    numBins = nu_config.bin_number;
    std::vector<double> empty_moments(N_MOMENTS*numReact,0.0);

    for ( const auto &ic : cell_inputs)
    {
        int cell_id = ic.first;
        cell_inputs[cell_id].inp_solution_vector.resize(numGas + N_MOMENTS * numReact + numReact * numBins);
        std::copy(cell_inputs[cell_id].inp_init_abund.begin(), cell_inputs[cell_id].inp_init_abund.end(), std::back_inserter(cell_inputs[cell_id].inp_solution_vector));
        std::copy(empty_moments.begin(), empty_moments.end(), std::back_inserter(cell_inputs[cell_id].inp_solution_vector));
        std::copy(cell_inputs[cell_id].inp_size_dist.begin(), cell_inputs[cell_id].inp_size_dist.end(), std::back_inserter(cell_inputs[cell_id].inp_solution_vector));
    }
    PLOGI << "generated solution vector";
}

// define the dump and restart file names. needed in order to check if a restart is needed for the cell.
void 
nuDust::load_outputFL_names()
{
    //creating the prefix to the output and restart files
    // todo::may need to rethink the names
    std::ostringstream binNumObj;

    binNumObj << nu_config.bin_number;
    std::string binNum = binNumObj.str();
    std::string modNum = nu_config.mod_number;

    name = "output/B"+binNum+"_"+net.network_label +"_";   
    nameRS = "restart/restart_B"+binNum+"_"+net.network_label +"_";  
    PLOGI << "data and restart file names defined";
}

// if a restart file exists for the cell, load data from file
void
nuDust::create_restart_cells(int cell_id)
{
    std::ifstream rs_file ( nameRS+std::to_string ( cell_id ) + ".dat" );
    std::string line_buffer;
    std::vector<std::string> line_tokens;

    if ( rs_file.is_open() )
    {
        std::getline ( rs_file, line_buffer );
        boost::split ( line_tokens, line_buffer, boost::is_any_of ( " \t" ), boost::token_compress_on );
        cell_inputs[cell_id].sim_start_time = boost::lexical_cast<double> ( line_tokens[0] );
        std::getline ( rs_file, line_buffer );
        boost::split ( line_tokens, line_buffer, boost::is_any_of ( "  \t" ), boost::token_compress_on );
        for ( auto it = line_tokens.begin(); it != line_tokens.end(); ++it )
        {
            try{cell_inputs[cell_id].inp_vd.push_back ( boost::lexical_cast<double> ( *it ));}
            catch(const std::exception& e){}
        }
        std::getline ( rs_file, line_buffer );
        boost::split ( line_tokens, line_buffer, boost::is_any_of ( "  \t" ), boost::token_compress_on );
        for ( auto it = line_tokens.begin(); it != line_tokens.end(); ++it )
        {
            try{cell_inputs[cell_id].inp_delSZ.push_back ( boost::lexical_cast<double> ( *it ));}
            catch(const std::exception& e){}
        }
        std::getline ( rs_file, line_buffer );
        boost::split ( line_tokens, line_buffer, boost::is_any_of ( "  \t" ), boost::token_compress_on );
        for ( auto it = line_tokens.begin(); it != line_tokens.end(); ++it )
        {
            try{cell_inputs[cell_id].inp_solution_vector.push_back ( boost::lexical_cast<double> ( *it ));}
            catch(const std::exception& e){}
        }
    }
    else
    {
        PLOGE << "Cannot open restart file " << nameRS+std::to_string ( cell_id ) + ".dat";
        exit(1);
    }
}

// creates the simulation cells. this checks if there is an output file or restart file. If there are no restart or output file, create cell. If there's a restart file, load that data instead. If there's an output file and no restart, assume that cell has completed integration.
void
nuDust::create_simulation_cells()
{
    for ( const auto &ic : cell_inputs)
    {
        if (not std::filesystem::exists(name+std::to_string ( ic.first ) + ".dat"))
        {
            if ( std::filesystem::exists(nameRS+std::to_string ( ic.first ) + ".dat"))
            {
                create_restart_cells(ic.first);
            }
            cells.emplace_back ( &net, &sputARR, &nu_config, ic.first, initial_elements, ic.second );
        }
    }
    PLOGI << "Created " << cell_inputs.size() << " cells";
}

// begin calculations
void
nuDust::run()
{  
    PLOGI << "Entering main integration loop";

    if(nu_config.do_destruction==1 && nu_config.do_nucleation==1)
    {
        std::cout << "! Starting Nucleation and Destruction\n";
    }
    else if(nu_config.do_nucleation==1)
    {
        std::cout << "! Starting Nucleation\n";
    }
    else if(nu_config.do_destruction==1)
    {
        std::cout << "! Starting Destruction\n";
    }
    else
    {
        std::cout << "! destruction and nucleation were both not selected. \n";
        std::cout << "! Try again\n";
    }

    // cells.size()
    #pragma omp parallel num_threads(3)
    {
        #pragma omp for nowait
        for (auto i = 0; i < cells.size(); ++i)
        {
            cells[i].solve(); 
        }
    }

    PLOGI << "Leaving main integration loop";
}
