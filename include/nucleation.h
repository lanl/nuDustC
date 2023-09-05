/*© 2023. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos
National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S.
Department of Energy/National Nuclear Security Administration. All rights in the program are.
reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear
Security Administration. The Government is granted for itself and others acting on its behalf a
nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare.
derivative works, distribute copies to the public, perform publicly and display publicly, and to permit.
others to do so.*/

#pragma once

#include "elements.h"
#include "network.h"
#include "sput_params.h"
#include "utilities.h"
#include "constants.h"
#include "cell.h"

#include <vector>
#include <string>
#include <map>
#include <chrono>
#include <fstream>
#include <boost/format.hpp>
#include <boost/math/interpolators/makima.hpp>
#include <boost/numeric/odeint.hpp>
#include <boost/serialization/array_wrapper.hpp>

class nucleation
{
    cell_state      cell_nuc;
    network*        net;
    params*         sputARR;
    xkin::element_list_t elm;

public:
    nucleation(cell_state cell_n,network* n);
    virtual ~nucleation(){};

    void nucleate(const std::vector<double>& x, const double time);
    void add_new_grn(const std::vector<double>& x);
};