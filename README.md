

# nuDustC
Nucleating Dust Code in C++

nuDustC++ is used to calculate the dust nucleation and destruction in gaseous systems. 

***All Units are in CGS***

# Installation
## Dependencies
**Required:** OpenMP, MPI, Boost, SunDials, Plog.
This build uses cmake. 

### Plog 
Plog is a header-only package. Running cmake fetches Plog, but it can also be acquired by:

```
git clone git@github.com:SergiusTheBest/plog.git
```

### Boost
If using conda, install boost as:

```
conda install conda-forge::boost
```

Further documentation and installation of Boost can be found at:

https://www.boost.org/

### Sundials
On MacOS using HomeBrew, run:

```
brew install sundials
```

Installation instructions for sundials can be found at: 

https://computing.llnl.gov/projects/sundials/faq#inst

## Building nuDustC++
To build nuDustC++, go to the head of the git repository (nuDustC/) and run:

```
$> mkdir build;
$> cd build;
$> cmake .. ;
$> make;
```

To build nuDustC++ in release mode to run faster, instead of the above, run:

```
$> mkdir build;
$> cd build;
$> cmake -DCMAKE_BUILD_TYPE=Release .. ;
$> make;
```

Run nuDustc with

```
$> ./nudustc++ -c data/inputs/test_config.ini
```

You may need to add the paths of the packages to *LD_LIBRARY_PATH*. If using conda nuDustC++ can then be run by:

```
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:~/miniconda3/lib ./nudustc++ -c data/inputs/test_config.ini
```

Or specifying the paths to the packages and libraries in *LD_LIBRARY_PATH*.

OpenMP (NUDUSTC_ENABLE_OPENMP), MPI (NUDUSTC_ENABLE_MPI), and sundials (NUDUSTC_USE_SUNDIALS) are turned off. This can be edited in the CMakeLists.txt file. 

If using MPI, run with

```
$> mpirun -N 'n' ./nudustc++ -c /data/inputs/test_config.txt
```

where 'n' is the number of processors.

### Inputs
Required: Config file. This lists the various input information such as data files, integration parameters, and calculation options.

The input file is specified after the executable such as,

```
$> ./nudustc++ -c /path/to/input/config/file/config_file.txt
```

Descriptions of the input configuration file structure and data files are below.

# Configuration File

### Integration Parameters 
*ode_dt_0*: Initial timestep for integration.

*ode_abs_err*: The Solver's absolute error.

*ode_rel_err*: The Solver's relative error.

*ode_dt_min*: The minimum allowed timestep.

*de_dt_max*: The maximum allowed timestep.

These determine the integrator's timesteps and allowed error. For a quick but less accurate run, increase the 'dt' and lower the '_err' parameters. Conversely, for a more time consuming but accurate run, lower 'dt' and '_err' parameters. The 'dt' is best determined by the timesteps used in the original hydrdynamical run of the trajectory data (described below as the *environment_file*).

    
### Data Files
*sizeDist_file*: This describes the size distribution for the model. Each cell is described in one line. Each line is an array of size distributions of grain species in the order specified in the header line.

*environment_file*: This contains the trajectory data for each timestep. The time is specified on a single line. Below, each cell is described in a single line: cell_ID, temperature (K), volume (cm^3), density(g/cm^3), pressure (Ba), velocity (cm/s), radius (cm).

*network_file*: This includes the chemical network of grain reactions. Each grain species takes up one line in this order: reactants, "->", products, "|", key species, Gibbs free energy 'A' term (A/10^4 K), Gibbs free energy 'B' term, surface energy of the condensate (ergs/cm^2), radius of condensate (angstroms). 

*abundance_file*: This lists the names of gas species in the header. Each cell has one line listing: cell ID and number density for each gas species. 

*shock_file*: This contains information on a shock. Each cell has one line: cell ID, the time of the shock, the shock temperature, the shock velocity.

Examples of each file is provided in the 'data' directory. The example files are prefaced by 'test_'.

### Size Distribution Parameters
*size_dist_min_rad_exponent_cm*: The exponent of the left edge of the distribution.

*size_dist_max_rad_exponent_cm*: The exponent of the right edge of the distribution.

*number_of_size_bins*: The number of size bins.
    
### Control Nucleation and Destruction
*do_destruction*: Set to 1 to enable destruction, 0 to disable.

*do_nucleation*: Set to 1 to enable nucleation, 0 to disable.

If both are set to '1', both destruction and nucleation are calculated. 

### Data Output Controls
*io_dump_n_steps* : Number of cycles until a dump file is updated.  

*io_restart_n_steps*: Number of cycles until a restart file is updated. 

### User Specified Shock Parameters
*pile_up_factor*: This is used to calculate the increase in density when a shock passes through. The density is multiplied by this number. 

*shock_velo*: This is the velocity of the shock that will be applied to all cells. 

*shock_temp*: This is the temperature of the shock that will be applied to all cells.

*sim_start_time*: This is the time of the shock. It is used to calculate the change in volume, assuming homologous expansion.

### Nomenclature
*mod_number*: This specifies the model used in order to differentiate multiple runs.

# Functionality
There are 5 main calculation paths including destruction or nucleation calculations. Destruction has 2 main branches based on if the user is specifying shock values or reading them in from a file. 

Nucleation & Destruction With User Input Shock Temperature & Velocity:
  Required Input Files: Hydrodynamical Trajectory file (Time, Temperature, Volumes, Density, Pressure, Velocity), Abundance file, & Network file, Shock Velocity, Shock Temperature, Shock Time, Pile up factor, Size Distribution File or Size Parameters

Nucleation & Destruction With Shock Times and Velocities from a file:
  Required Input Files: Hydrodynamical Trajectory file (Time, Temperature, Volumes, Density, Pressure, Velocity), Abundance file, & Network file, Shock Velocity, Shock Temperature, Shock Time, Pile up factor, Size Distribution file or Size Parameters

Nucleation:
  Required Input Files: Hydrodynamical Trajectory file (Time, Temperature, Volumes, Density), Abundance file, & Network file, Size Distribution file or Size Parameters

Destruction With Shock Times and Velocities from a file:
  Required Input Files: Shock file (Cell #, Time, Shock Temperature, Shock Velocity), Size Distribution file or Size Parameters, Abundance file, & Network file

Destruction With User Input Shock Temperature & Velocity:
  Required Input Files: Shock Velocity, Shock Temperature, Shock Time, Pile up factor, Size Distribution file or Size Parameters, Abundance file, & Network file

# Output Files
During a run, an 'output' and 'restart' directory are created. Saved in the folders are output data as a function of time and restart files used to restart a run for each individal cell.

The restart file is structured as:

```
time
array of velocities
array of dust size changes per bin
the integrator solution array (abundances of gases, moments of each dust grain, size bins of each dust grain)
```

The output data is structured as:

```
names of grain species
an array corresponding to the grain size of the size bin
the initial integrator solution array (abundances of gases, moments of each dust grain, size bins of each dust grain)
time 1
the integrator solution array at time 1 (abundances of gases, moments of each dust grain, size bins of each dust grain)
time 2
the integrator solution array at time 2 (abundances of gases, moments of each dust grain, size bins of each dust grain)
...
time end
the integrator solution array at time end (abundances of gases, moments of each dust grain, size bins of each dust grain)
```

# Selecting Integrators and Interpolators
The integrator is setup in *src/cell.cpp* in the *solve()* function. nuDustC++ comes defaulted with a Runge–Kutta–Dormand–Prince 5 integrator. Additional information on the available integrators offered by Boost can be found at:

https://www.boost.org/doc/libs/1_78_0/libs/numeric/odeint/doc/html/index.html

nuDustC++ currently uses a Makima 1-D interpolator. The interpolator is defined in *include/cell.h* in the 'cell' class declaration. Additional interpolators offered by Boost can be found at:

https://www.boost.org/doc/libs/1_78_0/libs/math/doc/html/interpolation.html

# Testing
To run a test of nudustc++,

```
./nudustc++ -c data/inputs/test_config.ini
```

nuDustC++ should complete the test run in under a minute or two in release mode. If it doesn't, try changing the configuration file to output data after more cycles by changing *io_dump_n_steps* and *io_restart_n_steps*. This will produce data files in the build directory's "output/" directory and restart data in the "restart/" directory.

# Restarting a Run
nuDustC++ automatically checks for restart files when creating each cell. If a restart file is found, that data is loaded into the cell object. If no restart file is found, the cell is initialized with data from the input files. Make sure the same config file used to start the run is selected when restarting. 

# Common Pitfalls
If the compiler cannot find required packages or libraries, make sure LD_LIBRARY_PATH is up to date and points to the location of each package or library.

MPI issues: make sure you have the installed location of MPI in your path. You might need to change CMakeLists.txt depending on your MPI build.

Ensure the configuration file points to the accessible location of each input file and contains the necessary parameters for the calculation path. 

Boost Math errors: Some older versions of Boost might not work with the Gauss Kronrod Quadrature. We haven't checked any builds older than 1.78.0. If you get math errors, try updating Boost. 

# Community guidelines

If you use nuDustC++ and need help, submit an issue to the nuDustC++ repository. If you'd like to contribute, just fork and submit a pull request. One of the developers will review your PR.

# Copyright
BSD 3-Clause License

Copyright (c) 2023, Los Alamos National Laboratory

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.



