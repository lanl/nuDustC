/*© 2023. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for
Los Alamos National Laboratory (LANL), which is operated by Triad National
Security, LLC for the U.S. Department of Energy/National Nuclear Security
Administration. All rights in the program are. reserved by Triad National
Security, LLC, and the U.S. Department of Energy/National Nuclear Security
Administration. The Government is granted for itself and others acting on its
behalf a nonexclusive, paid-up, irrevocable worldwide license in this material
to reproduce, prepare. derivative works, distribute copies to the public,
perform publicly and display publicly, and to permit. others to do so.*/

#include "configuration.h"
#include "logging.h"
#include "nudust.h"
#include "timer.h"

#include <plog/Initializers/RollingFileInitializer.h>
#include <plog/Log.h>
#include <string>
#include <vector>
//#ifdef NUDUSTC_ENABLE_MPI
#include <mpi.h>
//#endif

//#ifdef NUDUSTC_ENABLE_OPENMP
//#include <omp.h>
//#endif

// TODO: (maybe) this should be it's own module...not a lot done tho
#include <boost/program_options.hpp>
#include <boost/program_options/options_description.hpp>

// TODO: move these elsewhere
inline void banner() {
  std::cout << R"(
       +                                :
 .              .      `    |   o       !  '
   |    o   `       .      -O-         ,|.'
  -O-    .        _____     |   . -_--(-O-)----
   |  .       +  |  __ \     .    | |  `|'  
      ____  _   _| |  | |_   _ ___| |_' | 
     |  _ \| | | | |  | | | | / __| __| ! 
     | | | | |_| | |__| | |_| \__ \ |_ 
     |_| |_|\__,_|_____/ \__,_|___/\__|
===============================================
Computational modelling of astrophysical dust
and chemistry.

===============================================
| Developed by:
--------------
|     Sarah Stangl (sarahstangl@lanl.gov)
|     Christopher Mauney (mauneyc@lanl.gov)
===============================================     
)";
}

int main(int argc, char *argv[]) {
  int rank = 0, size = 1;
//#ifdef NUDUSTC_ENABLE_MPI
  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
//#endif
  namespace po = boost::program_options;

  // setup command-line options
  std::string config_filename;
  std::string log_filename;

  po::options_description desc("nuDust options");
  desc.add_options()("help", "print help message")(
      "config_file,c", po::value<std::string>(&config_filename),
      "filename with runtime parameters")(
      "log_file,l",
      po::value<std::string>(&log_filename)->default_value("log.txt"),
      "filename of log");

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);

  if (vm.count("help")) {
    std ::cout << desc << "\n";
    return 1;
  }

  if (!vm.count("config_file")) {
    std ::cout << "missing required configuration file!\n";
    std ::cout << "\tnudust++ -c data/inputs/default_config.ini\n";
    std ::cout << desc << "\n";
    return 1;
  }
  banner();

  std ::cout << "\n";
  std ::cout << "! log file = " << log_filename << "\n";
  std ::cout << "! configuration file = " << config_filename << "\n";
  

  if (rank == 0) {
    //banner();

    //std ::cout << "\n";
    //std ::cout << "! log file = " << log_filename << "\n";
    //std ::cout << "! configuration file = " << config_filename << "\n";
    std ::cout << "! pe = " << size << "\n";
  }
  std ::cout << "! starting ...\n";
#ifdef NUDUSTC_ENABLE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  plog::init(plog::debug, log_filename.c_str());
  plog::init<DetailLog_Root>(plog::info, &rootAppender);
  plog::init<DetailLog_Rank>(plog::info, &rankAppender);

  //SCLOG_0(rank) << "nuDust has started";
  //SCLOG_A(rank) << "starting ...\n";
  std ::cout << "! rank: " << rank << "\n";
  std ::cout << "! size: " << size << "\n";
  std::cout << "! nuDust has started\n";

  PLOGI << "nuDust has started";

  nuDust nd(config_filename, size, rank);
  std::cout << "! Setup has completed. Starting runs...\n";
  nd.run();
  std::cout << "! nuDust has finished!\n";

#ifdef ENABLE_BENCHMARK
  for (const auto &[k, v] : benchmark_map) {
    auto call_avg = std::accumulate(v.begin(), v.end(), 0) / v.size();
    SCLOG_A(rank) << "fn = " << k << " is called " << v.size()
                  << " times: avg/call = " << call_avg << " us\n";
  }
#endif

#ifdef NUDUSTC_ENABLE_MPI
  MPI_Finalize();
#endif
  return 0;
}