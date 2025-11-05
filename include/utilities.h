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

#ifdef NUDUSTC_ENABLE_MPI
#include <mpi.h>
#endif

#ifdef NUDUSTC_ENABLE_OPENMP
#include <omp.h>
#endif

namespace utilities
{

template<class T>
inline constexpr auto const square(const T& value){
  return value * value;
}

template<class T>
inline constexpr auto const cube(const T& value){
  return value * value * value;
}

inline auto const my_rank()
{
  int rank = 0;
#ifdef NUDUSTC_ENABLE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  return rank;
}

inline auto const my_tid()
{
  int tid=0;
#ifdef NUDUSTC_ENABLE_OPENMP
  tid = omp_get_thread_num();
#endif
  return tid;
}


} // namespace utilities
