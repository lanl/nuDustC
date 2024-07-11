# which cmake versions
cmake_minimum_required(VERSION 3.11...3.23)

# project name and properties
project(
  nuDust++
  VERSION 0.1
  DESCRIPTION "dust nucleation and kinetics"
  LANGUAGES CXX)

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
  message(
    FATAL_ERROR
      "You are trying to configure a build in the source directory. Move to another directory and rerun cmake.\nTIP: It is common to do (from the top of the source dir):\n\t \" $> mkdir build && cd build; cmake ..\" "
  )
endif()

option(NUDUSTC_ENABLE_OPENMP ON "Use OpenMP for cell/particle parallelization")
option(NUDUSTC_ENABLE_MPI ON "Use MPI for cell/particle parallelization")
option(NUDUSTC_USE_SUNDIALS OFF "Use sundials CVODE integrator")

# dependencies
list(APPEND BOOST_COMPONENTS program_options filesystem serialization)
find_package(Boost REQUIRED ${BOOST_COMPONENTS})

if(NUDUSTC_USE_SUNDIALS)
  find_package(SUNDIALS 5.8 REQUIRED)
endif()

if(NUDUSTC_ENABLE_OPENMP)
  find_package(OpenMP REQUIRED)
endif()

include(FetchContent)
FetchContent_Declare(
  plog
  GIT_REPOSITORY https://github.com/SergiusTheBest/plog.git
  GIT_TAG 1.1.9)

FetchContent_MakeAvailable(plog)

option(NUDUST_BENCHMARK "Enable benchmark wrappers" OFF)

set(NUD_EXE "nudustc++")

set(NUD_SRCS
    src/cell.cpp
    src/cellobserver.cpp
    src/configuration.cpp
    src/main.cpp
    src/network.cpp
    src/nudust.cpp
    src/reaction.cpp)

set(NUD_HEADERS
    include/axis.h
    include/cell.h
    include/cellobserver.h
    include/configuration.h
    include/constants.h
    include/elements.h
    include/makima.h
    include/network.h
    include/nudust.h
    include/reaction.h
    include/sput_params.h
    include/sputter.h
    include/utilities.h)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_executable(${NUD_EXE} ${NUD_SRCS} ${NUD_HEADERS})

target_include_directories(${NUD_EXE} PRIVATE ${PROJECT_SOURCE_DIR}/include)

set(with_mpi "$<BOOL:${NUDUSTC_ENABLE_MPI}>")
set(with_openmp "$<BOOL:${NUDUSTC_ENABLE_OPENMP}>")
set(use_sundials "$<BOOL:${NUDUSTC_USE_SUNDIALS}>")


if(NUDUSTC_ENABLE_MPI)
  find_package(MPI REQUIRED CXX)
  include_directories(${MPI_CXX_INCLUDE_DIRS})
  target_link_libraries(${NUD_EXE} PRIVATE ${MPI_CXX_LIBRARIES})
  target_link_libraries(${NUD_EXE} PUBLIC MPI::MPI_CXX)
  message("MPI Found")
endif()

target_link_libraries(
  ${NUD_EXE}
  PRIVATE Boost::headers
          Boost::program_options
          Boost::serialization
          Boost::filesystem
          plog::plog
          $<${use_sundials}:"SUNDIALS::cvode;SUNDIALS::nvecserial">
          #$<${with_mpi}:"MPI::MPI_CXX">
          $<${with_openmp}:"OpenMP::OpenMP_CXX">)

target_compile_definitions(
  ${NUD_EXE}
  PRIVATE $<$<BOOL:${NUDUST_BENCHMARK}>:ENABLE_BENCHMARK>
          $<$<BOOL:${NUDUST_ENABLE_MPI}>:NUDUSTC_ENABLE_MPI>
          # latest boost fails with gcc@12
          # https://github.com/boostorg/phoenix/issues/111
          BOOST_PHOENIX_STL_TUPLE_H_)

add_custom_command(
  TARGET ${NUD_EXE}
  POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_SOURCE_DIR}/data/
          $<TARGET_FILE_DIR:${NUD_EXE}>/data)

add_custom_command(
  TARGET ${NUD_EXE}
  POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E make_directory
          $<TARGET_FILE_DIR:${NUD_EXE}>/output)

add_custom_command(
  TARGET ${NUD_EXE}
  POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E make_directory
          $<TARGET_FILE_DIR:${NUD_EXE}>/restart)

# IDEs should put the headers in a nice place
source_group(
  TREE "${PROJECT_SOURCE_DIR}/include"
  PREFIX "Header Files"
  FILES ${NUD_HEADERS})
