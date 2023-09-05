# which cmake versions
cmake_minimum_required(VERSION 3.11...3.23)

# project name and properties
project(
  nuDust++
  VERSION 0.1
  DESCRIPTION "dust nucleation and kinetics"
  LANGUAGES CXX
)

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
  message(FATAL_ERROR "You are trying to configure a build in the source directory. Move to another directory and rerun cmake.\nTIP: It is common to do (from the top of the source dir):\n\t \" $> mkdir build && cd build; cmake ..\" ")
endif()

# dependencies
list(APPEND BOOST_COMPONENTS
  program_options
  filesystem
  serialization
)
find_package(Boost REQUIRED ${BOOST_COMPONENTS})

find_package(OpenMP REQUIRED)

find_package(MPI REQUIRED CXX)

include(FetchContent)
FetchContent_Declare(
  plog
  GIT_REPOSITORY https://github.com/SergiusTheBest/plog.git
  GIT_TAG 1.1.9)

FetchContent_MakeAvailable(plog)

option(NUDUST_BENCHMARK "Enable benchmark wrappers" OFF)

set(NUD_EXE "nudust++")

set(NUD_SRCS
  src/cell.cpp
  src/cellobserver.cpp
  src/configuration.cpp
  src/destruction.cpp
  src/main.cpp
  src/network.cpp
  src/nucleation.cpp
  src/nudust.cpp
  src/reaction.cpp
)

set(NUD_HEADERS
  include/axis.h
  include/cell.h
  include/cellobserver.h
  include/configuration.h
  include/constants.h
  include/destruction.h
  include/elements.h
  include/makima.h
  include/network.h
  include/nucleation.h
  include/nudust.h
  include/reaction.h
  include/sput_params.h
  include/sputter.h
  include/utilities.h
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_executable(${NUD_EXE} ${NUD_SRCS} ${NUD_HEADERS})

target_include_directories(${NUD_EXE}
  PRIVATE
    ${PROJECT_SOURCE_DIR}/include
)

target_link_libraries(${NUD_EXE}
  PRIVATE
    Boost::headers
    Boost::program_options
    Boost::serialization
    Boost::filesystem
    plog::plog
    MPI::MPI_CXX
    OpenMP::OpenMP_CXX
)

target_compile_definitions(${NUD_EXE}
  PRIVATE
    $<$<BOOL:${NUDUST_BENCHMARK}>:ENABLE_BENCHMARK>
  # latest boost fails with gcc@12 
  # https://github.com/boostorg/phoenix/issues/111
    BOOST_PHOENIX_STL_TUPLE_H_
)

add_custom_command(
  TARGET ${NUD_EXE} POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy_directory 
    ${CMAKE_SOURCE_DIR}/data/ $<TARGET_FILE_DIR:${NUD_EXE}>/data
)

add_custom_command(
  TARGET ${NUD_EXE} POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E make_directory
    $<TARGET_FILE_DIR:${NUD_EXE}>/output
)

add_custom_command(
  TARGET ${NUD_EXE} POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E make_directory
    $<TARGET_FILE_DIR:${NUD_EXE}>/restart
)


# IDEs should put the headers in a nice place
source_group(
  TREE "${PROJECT_SOURCE_DIR}/include"
  PREFIX "Header Files"
  FILES ${NUD_HEADERS}
)
