cmake_minimum_required(VERSION 2.7)
project(BIDON_MPI)

find_package(MPI REQUIRED)

set (CMAKE_CXX_STANDARD 11)
include_directories(${MPI_INCLUDE_PATH})

add_executable(mpi_and_exec mpi_and_exec.cpp)
target_link_libraries(mpi_and_exec ${MPI_LIBRARIES})

if(MPI_COMPILE_FLAGS)
  set_target_properties(mpi_and_exec PROPERTIES
    COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
  set_target_properties(mpi_and_exec PROPERTIES
    LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()
