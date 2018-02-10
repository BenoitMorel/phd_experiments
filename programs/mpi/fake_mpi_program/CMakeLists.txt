cmake_minimum_required(VERSION 2.7)
project(BIDON_MPI)

find_package(MPI REQUIRED)

set (CMAKE_CXX_STANDARD 11)
include_directories(${MPI_INCLUDE_PATH})

add_executable(fake_mpi_program fake_mpi_program.cpp)
target_link_libraries(fake_mpi_program ${MPI_LIBRARIES})

if(MPI_COMPILE_FLAGS)
  set_target_properties(fake_mpi_program PROPERTIES
    COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
  set_target_properties(fake_mpi_program PROPERTIES
    LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()
