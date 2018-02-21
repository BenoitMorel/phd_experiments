cmake_minimum_required(VERSION 2.7)
project(BIDON_MPI)

find_package(MPI REQUIRED)

set (CMAKE_CXX_STANDARD 11)
include_directories(${MPI_INCLUDE_PATH})

add_executable(print_proc_name print_proc_name.cpp)
target_link_libraries(print_proc_name ${MPI_LIBRARIES})

if(MPI_COMPILE_FLAGS)
  set_target_properties(print_proc_name PROPERTIES
    COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
  set_target_properties(print_proc_name PROPERTIES
    LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()
