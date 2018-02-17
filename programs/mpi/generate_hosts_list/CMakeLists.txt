cmake_minimum_required(VERSION 2.7)
project(BIDON_MPI)

find_package(MPI REQUIRED)

set (CMAKE_CXX_STANDARD 11)
include_directories(${MPI_INCLUDE_PATH})

add_executable(generate_hosts_list generate_hosts_list.cpp)
target_link_libraries(generate_hosts_list ${MPI_LIBRARIES})

if(MPI_COMPILE_FLAGS)
  set_target_properties(generate_hosts_list PROPERTIES
    COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
  set_target_properties(generate_hosts_list PROPERTIES
    LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()
