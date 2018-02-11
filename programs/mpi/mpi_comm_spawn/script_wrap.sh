#!/bin/bash
trap 'echo "EXIT"' EXIT

../../fake_mpi_program/build/fake_mpi_program $*
