#!/bin/bash 

#-------Copy configuration files------------
cp ../*_in.nc ./
cp ../make_input.m ./make_input.m
cp ../src_v1/csw.h ./csw.h

#-------Start process-----------------------
#time -o log mpirun -np 6 ../src_diag/cswexec > log &
{ time mpirun -np 4 ../src_v1/cswexec; } >> log 2>&1 &
# Run this job with "sh ../run_CSW.sh"


# Debugging options
# compile with -g in cflags and lflags
# mpirun -np 6 valgrind --leak-check=full --suppressions=/usr/share/openmpi/openmpi-valgrind.supp ../src_diag/cswexec
# mpirun -np 6 xterm -e gdb ../src_diag/cswexec
#
# compile with -pg in cflags and lflags
# run the program: mpirun -np 6 ../src_diag/cswexec
# view the report: gprof ../src_diag/cswexec
#
# compile with -convert in cflags and -lgcov  in lflags
# run program
# run gcov main.c to generate main.c.gcov
# examine number of calls: less main.c.gcov



