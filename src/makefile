# Simple makefile for CSW
#
########################################################################
# Set compiler options

# This set is for openmpi with gcc 
CC = mpicc
CFLAGS = -Wall -O3 -mcmodel=large -fpic -fcommon 

# For profiling with gprof and gcov
#CFLAGS = -Wall -O3 -mcmodel=large -fpic -fcommon -pg
#CFLAGS = -Wall -O3 -mcmodel=large -fpic -fcommon -fprofile-arcs -ftest-coverage

# This set is for openmpi with icc
#CC=mpiicc
#CFLAGS=-cc=icc -Wall -O3  -mcmodel=large -mkl -ipo -fpic -axAVX 

LIBS = -lm -lnetcdf -lmpi

# Flags are for profiling/debugging gprof and gcov
#LFLAGS=-g -pg # add -pg to CFLAGS too  
#LFAGS=-fprofile-arcs 

########################################################################
# Source files and executable name
OBJECTS = main.o calc_divergence.o calc_forces.o pass_p.o pass_uv.o read_grid.o read_tides.o read_wind.o write_diagnostics.o write_output.o

TARGET=cswexec
# To run, try: mpirun --mca orte_base_help_aggregate 0 -np 6 cswexec

########################################################################
# Make compiles
$(TARGET) : $(OBJECTS) 
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJECTS) $(LFLAGS) $(LIBS)

.PHONY: clean
clean:
	rm -f *.o *.d $(TARGET)



