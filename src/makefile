# Simple makefile for CSW
#
########################################################################
# Set compiler options

# This set is for openmpi with gcc 
CC = mpicc
CFLAGS = -Wall -O3 -mcmodel=large -fpic -fcommon

# NOTE: -01 seems to work find with gcc, but -O3 leads to problems in pass_UV.c that produces instabilities at the corners of the tiles. It could be a precision issue, or that the optimization is chopping the loops so that uninitialized values are being received. At any rate, putting print statements in the code to debug is a bad idea, because they may change the optimization methods. Some solutions might be (i) compiling pass_UV.c with O1 optimization, (ii) adding a precision flag or loop breaking flag, or (3) rewriting pass_UV.c to be "optimization safe" although, I'm not sure what that would entail. 

# This set is for openmpi with icc
#CC=mpiicc
#CFLAGS=-cc=icc -Wall -O3  -mcmodel=large -mkl -ipo -fpic -axAVX 

LIBS = -lm -lnetcdf -lmpi

# Flags are for profiling/debugging
# LFLAGS=-g -pg 

########################################################################
# Source files and executable name
OBJECTS = main.o calc_divergence.o calc_forces.o calc_ITGF.o pass_p.o pass_uv.o read_grid.o read_tides.o read_wind.o timestep_p.o  timestep_uv.o write_diagnostics.o write_output.o

TARGET=cswexec
# To run, try: mpirun --mca orte_base_help_aggregate 0 -np 6 cswexec

########################################################################
# Make compiles
$(TARGET) : $(OBJECTS) 
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJECTS) $(LFLAGS) $(LIBS)

.PHONY: clean
clean:
	rm -f *.o *.d $(TARGET)



