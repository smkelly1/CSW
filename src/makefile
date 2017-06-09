#
# 'make depend' uses makedepend to automatically generate dependencies 
#               (dependencies are added to end of Makefile)
# 'make'        build executable file 'mycc'
# 'make clean'  removes all .o and executable files
#

SHELL=/bin/bash

CC=mpicc
CFLAGS=-Wall -g -O3 -mcmodel=large

#CC=mpiicc
#CFLAGS=-Wall -O3 -xAVX -mkl -mcmodel=large


INCLUDES = -I./ 
#LFLAGS=-g -pg 
LIBS = -lm -lnetcdf -lmpi

SRCS = main.c calc_divergence.c calc_forces.c calc_ITGF.c init_output.c pass_p.c pass_uv.c read_grid.c read_input.c timestep_p.c  timestep_uv.c write_diagnostics.c write_output.c

OBJS = $(SRCS:.c=.o)

MAIN = cswexec

# To run, try: mpirun --mca orte_base_help_aggregate 0 -np 6 cswexec

# The following part of the makefile is generic; it can be used to 
# build any executable just by changing the definitions above and by
# deleting dependencies appended to the file from 'make depend'

.PHONY: depend clean

all:    $(MAIN)

$(MAIN): $(OBJS) 
	$(CC) $(CFLAGS) $(INCLUDES) -o $(MAIN) $(OBJS) $(LFLAGS) $(LIBS)

# this is a suffix replacement rule for building .o's from .c's
# it uses automatic variables $<: the name of the prerequisite of
# the rule(a .c file) and $@: the name of the target of the rule (a .o file) 
# (see the gnu make manual section about automatic variables)
.c.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c $<  -o $@

clean:
	$(RM) *.o *~ $(MAIN)

depend: $(SRCS)
	makedepend $(INCLUDES) $^

# DO NOT DELETE THIS LINE -- make depend needs it


