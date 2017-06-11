// Coupled-mode Shallow Water Model (CSW) 
//
// Release notes, 6-JUN-16:
// Built from the Matlab code CSW_v9.m
// This version does not include mean flow advection
// Written in C using openmpi
// Originally compiled with gcc on Ubuntu Linux 12.04
//
// Release notes, 5-MAY-17:
// Minor changes, removed some diagnostic lines
// Adjusted the passing functions so that E/W periodic boundaries can be turned off
// Need to add online calculations of F, div(F), and C to write-output
// 
// Copyright (c) 2016, Samuel M Kelly (smkelly@d.umn.edu)
//
// Publications and presentations that use this code should cite:
// S. M. Kelly, P. F. J. Lermusiaux, T. Duda, and P. J. Haley, Jr. (2015) A Coupled-mode Shallow Water model for tidal analysis: Internal-tide reflection and refraction by the Gulf Stream, J. Phys. Oceanogr., submitted. 
//

#include <mpi.h>
#include <csw.h>

////////////////////////////////////////////////////////////////////////
// Main program
int main(int argc, char *argv[])
{
	
	int rank, nrank, s;
	double t=0; // time
	int sW=0; // Write index
	
	#ifdef DIAGNOSTICS
		int sD=1; // Diagnostic index (start after one cycle)
		int Na=1; // Number of points for diagnostic average
	#endif
	
	////////////////////////////////////////////////////////////////
	// Start the MPI enviornment 
	MPI_Init(&argc,&argv);    
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&nrank);
	

	////////////////////////////////////////////////////////////////////		
	// Read input data
	read_input(rank);
	
	// Read grid
	read_grid(rank);

	// Initiate output file
	init_output(rank);
		

	//////////////////////////////////////////////////////////////
	// Begin forward integration
	for(s=1; s<(NT+1); s++) {
		
		////////////////////////////////////////////////////////////////
		// Momentum calls
		calc_forces(); 

		// Update momentum 
		timestep_uv();

		// Trade u & v at boundaries (only needed to compute diffusion and Coriolis)
		pass_uv(rank);

		
		////////////////////////////////////////////////////////////////
		// Divergence calls
 		calc_divergence();

		// Add internal-tide generating force
		#ifdef IT_FORCING
			calc_ITGF(t);
		#endif
				
		// Update pressure 
		timestep_p();

		// Trade p1 at boundaries (needed to compute pressure gradients and topographic coupling)
		pass_p(rank);


		////////////////////////////////////////////////////////////////
		// Update time and Write output
		
		t=s*DT; // U, V, and p1 now correspond to this time
		
		// Write ouput
		if (t >= (sW*DT_W)) {			
			write_output(sW,t,rank);
			++sW;
		}
				
		// Write diagnostics
		#ifdef DIAGNOSTICS
			if (t >= (sD*DT_D)) {			
				write_diagnostics(sD,Na,rank);
				++sD;
				Na=0; // Reset averaging counter
			}
			++Na;// Increase the averaging counter
		#endif
				
		// Print progress
		if (rank == 0) {
			printf ("Step %d \n",s);
		}
		
	} // End time-integration loop

		
	/////////////////////////////////////////////////////////////////////
	// Finalize MPI
	MPI_Finalize();
 
return 0;
}

