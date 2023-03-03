// Coupled-mode Shallow Water Model (CSW) 
//
// Copyright (c) 2016, Samuel M Kelly (smkelly@d.umn.edu)
//
// Publications and presentations that use this code should cite:
// S. M. Kelly, P. F. J. Lermusiaux, T. Duda, and P. J. Haley, Jr. (2016) A Coupled-mode Shallow Water model for tidal analysis: Internal-tide reflection and refraction by the Gulf Stream, J. Phys. Oceanogr., 3661-3679. 

#include <mpi.h>
#include "csw.h"

////////////////////////////////////////////////////////////////////////
// Main program
int main(int argc, char *argv[])
{

	int rank, nrank, s;
	double t=0; // time

	int sW=0; // Write index
	int sD=1; // Diagnostic index (start writing after one period)
	int sF=0; // Index for wind forcing
	int Na=1; // Number of points for diagnostic average (must define)

	////////////////////////////////////////////////////////////////
	// Start the MPI enviornment 
	MPI_Init(&argc,&argv);    
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&nrank);
	
	////////////////////////////////////////////////////////////////////
	// Read grid
	read_grid(rank);

	// Read forcing data
	#ifdef TIDE_FORCING
		read_tides(rank);
	#endif

	//////////////////////////////////////////////////////////////
	// Begin forward integration
	for(s=1; s<(NT+1); s++) {

		// Read WIND data
		#ifdef WIND_FORCING
			if (t >= (sF*DT_F)) {
				read_wind(sF,rank);
				++sF;
			}
		#endif

		////////////////////////////////////////////////////////////////
		// Momentum calls
		calc_forces(Na); 

		// Update momentum 
		timestep_uv();

		// Trade U & V at boundaries
		pass_uv(rank);

		////////////////////////////////////////////////////////////////
		// Write output
		#if defined(WRITE_ETA) || defined(WRITE_VELOCITY) || defined(WRITE_WIND)
			if (t >= (sW*DT_W)) {
				write_output(sW,t,rank);
				++sW;
			}
		#endif
		
		// Write diagnostics
		#if defined(ENERGY) || defined(FLUX) || defined(WORK) || defined(WRITE_SSH) || defined(WRITE_TRANSPORT)
			if (t >= (sD*DT_D)) {
				write_diagnostics(sD,Na,rank);
				++sD;
				Na=0; // Reset averaging counter
			}
			++Na; // Increase the averaging counter
		#endif
		
		////////////////////////////////////////////////////////////////
		// Divergence calls
		calc_divergence();

		// Add internal-tide generating force
		#ifdef TIDE_FORCING
			calc_ITGF(t);
		#endif

		// Update pressure 
		timestep_p();

		// Trade p at boundaries 
		pass_p(rank);

		////////////////////////////////////////////////////////////////
		// Update time
		t=s*DT; // U, V, and (p+p1)/2 now correspond to this time

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

