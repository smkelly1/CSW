#include <math.h>
#include <mpi.h>
#include "csw.h"

void pass_p(int rank)
{

	int i, j, m;
	int rank_e, rank_w, rank_n, rank_s, row;

	// NOTE: Both sent and received data is labeled by it's direction of travel 
			
	////////////////////////////////////////////////////////////////////
	// Step 0: implement periodic boundaries if we don't need MPI
	
	#if NPX==1 // Simple case for NPX==1
		for(m=0; m<NM; m++) {
			for (j=0; j<NY; j++) {
				p[m][j+1][NX+1]=p[m][j+1][1];   //Second point -> last point
				p[m][j+1][0]=p[m][j+1][NX];     //Second to last point -> first point
			}
		}
	#endif
	
	
	////////////////////////////////////////////////////////////////////
	// Step 1: Identify ranks
	
	// Identify east/west ranks
	#if NPX>1 // Only bother if there are multiple tiles in the X direction
		rank_e=rank+1;
		rank_w=rank-1;

		// Adjust for periodic boundaries
		if ((rank+1) % NPX == 0) {
			rank_e=(rank-NPX)+1;			
		}

		if (rank % NPX == 0) {
			rank_w=(rank+NPX)-1;			
		}
	#endif

	// Identify north/south ranks	
	#if NPY>1 // Only bother if there are multiple tiles in the Y direction
		rank_n=rank+NPX;
		rank_s=rank-NPX;
		row=rank/NPX;	
	#endif
	
	
	////////////////////////////////////////////////////////////////////
	// Step 2: Load-up the outgoing data
	
	#if NPX>1 // Only bother if there are multiple tiles in the X direction
		// Load points going east and west
		for (m=0; m<NM; m++) {
			for (j=0; j<NY; j++) {
				BOUNDARY.p_Se[m*NY+j]=p[m][j+1][NX];  // Second to last point goes east (Size is NX+1, last point is NX) 
				BOUNDARY.p_Sw[m*NY+j]=p[m][j+1][1];   // Second point goes west
			}
		}		
	#endif

	#if NPY>1 // Only bother if there are multiple tiles in the Y direction
		// Load points going north and south
		for (m=0; m<NM; m++) {
			for (i=0; i<NX; i++) {
				BOUNDARY.p_Sn[m*NX+i]=p[m][NY][i+1];  // Second to last point goes north (Size is NY+2, last point is NY+1) 
				BOUNDARY.p_Ss[m*NX+i]=p[m][1][i+1];   // Second point goes south
			}
		}		
	#endif	
		
	
	////////////////////////////////////////////////////////////////////
	// Step 3: Mail the envelopes
	
	#if NPX>1 // East-west messages: Only bother if there are multiple tiles in the X direction
	
		// First send data east, Odd ranks send first to avoid deadlock
		if (rank % 2) {
			if (rank_e > -1) {
				MPI_Send(BOUNDARY.p_Se,NM*NY,MPI_DOUBLE,rank_e,1,MPI_COMM_WORLD);
			}
			if (rank_w > -1) {
				MPI_Recv(BOUNDARY.p_Re,NM*NY,MPI_DOUBLE,rank_w,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			}
		}
		else {
			if (rank_w > -1) {
				MPI_Recv(BOUNDARY.p_Re,NM*NY,MPI_DOUBLE,rank_w,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			}
			if (rank_e > -1) {
				MPI_Send(BOUNDARY.p_Se,NM*NY,MPI_DOUBLE,rank_e,1,MPI_COMM_WORLD);
			}
		}
		
		// Second send data west, Odd ranks send first to avoid deadlock
		if (rank % 2) {
			if (rank_w > -1) {
				MPI_Send(BOUNDARY.p_Sw,NM*NY,MPI_DOUBLE,rank_w,2,MPI_COMM_WORLD);
			}
			if (rank_e > -1) {
				MPI_Recv(BOUNDARY.p_Rw,NM*NY,MPI_DOUBLE,rank_e,2,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			}
		}
		else {
			if (rank_e > -1) {
				MPI_Recv(BOUNDARY.p_Rw,NM*NY,MPI_DOUBLE,rank_e,2,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			}
			if (rank_w > -1) {
				MPI_Send(BOUNDARY.p_Sw,NM*NY,MPI_DOUBLE,rank_w,2,MPI_COMM_WORLD);
			}
		}
		
	#endif // Done with East-west messages

	#if NPY>1 // North-south messages: Only bother if there are multiple tiles in the Y direction

		// First send data north, Odd ranks send first to avoid deadlock
		if (row % 2) {
			if (rank_n < NPX*NPY) {
				MPI_Send(BOUNDARY.p_Sn,NM*NX,MPI_DOUBLE,rank_n,3,MPI_COMM_WORLD);
			}
			MPI_Recv(BOUNDARY.p_Rn,NM*NX,MPI_DOUBLE,rank_s,3,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		}
		else {
			if (rank_s >= 0){
				MPI_Recv(BOUNDARY.p_Rn,NM*NX,MPI_DOUBLE,rank_s,3,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			}
			if (rank_n < NPX*NPY) {
				MPI_Send(BOUNDARY.p_Sn,NM*NX,MPI_DOUBLE,rank_n,3,MPI_COMM_WORLD);
			}
		}
						
		// Second send data south, Odd ranks send first to avoid deadlock
		if (row % 2) {
			MPI_Send(BOUNDARY.p_Ss,NM*NX,MPI_DOUBLE,rank_s,4,MPI_COMM_WORLD);
			if (rank_n < NPX*NPY) {
				MPI_Recv(BOUNDARY.p_Rs,NM*NX,MPI_DOUBLE,rank_n,4,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			}
		}
		else {
			if (rank_n < NPX*NPY) {
				MPI_Recv(BOUNDARY.p_Rs,NM*NX,MPI_DOUBLE,rank_n,4,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			}
			if (rank_s >= 0) {
				MPI_Send(BOUNDARY.p_Ss,NM*NX,MPI_DOUBLE,rank_s,4,MPI_COMM_WORLD);
			}
		}
		
	#endif // Done with North-south messages


	////////////////////////////////////////////////////////////////////
	// Step 4: Unload the data

	#if NPX>1 // Only bother if there are multiple tiles in the X direction
		// Unload points going east and west
		if (rank_w > -1) {
			for (m=0; m<NM; m++) {
				for (j=0; j<NY; j++) {
					p[m][j+1][0]=BOUNDARY.p_Re[m*NY+j];	// Data moving east lands in first slot
				}
			}
		}
		
		if (rank_e > -1) {
			for (m=0; m<NM; m++) {
				for (j=0; j<NY; j++) {
					p[m][j+1][NX+1]=BOUNDARY.p_Rw[m*NY+j];	// Data moving west lands in last slot
				}
			}
		}
	#endif

	#if NPY>1 // North-south messages: Only bother if there are multiple tiles in the Y direction
		// Unload points going north and south	
		if (rank_s >= 0){
			for (m=0; m<NM; m++) {
				for (i=0; i<NX; i++) {
					p[m][0][i+1]=BOUNDARY.p_Rn[m*NX+i]; 	// Data moving north lands in first slot 
				}
			}
		}
		
		if (rank_n < NPX*NPY) {
			for (m=0; m<NM; m++) {
				for (i=0; i<NX; i++) {
					p[m][NY+1][i+1]=BOUNDARY.p_Rs[m*NX+i]; // Data moving south lands in last slot
				}
			}
		}		
	#endif	

}
