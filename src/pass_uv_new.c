#include <mpi.h>
#include "csw.h"

void pass_uv(int rank)
{
	int i, j, m;
	int rank_e, rank_w, rank_n, rank_s, row;

	// NOTE: Both sent and received data is labeled by it's direction of travel 
	
		
	////////////////////////////////////////////////////////////////////
	// Step 0: implement periodic boundaries if we don't need MPI
	#if NPX==1 && defined(PERIODICBC) //The simple case when NPX==1
		for(m=0; m<NM; m++) {
			for (j=0; j<NY; j++) {
				U[m][j+1][NX+2]=U[m][j+1][2]; //third point -> last point
				U[m][j+1][0]=U[m][j+1][NX+1]; //third to last point -> first point
			}
			for (j=0; j<NY+1; j++) {
				V[m][j+1][NX+1]=V[m][j+1][1]; //Second point -> last point
				V[m][j+1][0]=V[m][j+1][NX];   //Second to last point -> first point 
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
			#ifdef PERIODICBC
				rank_e=(rank-NPX)+1;
			#else
				rank_e=-1;
			#endif
		}

		if (rank % NPX == 0) {
			#ifdef PERIODICBC
				rank_w=(rank+NPX)-1;
			#else
				rank_w=-1;
			#endif
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
			for (j=0; j<NY+1; j++) {
				BOUNDARY.v_Se[m*(NY+1)+j]=V[m][j+1][NX]; // Second to last point goes east
				BOUNDARY.v_Sw[m*(NY+1)+j]=V[m][j+1][1];  // Second point goes west
			}
			
			#ifdef AX
				for (j=0; j<NY; j++) {
					BOUNDARY.u_Se[m*NY+j]=U[m][j+1][NX]; // Third to last point goes east
					BOUNDARY.u_Sw[m*NY+j]=U[m][j+1][2];  // Third point goes west
				}
			#endif
		}		
	#endif

	#if NPY>1 // Only bother if there are multiple tiles in the Y direction
		// Send end points north
		for (m=0; m<NM; m++) {
			for (i=0; i<NX+1; i++) {
				BOUNDARY.u_Sn[m*(NX+1)+i]=U[m][NY][i+1]; // Second to last point goes north
				BOUNDARY.u_Ss[m*(NX+1)+i]=U[m][1][i+1];  // Second point goes south

			}

			#ifdef AX
				for (i=0; i<NX; i++) {
					BOUNDARY.v_Sn[m*NX+i]=V[m][NY][i+1]; // Third to last point goes north
					BOUNDARY.v_Ss[m*NX+i]=V[m][2][i+1];  // Third point goes south
				}
			#endif
		}		
	#endif
	
	
	////////////////////////////////////////////////////////////////////
	// Step 3: Mail the envelopes
	
	#if NPX>1 // East-west messages: Only bother if there are multiple tiles in the X direction

		// First send data east, odd ranks send first to avoid deadlock
		if (rank % 2) {
			if (rank_e > -1) {
				MPI_Send(BOUNDARY.v_Se,NM*(NY+1),MPI_DOUBLE,rank_e,12,MPI_COMM_WORLD);
			}
			if (rank_w > -1) {
				MPI_Recv(BOUNDARY.v_Re,NM*(NY+1),MPI_DOUBLE,rank_w,12,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			}
			#ifdef AX
				if (rank_e > -1) {
					MPI_Send(BOUNDARY.u_Se,NM*NY,MPI_DOUBLE,rank_e,11,MPI_COMM_WORLD);
				}
				if (rank_w > -1) {
					MPI_Recv(BOUNDARY.u_Re,NM*NY,MPI_DOUBLE,rank_w,11,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				}
			#endif
		}
		else {
			if (rank_w > -1) {
				MPI_Recv(BOUNDARY.v_Re,NM*(NY+1),MPI_DOUBLE,rank_w,12,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			}
			if (rank_e > -1) {
				MPI_Send(BOUNDARY.v_Se,NM*(NY+1),MPI_DOUBLE,rank_e,12,MPI_COMM_WORLD);
			}
			#ifdef AX
				if (rank_w > -1) {
					MPI_Recv(BOUNDARY.u_Re,NM*NY,MPI_DOUBLE,rank_w,11,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				}
				if (rank_e > -1) {
					MPI_Send(BOUNDARY.u_Se,NM*NY,MPI_DOUBLE,rank_e,11,MPI_COMM_WORLD);
				}
			#endif
		}
		
		// Second send data west, odd ranks send first to avoid deadlock
		if (rank % 2) {
			if (rank_w > -1) {
				MPI_Send(BOUNDARY.v_Sw,NM*(NY+1),MPI_DOUBLE,rank_w,22,MPI_COMM_WORLD);
			}
			if (rank_e > -1) {
				MPI_Recv(BOUNDARY.v_Rw,NM*(NY+1),MPI_DOUBLE,rank_e,22,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			}
			#ifdef AX
				if (rank_w > -1) {
					MPI_Send(BOUNDARY.u_Sw,NM*NY,MPI_DOUBLE,rank_w,21,MPI_COMM_WORLD);
				}
				if (rank_e > -1) {
					MPI_Recv(BOUNDARY.u_Rw,NM*NY,MPI_DOUBLE,rank_e,21,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				}
			#endif
		}
		else {	
			if (rank_e > -1) {
				MPI_Recv(BOUNDARY.v_Rw,NM*(NY+1),MPI_DOUBLE,rank_e,22,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			}
			if (rank_w > -1) {
				MPI_Send(BOUNDARY.v_Sw,NM*(NY+1),MPI_DOUBLE,rank_w,22,MPI_COMM_WORLD);
			}
			#ifdef AX
				if (rank_e > -1) {
					MPI_Recv(BOUNDARY.u_Rw,NM*NY,MPI_DOUBLE,rank_e,21,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				}
				if (rank_w > -1) {
					MPI_Send(BOUNDARY.u_Sw,NM*NY,MPI_DOUBLE,rank_w,21,MPI_COMM_WORLD);
				}
			#endif
		}	
			
	#endif // Done with East-west messages


	#if NPY>1 // North-south messages: Only bother if there are multiple tiles in the Y direction

		// First send data north, odd rows send first to avoid deadlock
		if (row % 2) {
			if (rank_n < NPX*NPY) {
				MPI_Send(BOUNDARY.u_Sn,NM*(NX+1),MPI_DOUBLE,rank_n,31,MPI_COMM_WORLD);
				#ifdef AX
					MPI_Send(BOUNDARY.v_Sn,NM*NX,MPI_DOUBLE,rank_n,32,MPI_COMM_WORLD);
				#endif
			}
			MPI_Recv(BOUNDARY.u_Rn,NM*(NX+1),MPI_DOUBLE,rank_s,31,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			#ifdef AX
				MPI_Recv(BOUNDARY.v_Rn,NM*NX,MPI_DOUBLE,rank_s,32,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			#endif
		}
		else {
			if (rank_s >= 0){
				MPI_Recv(BOUNDARY.u_Rn,NM*(NX+1),MPI_DOUBLE,rank_s,31,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				#ifdef AX
					MPI_Recv(BOUNDARY.v_Rn,NM*NX,MPI_DOUBLE,rank_s,32,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				#endif
			}
			if (rank_n < NPX*NPY) {
				MPI_Send(BOUNDARY.u_Sn,NM*(NX+1),MPI_DOUBLE,rank_n,31,MPI_COMM_WORLD);
				#ifdef AX
					MPI_Send(BOUNDARY.v_Sn,NM*NX,MPI_DOUBLE,rank_n,32,MPI_COMM_WORLD);
				#endif
			}
		}
		
		// Second send data south, Odd rows send first to avoid deadlock
		if (row % 2) {
			MPI_Send(BOUNDARY.u_Ss,NM*(NX+1),MPI_DOUBLE,rank_s,41,MPI_COMM_WORLD);
			#ifdef AX
				MPI_Send(BOUNDARY.v_Ss,NM*NX,MPI_DOUBLE,rank_s,42,MPI_COMM_WORLD);
			#endif

			if (rank_n < NPX*NPY) {
				MPI_Recv(BOUNDARY.u_Rs,NM*(NX+1),MPI_DOUBLE,rank_n,41,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				#ifdef AX
					MPI_Recv(BOUNDARY.v_Rs,NM*NX,MPI_DOUBLE,rank_n,42,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				#endif
			}
		}
		else {
			if (rank_n < NPX*NPY) {
				MPI_Recv(BOUNDARY.u_Rs,NM*(NX+1),MPI_DOUBLE,rank_n,41,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				#ifdef AX
					MPI_Recv(BOUNDARY.v_Rs,NM*NX,MPI_DOUBLE,rank_n,42,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				#endif
			}
			if (rank_s >= 0){
				MPI_Send(BOUNDARY.u_Ss,NM*(NX+1),MPI_DOUBLE,rank_s,41,MPI_COMM_WORLD);
				#ifdef AX
					MPI_Send(BOUNDARY.v_Ss,NM*NX,MPI_DOUBLE,rank_s,42,MPI_COMM_WORLD);
				#endif
			}
		}
		
	#endif // Done with North-south messages


	////////////////////////////////////////////////////////////////////
	// Step 4: Unload the data
	#if NPX>1 // Only bother if there are multiple tiles in the X direction
		// Unload points going east and west
		if (rank_w > -1) {
			for (m=0; m<NM; m++) {						
				for (j=0; j<NY+1; j++) {
					V[m][j+1][0]=BOUNDARY.v_Re[m*(NY+1)+j];     // Data moving east lands in first slot
				}
				#ifdef AX
					for (j=0; j<NY; j++) {
						U[m][j+1][0]=BOUNDARY.u_Re[m*NY+j];     // Data moving east lands in first slot
					}
				#endif
			}
		}	
		
		if (rank_e > -1) {
			for (m=0; m<NM; m++) {						
				for (j=0; j<NY+1; j++) {
					V[m][j+1][NX+1]=BOUNDARY.v_Rw[m*(NY+1)+j];  // Data moving west lands in last slot
				}
				#ifdef AX
					for (j=0; j<NY; j++) {
						U[m][j+1][NX+2]=BOUNDARY.u_Rw[m*NY+j];  // Data moving west lands in last slot
					}
				#endif
			}
		}		
	#endif

	#if NPY>1 // Only bother if there are multiple tiles in the Y direction
		// Unload points going north and south
		if (rank_s >= 0){
			for (m=0; m<NM; m++) {
				for (i=0; i<NX+1; i++) {
					U[m][0][i+1]=BOUNDARY.u_Rn[m*(NX+1)+i];              // Data moving north lands in first slot
				}

				#ifdef AX
					for (i=0; i<NX; i++) {
						V[m][0][i+1]=BOUNDARY.v_Rn[m*NX+i];              // Data moving north lands in first slot
					}
				#endif
			}
		}
		
		if (rank_n < NPX*NPY) {
			for (m=0; m<NM; m++) {
				for (i=0; i<NX+1; i++) {
					U[m][NY+1][i+1]=BOUNDARY.u_Rs[m*(NX+1)+i];           // Data moving south lands in last slot
				}

				#ifdef AX
					for (i=0; i<NX; i++) {
						V[m][NY+2][i+1]=BOUNDARY.v_Rs[m*NX+i];           // Data moving south lands in last slot
					}
				#endif
			}
		}			
	#endif

}
