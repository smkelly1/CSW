#include <mpi.h>
#include "csw.h"

void pass_uv(int rank)
{
	int i, j, m;
	int rank_e, rank_w, rank_n, rank_s, row;

	// Need to send u north/south and v east/west for Coriolis
	double u_Sns[NM*(NX+1)], u_Rns[NM*(NX+1)];
	double v_Sew[NM*(NY+1)], v_Rew[NM*(NY+1)];

	// Only need to send u east/west and v north/south for diffusion
	#ifdef AX
		double u_Sew[NM*NY], u_Rew[NM*NY];
		double v_Sns[NM*NX], v_Rns[NM*NX];
	#endif

	////////////////////////////////////////////////////////////////////
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

		////////////////////////////////////////////////////////////////////
		// Send end points east
		for (m=0; m<NM; m++) {
			#ifdef AX
			for (j=0; j<NY; j++) {
				u_Sew[m*NY+j]=U[m][j+1][NX];     // Third to last point
			}
			#endif

			for (j=0; j<NY+1; j++) {
				v_Sew[m*(NY+1)+j]=V[m][j+1][NX]; // Second to last point
			}
		}

		// Odd points send first to avoid deadlock
		if (rank % 2) {
			#ifdef AX
				if (rank_e > -1) {
					MPI_Send(u_Sew,NM*NY,MPI_DOUBLE,rank_e,11,MPI_COMM_WORLD);
				}
				if (rank_w > -1) {
					MPI_Recv(u_Rew,NM*NY,MPI_DOUBLE,rank_w,11,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				}
			#endif

			if (rank_e > -1) {
				MPI_Send(v_Sew,NM*NY,MPI_DOUBLE,rank_e,12,MPI_COMM_WORLD);
			}
			if (rank_w > -1) {
				MPI_Recv(v_Rew,NM*NY,MPI_DOUBLE,rank_w,12,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			}
		}
		else {
			#ifdef AX
				if (rank_w > -1) {
					MPI_Recv(u_Rew,NM*NY,MPI_DOUBLE,rank_w,11,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				}
				if (rank_e > -1) {
					MPI_Send(u_Sew,NM*NY,MPI_DOUBLE,rank_e,11,MPI_COMM_WORLD);
				}
			#endif

			if (rank_w > -1) {
				MPI_Recv(v_Rew,NM*NY,MPI_DOUBLE,rank_w,12,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			}
			if (rank_e > -1) {
				MPI_Send(v_Sew,NM*NY,MPI_DOUBLE,rank_e,12,MPI_COMM_WORLD);
			}
		}

		// Receive first point from west
		if (rank_w > -1) {
			for (m=0; m<NM; m++) {
				#ifdef AX
					for (j=0; j<NY; j++) {
						U[m][j+1][0]=u_Rew[m*NY+j]; 
					}
				#endif
				
				for (j=0; j<NY+1; j++) {
					V[m][j+1][0]=v_Rew[m*(NY+1)+j];
				}
			}
		}

		////////////////////////////////////////////////////////////////
		// Send beginning points west
		for (m=0; m<NM; m++) {
			#ifdef AX
				for (j=0; j<NY; j++) {
					u_Sew[m*NY+j]=U[m][j+1][2]; // Send third point
				}
			#endif

			for (j=0; j<NY+1; j++) {
				v_Sew[m*(NY+1)+j]=V[m][j+1][1]; // Send second point
			}
		}

		// Odd points send first to avoid deadlock
		if (rank % 2) {
			#ifdef AX
				if (rank_w > -1) {
					MPI_Send(u_Sew,NM*NY,MPI_DOUBLE,rank_w,21,MPI_COMM_WORLD);
				}
				if (rank_e > -1) {
					MPI_Recv(u_Rew,NM*NY,MPI_DOUBLE,rank_e,21,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				}
			#endif

			if (rank_w > -1) {
				MPI_Send(v_Sew,NM*NY,MPI_DOUBLE,rank_w,22,MPI_COMM_WORLD);
			}
			if (rank_e > -1) {
				MPI_Recv(v_Rew,NM*NY,MPI_DOUBLE,rank_e,22,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			}
		}
		else {
			#ifdef AX
				if (rank_e > -1) {
					MPI_Recv(u_Rew,NM*NY,MPI_DOUBLE,rank_e,21,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				}
				if (rank_w > -1) {
					MPI_Send(u_Sew,NM*NY,MPI_DOUBLE,rank_w,21,MPI_COMM_WORLD);
				}
			#endif
			
			if (rank_e > -1) {
				MPI_Recv(v_Rew,NM*NY,MPI_DOUBLE,rank_e,22,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			}
			if (rank_w > -1) {
				MPI_Send(v_Sew,NM*NY,MPI_DOUBLE,rank_w,22,MPI_COMM_WORLD);
			}
		}

		// Receive last point from east
		if (rank_e > -1) {
			for (m=0; m<NM; m++) {
				#ifdef AX
					for (j=0; j<NY; j++) {
						U[m][j+1][NX+2]=u_Rew[m*NY+j]; 
					}
				#endif

				for (j=0; j<NY+1; j++) {
					V[m][j+1][NX+1]=v_Rew[m*(NY+1)+j];
				}
			}
		}
	#endif


	// The simple case when NPX==1
	#if NPX==1 && defined(PERIODICBC)
		for(m=0; m<NM; m++) {
			for (j=0; j<NY; j++) {
				U[m][j+1][NX+2]=U[m][j+1][2]; //third point -> last point
				U[m][j+1][0]=U[m][j+1][NX+1]; //third to last point -> first point
			}
			for (j=0; j<NY+1; j++) {
				V[m][j+1][NX+1]=V[m][j+1][1]; //Second point -> last point
				V[m][j+1][0]=V[m][j+1][NX]; //Second to last point -> first point 
			}
		}
	#endif


	////////////////////////////////////////////////////////////////////
	//// Identify north/south ranks	
	#if NPY>1 // Only bother if there are multiple tiles in the Y direction

		rank_n=rank+NPX;
		rank_s=rank-NPX;
		row=rank/NPX;


		////////////////////////////////////////////////////////////////////
		// Send end points north
		for (m=0; m<NM; m++) {
			for (i=0; i<NX+1; i++) {
				u_Sns[m*(NX+1)+i]=U[m][NY][i+1]; // Second to last point
			}

			#ifdef AX
				for (i=0; i<NX; i++) {
					v_Sns[m*NX+i]=V[m][NY][i+1]; // Third to last point
				}
			#endif
		}

		// Odd rows send first to avoid deadlock
		if (row % 2) {
			if (rank_n < NPX*NPY) {
				MPI_Send(u_Sns,NM*NX,MPI_DOUBLE,rank_n,31,MPI_COMM_WORLD);
				#ifdef AX
					MPI_Send(v_Sns,NM*NX,MPI_DOUBLE,rank_n,32,MPI_COMM_WORLD);
				#endif
			}
			MPI_Recv(u_Rns,NM*NX,MPI_DOUBLE,rank_s,31,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			#ifdef AX
				MPI_Recv(v_Rns,NM*NX,MPI_DOUBLE,rank_s,32,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			#endif
		}
		else {
			if (rank_s >= 0){
				MPI_Recv(u_Rns,NM*NX,MPI_DOUBLE,rank_s,31,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				#ifdef AX
					MPI_Recv(v_Rns,NM*NX,MPI_DOUBLE,rank_s,32,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				#endif
			}
			if (rank_n < NPX*NPY) {
				MPI_Send(u_Sns,NM*NX,MPI_DOUBLE,rank_n,31,MPI_COMM_WORLD);
				#ifdef AX
					MPI_Send(v_Sns,NM*NX,MPI_DOUBLE,rank_n,32,MPI_COMM_WORLD);
				#endif
			}
		}

		// Receive first point from the south
		if (rank_s >= 0){
			for (m=0; m<NM; m++) {
				for (i=0; i<NX+1; i++) {
					U[m][0][i+1]=u_Rns[m*(NX+1)+i]; 
				}

				#ifdef AX
					for (i=0; i<NX; i++) {
						V[m][0][i+1]=v_Rns[m*NX+i];
					}
				#endif
			}
		}

		//////////////////////////////////////////////////////////////////////
		// Send beginning points south
		for (m=0; m<NM; m++) {
			for (i=0; i<NX+1; i++) {
				u_Sns[m*(NX+1)+i]=U[m][1][i+1]; // Send second point
			}

			#ifdef AX
				for (i=0; i<NX; i++) {
					v_Sns[m*NX+i]=V[m][2][i+1]; // Send third point
				}
			#endif
		}

		// Odd rows send first to avoid deadlock
		if (row % 2) {
			MPI_Send(u_Sns,NM*NX,MPI_DOUBLE,rank_s,41,MPI_COMM_WORLD);
			#ifdef AX
				MPI_Send(v_Sns,NM*NX,MPI_DOUBLE,rank_s,42,MPI_COMM_WORLD);
			#endif

			if (rank_n < NPX*NPY) {
				MPI_Recv(u_Rns,NM*NX,MPI_DOUBLE,rank_n,41,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				#ifdef AX
					MPI_Recv(v_Rns,NM*NX,MPI_DOUBLE,rank_n,42,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				#endif
			}
		}
		else {
			if (rank_n < NPX*NPY) {
				MPI_Recv(u_Rns,NM*NX,MPI_DOUBLE,rank_n,41,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				#ifdef AX
					MPI_Recv(v_Rns,NM*NX,MPI_DOUBLE,rank_n,42,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				#endif
			}
			if (rank_s >= 0){
				MPI_Send(u_Sns,NM*NX,MPI_DOUBLE,rank_s,41,MPI_COMM_WORLD);
				#ifdef AX
					MPI_Send(v_Sns,NM*NX,MPI_DOUBLE,rank_s,42,MPI_COMM_WORLD);
				#endif
			}
		}

		// Receive last point from the north
		if (rank_n < NPX*NPY) {
			for (m=0; m<NM; m++) {
				for (i=0; i<NX+1; i++) {
					U[m][NY+1][i+1]=u_Rns[m*(NX+1)+i]; 
				}

				#ifdef AX
					for (i=0; i<NX; i++) {
						V[m][NY+2][i+1]=v_Rns[m*NX+i];
					}
				#endif
			}
		}
	#endif

}
