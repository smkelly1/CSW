#include <math.h>
#include <mpi.h>
#include <csw.h>

void pass_p(int rank)
{

	int i, j, m;
	int rank_e, rank_w, rank_n, rank_s, row;
	double p_Sew[NM*NY], p_Rew[NM*NY], p_Sns[NM*NX], p_Rns[NM*NX];
	
	
	////////////////////////////////////////////////////////////////////
	// Identify east/west ranks
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
	// Send second to last point to east
	for (m=0; m<NM; m++) {
		for (j=0; j<NY; j++) {
			p_Sew[m*NY+j]=p1[m][j+1][NX]; // Length is NX+2, last point is NX+1 
		}
	}

	// Odd ranks send first to avoid deadlock
	if (rank % 2) {
		//printf("(1) Rank %d is sending east to %d\n",rank,rank_e);
		if (rank_e > -1) {
			MPI_Send(p_Sew,NM*NY,MPI_DOUBLE,rank_e,1,MPI_COMM_WORLD);
		}
		if (rank_w > -1) {
			MPI_Recv(p_Rew,NM*NY,MPI_DOUBLE,rank_w,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		}
	}
	else {
		//printf("(2) Rank %d is sending east to %d\n",rank,rank_e);
		if (rank_w > -1) {
			MPI_Recv(p_Rew,NM*NY,MPI_DOUBLE,rank_w,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		}
		if (rank_e > -1) {
			MPI_Send(p_Sew,NM*NY,MPI_DOUBLE,rank_e,1,MPI_COMM_WORLD);
		}
	}
		
	// Receive first point from west
	if (rank_w > -1) {
		for (m=0; m<NM; m++) {
			for (j=0; j<NY; j++) {
				p1[m][j+1][0]=p_Rew[m*NY+j]; // Write first point
			}
		}	
	}

	////////////////////////////////////////////////////////////////
	// Send second point to west
	for (m=0; m<NM; m++) {
		for (j=0; j<NY; j++) {
			p_Sew[m*NY+j]=p1[m][j+1][1]; // Send second point
		}
	}
		
	// Even ranks send first to avoid deadlock
	if (rank % 2) {
		//printf("(1a) Rank %d is sending west to %d\n",rank,rank_w);
		if (rank_w > -1) {
			MPI_Send(p_Sew,NM*NY,MPI_DOUBLE,rank_w,2,MPI_COMM_WORLD);
		}
		if (rank_e > -1) {
			MPI_Recv(p_Rew,NM*NY,MPI_DOUBLE,rank_e,2,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		}
		//printf("(2b) Rank %d has received from %d\n",rank,rank_e);
	}
	else {
		if (rank_e > -1) {
			MPI_Recv(p_Rew,NM*NY,MPI_DOUBLE,rank_e,2,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		}
		//printf("(1b) Rank %d has received from %d\n",rank,rank_e);
		//printf("(2a) Rank %d is sending west to %d\n",rank,rank_w);
		if (rank_w > -1) {
			MPI_Send(p_Sew,NM*NY,MPI_DOUBLE,rank_w,2,MPI_COMM_WORLD);
		}
	}
		
	// Receive last point from east
	if (rank_e > -1) {
		for (m=0; m<NM; m++) {
			for (j=0; j<NY; j++) {
				p1[m][j+1][NX+1]=p_Rew[m*NY+j]; 
			}
		}
	}


	////////////////////////////////////////////////////////////////////////
	//// Identify north/south ranks	
	rank_n=rank+NPX;
	rank_s=rank-NPX;
	row=rank/NPX;	

	
	////////////////////////////////////////////////////////////////////
	// Send second to last point to north
	for (m=0; m<NM; m++) {
		for (i=0; i<NX; i++) {
			p_Sns[m*NX+i]=p1[m][NY][i+1];  
		}
	}
		
	// Odd rows send first to avoid deadlock
	if (row % 2) {
		if (rank_n < NPX*NPY) {
			MPI_Send(p_Sns,NM*NX,MPI_DOUBLE,rank_n,3,MPI_COMM_WORLD);
		}
		MPI_Recv(p_Rns,NM*NX,MPI_DOUBLE,rank_s,3,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	}
	else {
		if (rank_s >= 0){
			MPI_Recv(p_Rns,NM*NX,MPI_DOUBLE,rank_s,3,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		}
		if (rank_n < NPX*NPY) {
			MPI_Send(p_Sns,NM*NX,MPI_DOUBLE,rank_n,3,MPI_COMM_WORLD);
		}
	}
		
	// Receive first point (if there is a cell to the south)
	if (rank_s >= 0){
		for (m=0; m<NM; m++) {
			for (i=0; i<NX; i++) {
				p1[m][0][i+1]=p_Rns[m*NX+i];  
			}
		}
	}	
	
	////////////////////////////////////////////////////////////////////
	// Send second point to south
		
	// Send second point
	for (m=0; m<NM; m++) {
		for (i=0; i<NX; i++) {
			p_Sns[m*NX+i]=p1[m][1][i+1]; 
		}
	}
		
	// Odd rows send first to avoid deadlock
	if (row % 2) {
		MPI_Send(p_Sns,NM*NX,MPI_DOUBLE,rank_s,4,MPI_COMM_WORLD);
		if (rank_n < NPX*NPY) {
			MPI_Recv(p_Rns,NM*NX,MPI_DOUBLE,rank_n,4,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		}
	}
	else {
		if (rank_n < NPX*NPY) {
			MPI_Recv(p_Rns,NM*NX,MPI_DOUBLE,rank_n,4,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		}
		if (rank_s >= 0) {
			MPI_Send(p_Sns,NM*NX,MPI_DOUBLE,rank_s,4,MPI_COMM_WORLD);
		}
	}
		
	// Receive last point (if there is a cell to the north)
	if (rank_n < NPX*NPY) {
		for (m=0; m<NM; m++) {
			for (i=0; i<NX; i++) {
				p1[m][NY+1][i+1]=p_Rns[m*NX+i]; 
			}
		}
	}
	
}
