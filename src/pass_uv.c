#include <mpi.h>
#include <csw.h>

void pass_uv(int rank)
{
	int i, j, m;
	int rank_e, rank_w, rank_n, rank_s, row;
	double u_Sew[NM*NY], u_Rew[NM*NY], u_Sns[NM*(NX+1)], u_Rns[NM*(NX+1)];
	double v_Sew[NM*(NY+1)], v_Rew[NM*(NY+1)], v_Sns[NM*NX], v_Rns[NM*NX];
	
	
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
			for (j=0; j<NY; j++) {
				u_Sew[m*NY+j]=U[m][j+1][NX];     // Third to last point
				v_Sew[m*(NY+1)+j]=V[m][j+1][NX]; // Second to last point				
			}
			v_Sew[m*(NY+1)+NY]=V[m][NY+1][NX+1]; // Add extra point on V boundary
		}
	
		// Odd points send first to avoid deadlock
		if (rank % 2) {
			if (rank_e > -1) {
				MPI_Send(u_Sew,NM*NY,MPI_DOUBLE,rank_e,11,MPI_COMM_WORLD);
			}
			if (rank_w > -1) {
				MPI_Recv(u_Rew,NM*NY,MPI_DOUBLE,rank_w,11,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			}
			
			if (rank_e > -1) {
				MPI_Send(v_Sew,NM*NY,MPI_DOUBLE,rank_e,12,MPI_COMM_WORLD);
			}
			if (rank_w > -1) {
				MPI_Recv(v_Rew,NM*NY,MPI_DOUBLE,rank_w,12,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			}
		}
		else {
			if (rank_w > -1) {
				MPI_Recv(u_Rew,NM*NY,MPI_DOUBLE,rank_w,11,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			}
			if (rank_e > -1) {
				MPI_Send(u_Sew,NM*NY,MPI_DOUBLE,rank_e,11,MPI_COMM_WORLD);
			}
	
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
				for (j=0; j<NY; j++) {
					U[m][j+1][0]=u_Rew[m*NY+j]; 
					V[m][j+1][0]=v_Rew[m*(NY+1)+j];				
				}
				V[m][NY+1][0]=v_Rew[m*(NY+1)+NY]; 
			}
		}
		
		////////////////////////////////////////////////////////////////
		// Send beginning points west
		for (m=0; m<NM; m++) {
			for (j=0; j<NY; j++) {
				u_Sew[m*NY+j]=U[m][j+1][2];     // Send third point
				v_Sew[m*(NY+1)+j]=V[m][j+1][1];	// Send second point			
			}
			v_Sew[m*(NY+1)+NY]=V[m][NY+1][1]; 
		}
			
		// Odd points send first to avoid deadlock
		if (rank % 2) {
			if (rank_w > -1) {
				MPI_Send(u_Sew,NM*NY,MPI_DOUBLE,rank_w,21,MPI_COMM_WORLD);
			}
			if (rank_e > -1) {
				MPI_Recv(u_Rew,NM*NY,MPI_DOUBLE,rank_e,21,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			}
			
			if (rank_w > -1) {
				MPI_Send(v_Sew,NM*NY,MPI_DOUBLE,rank_w,22,MPI_COMM_WORLD);
			}
			if (rank_e > -1) {
				MPI_Recv(v_Rew,NM*NY,MPI_DOUBLE,rank_e,22,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			}
		}
		else {
			if (rank_e > -1) {
				MPI_Recv(u_Rew,NM*NY,MPI_DOUBLE,rank_e,21,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			}
			if (rank_w > -1) {
				MPI_Send(u_Sew,NM*NY,MPI_DOUBLE,rank_w,21,MPI_COMM_WORLD);
			}
			
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
				for (j=0; j<NY; j++) {
					U[m][j+1][NX+2]=u_Rew[m*NY+j]; 
					V[m][j+1][NX+1]=v_Rew[m*(NY+1)+j];				
				}
				V[m][NY+1][NX+1]=v_Rew[m*(NY+1)+NY]; 
			}
		}
	#endif
	
	
	// The simple case when NPX==1
	#if NPX==1 && defined(PERIODICBC)
		for(m=0; m<NM; m++) {
			for (j=0; j<NY; j++) {
				U[m][j+1][NX+2]=U[m][j+1][2]; //third point -> last point
				U[m][j+1][0]=U[m][j+1][NX+1]; //third to last point -> first point 
				
				V[m][j+1][NX+1]=V[m][j+1][1]; //Second point -> last point
				V[m][j+1][0]=V[m][j+1][NX];   //Second to last point -> first point 
			}
			// Need to pass one more value for V
			V[m][NY+1][NX+1]=V[m][NY+1][1]; //Second point -> last point
			V[m][NY+1][0]=V[m][NY+1][NX];   //Second to last point -> first point 
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
			for (i=0; i<NX; i++) {
				u_Sns[m*(NX+1)+i]=U[m][NY][i+1]; // Second to last point
				v_Sns[m*NX+i]=V[m][NY][i+1];     // Third to last point		
			}
			u_Sns[m*(NX+1)+NX]=U[m][NY][NX+1]; 
		}
			
		// Odd rows send first to avoid deadlock
		if (row % 2) {
			if (rank_n < NPX*NPY) {
				MPI_Send(u_Sns,NM*NX,MPI_DOUBLE,rank_n,31,MPI_COMM_WORLD);
				MPI_Send(v_Sns,NM*NX,MPI_DOUBLE,rank_n,32,MPI_COMM_WORLD);
			}
			MPI_Recv(u_Rns,NM*NX,MPI_DOUBLE,rank_s,31,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			MPI_Recv(v_Rns,NM*NX,MPI_DOUBLE,rank_s,32,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		}
		else {
			if (rank_s >= 0){
				MPI_Recv(u_Rns,NM*NX,MPI_DOUBLE,rank_s,31,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				MPI_Recv(v_Rns,NM*NX,MPI_DOUBLE,rank_s,32,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			}
			if (rank_n < NPX*NPY) {
				MPI_Send(u_Sns,NM*NX,MPI_DOUBLE,rank_n,31,MPI_COMM_WORLD);
				MPI_Send(v_Sns,NM*NX,MPI_DOUBLE,rank_n,32,MPI_COMM_WORLD);
			}
		}
				
		// Receive first point from the south
		if (rank_s >= 0){
			for (m=0; m<NM; m++) {
				for (i=0; i<NX; i++) {
					U[m][0][i+1]=u_Rns[m*(NX+1)+i]; 
					V[m][0][i+1]=v_Rns[m*NX+i];					
				}
				U[m][0][NX+1]=u_Rns[m*(NX+1)+NX]; 
			}		
		}		
		
		//////////////////////////////////////////////////////////////////////
		// Send beginning points south
		for (m=0; m<NM; m++) {
			for (i=0; i<NX; i++) {
				u_Sns[m*(NX+1)+i]=U[m][1][i+1]; // Send second point
				v_Sns[m*NX+i]=V[m][2][i+1]; 	// Send third point	
			}
			u_Sns[m*(NX+1)+NX]=U[m][1][NX+1]; 
		}
			
		// Odd rows send first to avoid deadlock
		if (row % 2) {
			MPI_Send(u_Sns,NM*NX,MPI_DOUBLE,rank_s,41,MPI_COMM_WORLD);
			MPI_Send(v_Sns,NM*NX,MPI_DOUBLE,rank_s,42,MPI_COMM_WORLD);
	
			if (rank_n < NPX*NPY) {
				MPI_Recv(u_Rns,NM*NX,MPI_DOUBLE,rank_n,41,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				MPI_Recv(v_Rns,NM*NX,MPI_DOUBLE,rank_n,42,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			}
		}
		else {
			if (rank_n < NPX*NPY) {
				MPI_Recv(u_Rns,NM*NX,MPI_DOUBLE,rank_n,41,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				MPI_Recv(v_Rns,NM*NX,MPI_DOUBLE,rank_n,42,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			}
			if (rank_s >= 0){
				MPI_Send(u_Sns,NM*NX,MPI_DOUBLE,rank_s,41,MPI_COMM_WORLD);
				MPI_Send(v_Sns,NM*NX,MPI_DOUBLE,rank_s,42,MPI_COMM_WORLD);
			}
		}
						
		// Receive last point from the north
		if (rank_n < NPX*NPY) {
			for (m=0; m<NM; m++) {
				for (i=0; i<NX; i++) {
					U[m][NY+1][i+1]=u_Rns[m*(NX+1)+i]; 
					V[m][NY+2][i+1]=v_Rns[m*NX+i];
				}
				U[m][NY+1][NX+1]=u_Rns[m*(NX+1)+NX]; 
			}		
		}
	#endif
	
}
