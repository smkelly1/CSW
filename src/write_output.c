#include <netcdf.h>
#include "csw.h"

#ifdef WRITE_OUTPUT
void write_output(int sW, int rank)
{
	int i, j, n; 
	int ncid, status, varid, dimid[4], dimid2D[3];
	char name[100];
	double day;

	size_t start[]={sW, 0, 0, 0};
	size_t count[]={1, NMW, NY, NX};

	size_t start2D[]={sW, 0, 0};
	size_t count2D[]={1, NY, NX};

	// Get file name
	sprintf(name,FILE_OUT ".%03d.nc",rank);    


	////////////////////////////////////////////////////////////////////
	// Initiate the output file if this is the first call
	if (sW==0) {

		// Create file 
		if ((status = nc_create(name, NC_CLOBBER, &ncid)))
			ERR(status);

		// Create dimensions  
		if ((status = nc_def_dim(ncid, "x", NX, &dimid[3])))
			ERR(status);

		if ((status = nc_def_dim(ncid, "y", NY, &dimid[2])))
			ERR(status);

		if ((status = nc_def_dim(ncid, "mode", NMW, &dimid[1])))
			ERR(status);

		if ((status = nc_def_dim(ncid, "time", NC_UNLIMITED, &dimid[0])))
			ERR(status);

		dimid2D[0]=dimid[0];
		dimid2D[1]=dimid[2];
		dimid2D[2]=dimid[3];

		// Define the variables						
		if ((status = nc_def_var(ncid, "yday", NC_DOUBLE, 1, &dimid[0], &varid)))
			ERR(status);

		#ifdef WRITE_ETA
			if ((status = nc_def_var(ncid, "eta", NC_FLOAT, 4, dimid, &varid)))
				ERR(status);			
		#endif
		
		#ifdef DAMP_GROWTH			
			if ((status = nc_def_var(ncid, "eta_amp", NC_FLOAT, 3, dimid2D, &varid)))
				ERR(status);
				
			if ((status = nc_def_var(ncid, "eta_high", NC_FLOAT, 3, dimid2D, &varid)))
				ERR(status);
				
			if ((status = nc_def_var(ncid, "eta_low", NC_FLOAT, 3, dimid2D, &varid)))
				ERR(status);
				
			if ((status = nc_def_var(ncid, "flag_growth", NC_FLOAT, 3, dimid2D, &varid)))
				ERR(status);
		#endif

		#ifdef WRITE_MIX
			if ((status = nc_def_var(ncid, "u_mix", NC_FLOAT, 3, dimid2D, &varid)))
				ERR(status);
			if ((status = nc_def_var(ncid, "v_mix", NC_FLOAT, 3, dimid2D, &varid)))
				ERR(status);			
		#endif

		#ifdef WRITE_VELOCITY
			if ((status = nc_def_var(ncid, "u", NC_FLOAT, 4, dimid, &varid)))
				ERR(status);
		
			if ((status = nc_def_var(ncid, "v", NC_FLOAT, 4, dimid, &varid)))
				ERR(status);
		#endif

		#ifdef WRITE_WIND
			if ((status = nc_def_var(ncid, "TAUX", NC_FLOAT, 3, dimid2D, &varid)))
				ERR(status);
				
			if ((status = nc_def_var(ncid, "TAUY", NC_FLOAT, 3, dimid2D, &varid)))
				ERR(status);
		#endif
				
		// Close the file
		if ((status = nc_close(ncid)))
			ERR(status);
			
	} // Done initiating file 


	////////////////////////////////////////////////////////////////////
	// Open the file for writing
	if ((status = nc_open(name, NC_WRITE, &ncid)))
		ERR(status);

	#ifdef WRITE_VELOCITY
		// Write u velocity
		for(n=0; n<NMW; n++){
			for(j=0; j<NY; j++){
				for(i=0; i<NX; i++){
					tmp[n][j][i]=(float)((U[n][j+1][i+1]+U[n][j+1][i+2])/(2*H[j+1][i+1]));
				}
			}
		}

		if ((status = nc_inq_varid(ncid, "u", &varid)))
			ERR(status);

		if ((status = nc_put_vara_float(ncid, varid, start, count, &tmp[0][0][0])))
			ERR(status);
			
		// Write v velocity
		for(n=0; n<NMW; n++){
			for(j=0; j<NY; j++){
				for(i=0; i<NX; i++){
					tmp[n][j][i]=(float)((V[n][j+1][i+1]+V[n][j+2][i+1])/(2*H[j+1][i+1]));
				}
			}
		}

		if ((status = nc_inq_varid(ncid, "v", &varid)))
			ERR(status);

		if ((status = nc_put_vara_float(ncid, varid, start, count, &tmp[0][0][0])))
			ERR(status);
	#endif // WRITE_VELOCITY


	#ifdef WRITE_ETA
		// Write pressure
		for(n=0; n<NMW; n++){
			for(j=0; j<NY; j++){
				for(i=0; i<NX; i++){
					tmp[n][j][i]=(float)(p[n][j+1][i+1]*phi_surf[n][j+1][i+1]/9.81);
				}
			}
		}

		if ((status = nc_inq_varid(ncid, "eta", &varid)))
			ERR(status);

		if ((status = nc_put_vara_float(ncid, varid, start, count, &tmp[0][0][0])))
			ERR(status);			
	#endif // WRITE_ETA

	#ifdef WRITE_MIX
		// Write u_mix 
		for(j=0; j<NY; j++){
			for(i=0; i<NX; i++){
				tmp2D[j][i]=0;
				for(n=0; n<NMW; n++){
					tmp2D[j][i]=tmp2D[j][i]+(float)((U[n][j+1][i+1]+U[n][j+2][i+1])*phi_surf[n][j][i]/(2*H[j+1][i+1]));	
				}											
			}
		}

		if ((status = nc_inq_varid(ncid, "u_mix", &varid)))
			ERR(status);

		if ((status = nc_put_vara_float(ncid, varid, start2D, count2D, &tmp2D[0][0])))
			ERR(status);
			
		// Write v_mix	
		for(j=0; j<NY; j++){
			for(i=0; i<NX; i++){
				tmp2D[j][i]=0;
				for(n=0; n<NMW; n++){
					tmp2D[j][i]=tmp2D[j][i]+(float)((V[n][j+1][i+1]+V[n][j+2][i+1])*phi_surf[n][j][i]/(2*H[j+1][i+1]));	
				}											
			}
		}

		if ((status = nc_inq_varid(ncid, "v_mix", &varid)))
			ERR(status);

		if ((status = nc_put_vara_float(ncid, varid, start2D, count2D, &tmp2D[0][0])))
			ERR(status);		
	#endif	
	

	#ifdef DAMP_GROWTH
		// Write total SSH
		for(j=0; j<NY; j++){
			for(i=0; i<NX; i++){
				tmp2D[j][i]=(float)(etaA[j+1][i+1]);												
			}
		}

		if ((status = nc_inq_varid(ncid, "eta_amp", &varid)))
			ERR(status);

		if ((status = nc_put_vara_float(ncid, varid, start2D, count2D, &tmp2D[0][0])))
			ERR(status);
			
		// Write high frequency SSH
		for(j=0; j<NY; j++){
			for(i=0; i<NX; i++){
				tmp2D[j][i]=(float)(eta[j+1][i+1]-etaG[j+1][i+1]);												
			}
		}

		if ((status = nc_inq_varid(ncid, "eta_high", &varid)))
			ERR(status);

		if ((status = nc_put_vara_float(ncid, varid, start2D, count2D, &tmp2D[0][0])))
			ERR(status);
			
		// Write low frequency SSH
		for(j=0; j<NY; j++){
			for(i=0; i<NX; i++){
				tmp2D[j][i]=(float)(etaG[j+1][i+1]);												
			}
		}

		if ((status = nc_inq_varid(ncid, "eta_low", &varid)))
			ERR(status);

		if ((status = nc_put_vara_float(ncid, varid, start2D, count2D, &tmp2D[0][0])))
			ERR(status);
			
		// Write flag for extra damping
		for(j=0; j<NY; j++){
			for(i=0; i<NX; i++){
				tmp2D[j][i]=(float)(flag_growth[j+1][i+1]);												
			}
		}

		if ((status = nc_inq_varid(ncid, "flag_growth", &varid)))
			ERR(status);

		if ((status = nc_put_vara_float(ncid, varid, start2D, count2D, &tmp2D[0][0])))
			ERR(status);		
	#endif	


	#ifdef WRITE_WIND
		// Write TAUX
		for(j=0; j<NY; j++){
			for(i=0; i<NX; i++){
				tmp2D[j][i]=(float)(tau_x[j+1][i+1]);
			}
		}

		if ((status = nc_inq_varid(ncid, "TAUX", &varid)))
			ERR(status);

		if ((status = nc_put_vara_float(ncid, varid, start2D, count2D, &tmp2D[0][0])))
			ERR(status);
			
		// Write TAUY
		for(j=0; j<NY; j++){
			for(i=0; i<NX; i++){
				tmp2D[j][i]=(float)(tau_y[j+1][i+1]);
			}
		}

		if ((status = nc_inq_varid(ncid, "TAUY", &varid)))
			ERR(status);

		if ((status = nc_put_vara_float(ncid, varid, start2D, count2D, &tmp2D[0][0])))
			ERR(status);
	#endif // WRITE_WIND


	// Write time
	day=t/(24*3600);

	if ((status = nc_inq_varid(ncid, "yday", &varid)))
		ERR(status);

	if ((status = nc_put_vara_double(ncid, varid, &start[0], &count[0], &day)))
		ERR(status);


	// Close the file
	if ((status = nc_close(ncid)))
		ERR(status);
}
#endif
