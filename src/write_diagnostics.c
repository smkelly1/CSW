#include <netcdf.h>
#include <math.h>
#include "csw.h"

void write_diagnostics(int sD, int Na, int rank)
{
	int i, j, n; 
	int ncid, status, varid, dimid[4];
	char name[100];

	size_t start[]={sD-1, 0, 0, 0};
	size_t count[]={1, NMW, NY, NX};

	// Get file name
	sprintf(name,FILE_DIAG ".%03d.nc",rank); 

	////////////////////////////////////////////////////////////////////
	// Initiate the output file if this is the first call
	if (sD==1) {
		
		// Create file 
		if ((status = nc_create(name, NC_CLOBBER, &ncid)))
			ERR(status);
			
		// Create dimensions  
		if ((status = nc_def_dim(ncid, "longitude", NX, &dimid[3])))
			ERR(status);

		if ((status = nc_def_dim(ncid, "latitude", NY, &dimid[2])))
			ERR(status);

		if ((status = nc_def_dim(ncid, "mode", NMW, &dimid[1])))
			ERR(status);

		if ((status = nc_def_dim(ncid, "time", NC_UNLIMITED, &dimid[0])))
			ERR(status);
			
		// Create variables							
		if ((status = nc_def_var(ncid, "period", NC_DOUBLE, 1, &dimid[0], &varid)))
			ERR(status);

		#if defined(FLAG_GROWTH) 
			if ((status = nc_def_var(ncid, "eta_low", NC_FLOAT, 4, dimid, &varid)))
				ERR(status);
		#endif

		#ifdef WRITE_SSH
			if ((status = nc_def_var(ncid, "SSH_amp", NC_FLOAT, 4, dimid, &varid)))
				ERR(status);

			if ((status = nc_def_var(ncid, "SSH_phase", NC_FLOAT, 4, dimid, &varid)))
				ERR(status);
		#endif

		#ifdef WRITE_TRANSPORT
			if ((status = nc_def_var(ncid, "U_amp", NC_FLOAT, 4, dimid, &varid)))
				ERR(status);

			if ((status = nc_def_var(ncid, "U_phase", NC_FLOAT, 4, dimid, &varid)))
				ERR(status);
			
			if ((status = nc_def_var(ncid, "V_amp", NC_FLOAT, 4, dimid, &varid)))
				ERR(status);

			if ((status = nc_def_var(ncid, "V_phase", NC_FLOAT, 4, dimid, &varid)))
				ERR(status);
		#endif

		#ifdef ENERGY
			if ((status = nc_def_var(ncid, "KE", NC_FLOAT, 4, dimid, &varid)))
				ERR(status);
		
			if ((status = nc_def_var(ncid, "PE", NC_FLOAT, 4, dimid, &varid)))
				ERR(status);
		#endif
		
		#ifdef FLUX
			if ((status = nc_def_var(ncid, "up", NC_FLOAT, 4, dimid, &varid)))
				ERR(status);

			if ((status = nc_def_var(ncid, "vp", NC_FLOAT, 4, dimid, &varid)))
				ERR(status);
		#endif

		#ifdef WORK  
			if ((status = nc_def_var(ncid, "W", NC_FLOAT, 4, dimid, &varid)))
				ERR(status);
				
			if ((status = nc_def_var(ncid, "C0", NC_FLOAT, 4, dimid, &varid)))
				ERR(status);

			if ((status = nc_def_var(ncid, "Cn", NC_FLOAT, 4, dimid, &varid)))
				ERR(status);

			if ((status = nc_def_var(ncid, "D", NC_FLOAT, 4, dimid, &varid)))
				ERR(status);

			if ((status = nc_def_var(ncid, "divF", NC_FLOAT, 4, dimid, &varid)))
				ERR(status);

			if ((status = nc_def_var(ncid, "error", NC_FLOAT, 4, dimid, &varid)))
				ERR(status);
		#endif
		
		// Close the file
		if ((status = nc_close(ncid)))
			ERR(status);
			
	} // Done initiating file 


	////////////////////////////////////////////////////////////////////
   	// Open the file and write data
	if ((status = nc_open(name, NC_WRITE, &ncid)))
		ERR(status);

	// Write Flag for exponential growth
	#ifdef FLAG_GROWTH
		for(i=0; i<NX; i++){
			for(j=0; j<NY; j++){
				for(n=0; n<NMW; n++){					
					tmp[n][j][i]=(float)(p_low[n][j+1][i+1]*phi_surf[n][j+1][i+1]/9.81);					
				}
			}
		}

		if ((status = nc_inq_varid(ncid, "eta_low", &varid)))
			ERR(status);

		if ((status = nc_put_vara_float(ncid, varid, start, count, &tmp[0][0][0])))
			ERR(status);		
	#endif

	////////////////////////////////////////////////////////////////////
	// Average and write diagnostic fields  
	for(n=0; n<NMW; n++){
		for(j=0; j<NY; j++){
			for(i=0; i<NX; i++){

				#ifdef ENERGY
					KE[n][j][i]=KE[n][j][i]/Na;
					PE[n][j][i]=PE[n][j][i]/Na;
				#endif

				#ifdef FLUX
					up[n][j][i]=up[n][j][i]/Na;
					vp[n][j][i]=vp[n][j][i]/Na;
				#endif

				#ifdef WORK
					W[n][j][i]=W[n][j][i]/Na;
					C0[n][j][i]=C0[n][j][i]/Na;
					Cn[n][j][i]=Cn[n][j][i]/Na;
					D[n][j][i]=D[n][j][i]/Na;
					divF[n][j][i]=divF[n][j][i]/Na;

					tmp[n][j][i]=W[n][j][i]+C0[n][j][i]+Cn[n][j][i]-divF[n][j][i]-D[n][j][i]; // residual/error
				#endif

				#ifdef WRITE_SSH
					SSH_amp[n][j][i]=SSH_amp[n][j][i]*phi_surf[n][j+1][i+1]/9.81; // Convert from pressure to SSH
					SSH_phase[n][j][i]=SSH_phase[n][j][i]/Na*360; // This is time of high tide
				#endif

				#ifdef WRITE_TRANSPORT
					U_phase[n][j][i]=U_phase[n][j][i]/Na*360; // This is time of strongest (positive) flow
					V_phase[n][j][i]=V_phase[n][j][i]/Na*360; // This is time of strongest (positive) flow
				#endif
			}
		}
	}
	
	#ifdef ENERGY
		if ((status = nc_inq_varid(ncid, "KE", &varid)))
			ERR(status);

		if ((status = nc_put_vara_float(ncid, varid, start, count, &KE[0][0][0])))
			ERR(status);

		if ((status = nc_inq_varid(ncid, "PE", &varid)))
			ERR(status);

		if ((status = nc_put_vara_float(ncid, varid, start, count, &PE[0][0][0])))
			ERR(status);
	#endif

	#ifdef FLUX
		if ((status = nc_inq_varid(ncid, "up", &varid)))
			ERR(status);

		if ((status = nc_put_vara_float(ncid, varid, start, count, &up[0][0][0])))
			ERR(status);

		if ((status = nc_inq_varid(ncid, "vp", &varid)))
			ERR(status);

		if ((status = nc_put_vara_float(ncid, varid, start, count, &vp[0][0][0])))
			ERR(status);
	#endif

	#ifdef WORK
		if ((status = nc_inq_varid(ncid, "W", &varid)))
			ERR(status);

		if ((status = nc_put_vara_float(ncid, varid, start, count, &W[0][0][0])))
			ERR(status);
			
		if ((status = nc_inq_varid(ncid, "C0", &varid)))
			ERR(status);

		if ((status = nc_put_vara_float(ncid, varid, start, count, &C0[0][0][0])))
			ERR(status);

		if ((status = nc_inq_varid(ncid, "Cn", &varid)))
			ERR(status);

		if ((status = nc_put_vara_float(ncid, varid, start, count, &Cn[0][0][0])))
			ERR(status);

		if ((status = nc_inq_varid(ncid, "D", &varid)))
			ERR(status);

		if ((status = nc_put_vara_float(ncid, varid, start, count, &D[0][0][0])))
			ERR(status);

		if ((status = nc_inq_varid(ncid, "divF", &varid)))
			ERR(status);

		if ((status = nc_put_vara_float(ncid, varid, start, count, &divF[0][0][0])))
			ERR(status);

		if ((status = nc_inq_varid(ncid, "error", &varid)))
			ERR(status);

		if ((status = nc_put_vara_float(ncid, varid, start, count, &tmp[0][0][0])))
			ERR(status);
	#endif

	#ifdef WRITE_SSH
		if ((status = nc_inq_varid(ncid, "SSH_amp", &varid)))
			ERR(status);
	
		if ((status = nc_put_vara_float(ncid, varid, start, count, &SSH_amp[0][0][0])))
			ERR(status);

		if ((status = nc_inq_varid(ncid, "SSH_phase", &varid)))
			ERR(status);

		if ((status = nc_put_vara_float(ncid, varid, start, count, &SSH_phase[0][0][0])))
			ERR(status);
	#endif
	
	#ifdef WRITE_TRANSPORT
		if ((status = nc_inq_varid(ncid, "U_amp", &varid)))
			ERR(status);
	
		if ((status = nc_put_vara_float(ncid, varid, start, count, &U_amp[0][0][0])))
			ERR(status);

		if ((status = nc_inq_varid(ncid, "U_phase", &varid)))
			ERR(status);

		if ((status = nc_put_vara_float(ncid, varid, start, count, &U_phase[0][0][0])))
			ERR(status);
			
		if ((status = nc_inq_varid(ncid, "V_amp", &varid)))
			ERR(status);
	
		if ((status = nc_put_vara_float(ncid, varid, start, count, &V_amp[0][0][0])))
			ERR(status);

		if ((status = nc_inq_varid(ncid, "V_phase", &varid)))
			ERR(status);

		if ((status = nc_put_vara_float(ncid, varid, start, count, &V_phase[0][0][0])))
			ERR(status);
	#endif

	// Write time  
	if ((status = nc_inq_varid(ncid, "period", &varid)))
		ERR(status);

	if ((status = nc_put_vara_int(ncid, varid, &start[0], &count[0], &sD)))
		ERR(status);


	// Close the file
	if ((status = nc_close(ncid)))
		ERR(status);


	////////////////////////////////////////////////////////////////////	
	// Clear diagnostic fields
	for(n=0; n<NMW; n++){
		for(j=0; j<NY; j++){
			for(i=0; i<NX; i++){

				#ifdef ENERGY
					KE[n][j][i]=0;
					PE[n][j][i]=0;
				#endif

				#ifdef FLUX
					up[n][j][i]=0;
					vp[n][j][i]=0;
				#endif

				#ifdef WORK
					W[n][j][i]=0;
					C0[n][j][i]=0;
					Cn[n][j][i]=0;
					divF[n][j][i]=0;
					D[n][j][i]=0;
				#endif

				#ifdef WRITE_SSH
					SSH_amp[n][j][i]=0;
					SSH_phase[n][j][i]=0;
				#endif

				#ifdef WRITE_TRANSPORT
					U_amp[n][j][i]=0;
					U_phase[n][j][i]=0;
					V_amp[n][j][i]=0;
					V_phase[n][j][i]=0;
				#endif

			}
		}
	}


}
