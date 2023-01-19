#include <netcdf.h>
#include "csw.h"


void write_output(int sW, double t, int rank)
{
	#if defined(WRITE_PRESSURE) || defined(WRITE_VELOCITY)
		int i, j, n; 
	#endif
	
	int ncid, status, varid;
	char name[100];
	double day;

	size_t start[]={sW, 0, 0, 0};
	size_t count[]={1, NMW, NY, NX};

	// Get file name
	sprintf(name,FILE_OUT ".%03d.nc",rank);    

	// Open the file
	if ((status = nc_open(name, NC_WRITE, &ncid)))
		ERR(status);

	#ifdef WRITE_VELOCITY

		// Write u velocity
		for(i=0; i<NX; i++){
			for(j=0; j<NY; j++){
				for(n=0; n<NMW; n++){
					tmp[n][j][i]=(float)((U[n][j+1][i+1]+U[n][j+1][i+2])/(2*H[j+1][i+1]));
				}
			}
		}

		if ((status = nc_inq_varid(ncid, "u", &varid)))
			ERR(status);

		if ((status = nc_put_vara_float(ncid, varid, start, count, &tmp[0][0][0])))
			ERR(status);
			
		// Write v velocity
		for(i=0; i<NX; i++){
			for(j=0; j<NY; j++){
				for(n=0; n<NMW; n++){
					tmp[n][j][i]=(float)((V[n][j+1][i+1]+V[n][j+2][i+1])/(2*H[j+1][i+1]));
				}
			}
		}

		if ((status = nc_inq_varid(ncid, "v", &varid)))
			ERR(status);

		if ((status = nc_put_vara_float(ncid, varid, start, count, &tmp[0][0][0])))
			ERR(status);

	#endif // WRITE_VELOCITY


	#ifdef WRITE_PRESSURE

		// Write pressure
		for(i=0; i<NX; i++){
			for(j=0; j<NY; j++){
				for(n=0; n<NMW; n++){
					tmp[n][j][i]=(float)(p1[n][j+1][i+1]*phi_surf[n][j+1][i+1]/9.81);
				}
			}
		}

		if ((status = nc_inq_varid(ncid, "eta", &varid)))
			ERR(status);

		if ((status = nc_put_vara_float(ncid, varid, start, count, &tmp[0][0][0])))
			ERR(status);

	#endif // WRITE_PRESSURE


	#ifdef WRITE_WIND

		// Write TAUX
		for(n=0; n<NMW; n++){
			for(j=0; j<NY; j++){
				for(i=0; i<NX; i++){
					//tmp[n][j][i]=(float)(RHO*tau_x[j+1][i+1]*phi_surf[n][j+1][i+1]);
					tmp[n][j][i]=(float)(tau_x[j+1][i+1]);
				}
			}
		}

		if ((status = nc_inq_varid(ncid, "TAUX", &varid)))
			ERR(status);

		if ((status = nc_put_vara_float(ncid, varid, start, count, &tmp[0][0][0])))
			ERR(status);
			
		// Write TAUY
		for(n=0; n<NMW; n++){
			for(j=0; j<NY; j++){
				for(i=0; i<NX; i++){
					//tmp[n][j][i]=(float)(RHO*tau_y[j+1][i+1]*phi_surf[n][j+1][i+1]);
					tmp[n][j][i]=(float)(tau_y[j+1][i+1]);
				}
			}
		}

		if ((status = nc_inq_varid(ncid, "TAUY", &varid)))
			ERR(status);

		if ((status = nc_put_vara_float(ncid, varid, start, count, &tmp[0][0][0])))
			ERR(status);

	#endif // WRITE_PRESSURE

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
