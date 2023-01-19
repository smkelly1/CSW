#include <math.h>
#include <netcdf.h>
#include "csw.h"

void read_grid(int rank)
{

	int x0, y0, ncid, varid, status;
	int i, j;

	#if defined(FILE_R) || defined(FILE_KAPPA) || defined(FILE_NU)
		int n;
	#endif

	// Define start points 
	y0=(int)(floor(rank/NPX)*NY); // y start
	x0=(int)(floor(rank-floor(rank/NPX)*NPX)*NX); // x start

	size_t start_H[]={y0, x0};
	size_t start_c[]={0, y0, x0};

	// Define counts to read 
	size_t count_H[]={NY+2, NX+2};
	size_t count_c[]={NM, NY+2, NX+2};

	// Define variables that might be needed
	#if defined(SPHERE) || defined(CORIOLIS)
		size_t start_lon[]={x0};
		size_t count_lon[]={NX+2};
		
		size_t start_lat[]={y0};
		size_t count_lat[]={NY+2};
	#endif

	// Open the grid file
	if((status = nc_open(FILE_GRID, NC_NOWRITE, &ncid)))
		ERR(status);

	// Read the grid
	if((status = nc_inq_varid(ncid, "H", &varid)))
		ERR(status);

	if((status = nc_get_vara_double(ncid, varid, start_H, count_H, &H[0][0])))
		ERR(status);

	if((status = nc_inq_varid(ncid, "c", &varid)))
		ERR(status);
		
	if((status = nc_get_vara_double(ncid, varid, start_c, count_c, &c[0][0][0])))
		ERR(status);

	if((status = nc_inq_varid(ncid, "phi_bott", &varid)))
		ERR(status);
	
	if((status = nc_get_vara_double(ncid, varid, start_c, count_c, &phi_bott[0][0][0])))
		ERR(status);
		
	if((status = nc_inq_varid(ncid, "phi_surf", &varid)))
		ERR(status);
	
	if((status = nc_get_vara_double(ncid, varid, start_c, count_c, &phi_surf[0][0][0])))
		ERR(status);

	// Read latitude in degrees and convert to radians
	#if defined(SPHERE) || defined(CORIOLIS)
	
		if((status = nc_inq_varid(ncid, "lon", &varid)))
			ERR(status);
	
		if((status = nc_get_vara_double(ncid, varid, start_lon, count_lon, &lon[0])))
			ERR(status);
	
		if((status = nc_inq_varid(ncid, "lat", &varid)))
			ERR(status);

		if((status = nc_get_vara_double(ncid, varid, start_lat, count_lat, &lat[0])))
			ERR(status);
			
		// Convert degrees to radians and calculate the inertial frequency
		for(i=0; i<NX+2; i++) {
			lon[i]=lon[i]*M_PI/180;
		}
		
		for(j=0; j<NY+2; j++){
			lat[j]=lat[j]*M_PI/180;
			f[j]=2*(2*M_PI)/(24*3600)*sin(lat[j]);	
		}
		
	#endif

	// Open the damping file
	#ifdef FILE_R
	if((status = nc_open(FILE_R, NC_NOWRITE, &ncid)))
		ERR(status);

	if((status = nc_inq_varid(ncid, "r", &varid)))
		ERR(status);
		
	if((status = nc_get_vara_double(ncid, varid, start_c, count_c, &r[0][0][0])))
		ERR(status);

		#ifdef NO_ANTARCTIC
			for(n=0; n<NM; n++){
				for(j=0; j<NY+2; j++){
					if(lat[j]<-60*M_PI/180) { // Recall: lat is in radians
						for(i=0; i<NX+2; i++){
							r[n][j][i]=0.0;
						}
					}
				}
			}
		#endif

		#ifdef R
			for(n=0; n<NM; n++){
				for(j=0; j<NY+2; j++){
					for(i=0; i<NX+2; i++){
						r[n][j][i]=fmax(r[n][j][i],R);
					}
				}
			}
		#endif
	#endif // end FILE_R

	// Open the viscosity file
	#ifdef FILE_NU
	if((status = nc_open(FILE_NU, NC_NOWRITE, &ncid)))
		ERR(status);
		
	if((status = nc_inq_varid(ncid, "nu", &varid)))
		ERR(status);
		
	if((status = nc_get_vara_double(ncid, varid, start_c, count_c, &nu[0][0][0])))
		ERR(status);

		#ifdef NO_ANTARCTIC
			for(n=0; n<NM; n++){
				for(j=0; j<NY+2; j++){
					if(lat[j]<-60*M_PI/180) { // Recall: lat is in radians
						for(i=0; i<NX+2; i++){
							nu[n][j][i]=0.0;
						}
					}
				}
			}
		#endif

		#ifdef NU
			for(n=0; n<NM; n++){
				for(j=0; j<NY+2; j++){
					for(i=0; i<NX+2; i++){
						nu[n][j][i]=fmax(nu[n][j][i],NU);
					}
				}
			}
		#endif
	#endif // end FILE_NU

	// Open the diffusivity file
	#ifdef FILE_KAPPA
	if((status = nc_open(FILE_KAPPA, NC_NOWRITE, &ncid)))
		ERR(status);
		
	if((status = nc_inq_varid(ncid, "kappa", &varid)))
		ERR(status);
		
	if((status = nc_get_vara_double(ncid, varid, start_c, count_c, &kappa[0][0][0])))
		ERR(status);

		#ifdef NO_ANTARCTIC
			for(n=0; n<NM; n++){
				for(j=0; j<NY+2; j++){
					if(lat[j]<-60*M_PI/180) { // Recall: lat is in radians
						for(i=0; i<NX+2; i++){
							kappa[n][j][i]=0.0;
						}
					}
				}
			}
		#endif

		#ifdef KAPPA
			for(n=0; n<NM; n++){
				for(j=0; j<NY+2; j++){
					for(i=0; i<NX+2; i++){
						kappa[n][j][i]=fmax(kappa[n][j][i],KAPPA);
					}
				}
			}
		#endif
	#endif // end FILE_KAPPA
	
	// Compute topographic gradients
	for(j=0; j<NY; j++){
		for(i=0; i<NX; i++){			
			dHdx[j][i]=(H[j+1][i+2]-H[j+1][i])/(2*A*cos(lat[j+1])*DX);
			dHdy[j][i]=(H[j+2][i+1]-H[j][i+1])/(2*A*DX);

			if((H[j+1][i+2]-H[j+1][i])/(2*H[j+1][i+1]) > 0.25){dHdx[j][i]=0.25*H[j+1][i+1]/(A*cos(lat[j+1])*DX);}
			if((H[j+1][i+2]-H[j+1][i])/(2*H[j+1][i+1]) < -0.25){dHdx[j][i]=-0.25*H[j+1][i+1]/(A*cos(lat[j+1])*DX);}

			if((H[j+2][i+1]-H[j][i+1])/(2*H[j+1][i+1]) > 0.25){dHdy[j][i]=0.25*H[j+1][i+1]/(A*DX);}
			if((H[j+2][i+1]-H[j][i+1])/(2*H[j+1][i+1]) < -0.25){dHdy[j][i]=-0.25*H[j+1][i+1]/(A*DX);}
		}
	}
	
	for(j=0; j<NY; j++){
		for(i=0; i<NX+1; i++){
			dHdx_u[j][i]=(H[j+1][i+1]-H[j+1][i])/(A*cos(lat[j+1])*DX);

			if((H[j+1][i+1]-H[j+1][i])/((H[j+1][i]+H[j+1][i+1])/2) > 0.25){dHdx_u[j][i]=0.25*((H[j+1][i]+H[j+1][i+1])/2)/(A*DX);}
			if((H[j+1][i+1]-H[j+1][i])/((H[j+1][i]+H[j+1][i+1])/2) < -0.25){dHdx_u[j][i]=-0.25*((H[j+1][i]+H[j+1][i+1])/2)/(A*DX);}
		}
	}
	
	for(j=0; j<NY+1; j++){
		for(i=0; i<NX; i++){
			dHdy_v[j][i]=(H[j+1][i+1]-H[j][i+1])/(A*DX);
			if((H[j+1][i+1]-H[j][i+1])/((H[j][i+1]+H[j+1][i+1])/2) > 0.25){dHdy_v[j][i]=0.25*((H[j][i+1]+H[j+1][i+1])/2)/(A*DX);}
			if((H[j+1][i+1]-H[j][i+1])/((H[j][i+1]+H[j+1][i+1])/2) < -0.25){dHdy_v[j][i]=-0.25*((H[j][i+1]+H[j+1][i+1])/2)/(A*DX);}		
		}
	}


}
