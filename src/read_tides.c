#include <math.h>
#include <netcdf.h>
#include <csw.h>

void read_tides(int rank)
{
	
	int i, j, k, x0, y0, ncid, varid, status;
	double dHdx, dHdy;
	
	// Define start points
	y0=(int)(floor(rank/NPX)*NY); // y start
	x0=(int)(floor(rank-floor(rank/NPX)*NPX)*NX); // x start
		
	size_t start_ITGF[]={0, y0, x0};
	size_t count_ITGF[]={NC, NY, NX};


	// Open the file
	if ((status = nc_open(FILE_TIDES, NC_NOWRITE, &ncid)))
		ERR(status);
		
			
	// Load surface tides				
	if ((status = nc_inq_varid(ncid, "omega", &varid)))
		ERR(status);

	if ((status = nc_get_var_double(ncid, varid, &ITGF.omega[0])))
		ERR(status);		
	

	if ((status = nc_inq_varid(ncid, "Ur", &varid)))
		ERR(status);

	if ((status = nc_get_vara_float(ncid, varid, start_ITGF, count_ITGF, &ITGF.Ur[0][0][0])))
		ERR(status);	
		
	if ((status = nc_inq_varid(ncid, "Ui", &varid)))
		ERR(status);

	if ((status = nc_get_vara_float(ncid, varid, start_ITGF, count_ITGF, &ITGF.Ui[0][0][0])))
		ERR(status);
			

	if ((status = nc_inq_varid(ncid, "Vr", &varid)))
		ERR(status);

	if ((status = nc_get_vara_float(ncid, varid, start_ITGF, count_ITGF, &ITGF.Vr[0][0][0])))
		ERR(status);	
		
	if ((status = nc_inq_varid(ncid, "Vi", &varid)))
		ERR(status);

	if ((status = nc_get_vara_float(ncid, varid, start_ITGF, count_ITGF, &ITGF.Vi[0][0][0])))
		ERR(status);
		
		
	// Compute forcing
	for(k=0; k<NC; k++){
		for(j=0; j<NY; j++){
			for(i=0; i<NX; i++){			
	
				if (H[j+1][i+1]>H_MIN) {	
									
					dHdx=(H[j+1][i+2]-H[j+1][i])/(2*A*cos(lat[j+1])*DX);
					dHdy=(H[j+2][i+1]-H[j][i+1])/(2*A*DX);
					
					ITGF.F[k][j][i]=(((double)(ITGF.Ur[k][j][i])+I*(double)(ITGF.Ui[k][j][i]))*dHdx+((double)(ITGF.Vr[k][j][i])+I*(double)(ITGF.Vi[k][j][i]))*dHdy)/H[j+1][i+1];
					
				}
				
			}
		}
	}
		
}
