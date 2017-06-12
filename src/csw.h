// CSW header file
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>

// Input/output files
#define FILE_GRID  "../../17-6_global_grids/50th_deg_grid.nc"
#define FILE_IN    "../50th_deg_in.nc"
#define FILE_OUT   "50th_deg_out"

// Grid spacing
#define DX ((1.0/50)*M_PI/180) // Grid spacing in m or radians

// Grid size 
#define NPX 8  				// Number of processors in X
#define NPY 8				// Number of processors in X

#define NX (18000/NPX)       // Grid size, this must be an integer
#define NY (7296/NPY)       // This must be an integer 
							// Note: reducing the total y-grid size will eliminate arctic cells
#define NM 2			    // Number of modes
#define NC 1                // Number of tidal frequencies

// Time steps
#define DT   (12.42*3600/400)    // Forward model time step [sec]
#define DT_W (12.42*3600/2)      // Pressure write time step
#define DT_D (12.42*3600/1)      // Diagnostics write time step
#define NT   (100*12.42*3600/DT) // Simulation duration (time steps)

// Dissipation (commenting these parameters removes the relevant code)
#define R	(1.0/(10*12.42*3600)) // Linear "Rayleigh" damping
//#define CD	0.0025 	           // Quadratic bottom drag (CD=0.0025 is standard)
//#define AX	10.0                   // Horizontal Laplacian viscosity Bryan (1975) uses Ax=u*DX/2 
						 	       // Quick reference for U=1 cm/s: 1/10 deg = 50, 1/25 deg = 20, 1/50 deg = 10, 1/100 deg = 5

// Flags (These could be written to the input file, but it's quicker to compile than run MATLAB)
#define CORIOLIS               // Include the Coriolis force
#define MODECOUPLE             // Include topograhic coupling 
#define SPHERE                 // Use spherical coordinates (Cartesian is default)
#define PERIODICBC             // Use periodic boundary conditions in longitude
#define IT_FORCING             // Use an Internal-Tide Generating Function
//#define NO_FORCING_MASK        // Turn off internal-tide generation in masked regions (i.e., high latitude shelves)
#define DIAGNOSTICS            // Compute and write energy diagnostics

// Constants 
#define A   6371000.0          // radius of Earth
#define RHO 1000.0		       // Reference density

// NetCDF stuff
//#define WRITE_VELOCITY
#define ERRCODE 2
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}


////////////////////////////////////////////////////////////////////////
// Structure definitions
struct type_N {
	int x, y, m, t, c;
};

struct type_mask {
	double p[NM][NY][NX], U[NM][NY][NX+1], V[NM][NY+1][NX];
};

struct type_ITGF {
	double omega[NC];
	double pr[NC][NM][NY][NX], pi[NC][NM][NY][NX];
};

struct type_T {
	double x[NM][NM][NY+2][NX+2], y[NM][NM][NY+2][NX+2];
};

////////////////////////////////////////////////////////////////////////
// Function prototypes
void read_input(int);
void read_grid(int);

void init_output(int);

void write_output(int, double, int);
void write_diagnostics(int, int, int);
				
void calc_forces(void); 		 
void calc_divergence(void);
void calc_ITGF(double);

void pass_p(int);
void pass_uv(int);

void timestep_p(void);
void timestep_uv(void);


////////////////////////////////////////////////////////////////////////
// Declare variables (we're using global variables here, because the biguns wont fit on the stack. This is a bummer because the size has to be written in the header as a macro instead of read off the grid file during runtime)

// Input
struct type_mask mask; // Mask
struct type_ITGF ITGF; // Forcing

// Grid
double lat[NY+2];
double H[NY+2][NX+2];
double c[NM][NY+2][NX+2];
double f[NY+2][NX+2];
struct type_T T; // Coupling coefficients 

// Variables and forces
double U[NM][NY+2][NX+3];
double U1[NM][NY+2][NX+3];
double Fu[NM][NY][NX+1]; 
double Fu_1[NM][NY][NX+1];
double Fu_2[NM][NY][NX+1];
double Fu_eps[NM][NY][NX+1];

double V[NM][NY+3][NX+2];
double V1[NM][NY+3][NX+2];
double Fv[NM][NY+1][NX]; 
double Fv_1[NM][NY+1][NX]; 
double Fv_2[NM][NY+1][NX];
double Fv_eps[NM][NY+1][NX]; 
	
double p[NM][NY+2][NX+2]; 
double p1[NM][NY+2][NX+2];
double Fp[NM][NY][NX];
double Fp_1[NM][NY][NX];
double Fp_2[NM][NY][NX];

// Energy diagnostics and temporary storage for writing output
#ifdef DIAGNOSTICS
	float KE[NM][NY][NX];
	float PE[NM][NY][NX];
	float up[NM][NY][NX];
	float vp[NM][NY][NX];
	
	float C0[NM][NY][NX];
	float Cn[NM][NY][NX];
	float D[NM][NY][NX];
	float divF[NM][NY][NX];
	float residual[NM][NY][NX];
	
	float tmp[NM][NY][NX];
#endif
