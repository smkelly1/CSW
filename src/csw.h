// CSW header file
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>

// Input/output files
#define FILE_GRID  "../../22-12_grid/5th_deg_grid.nc"
#define FILE_TIDES "../TPXO.nc"
#define FILE_WIND  "../1987_winter_CSW.nc"
#define FILE_OUT   "out"
//#define FILE_R     "../r.nc"      // Spatially variable linear damping
//#define FILE_NU    "../nu.nc"     // Spatially variable horizontal viscosity
//#define FILE_KAPPA "../kappa.nc"     // Spatially variable horizontal diffusivity
//#define GAMMA 0.25                 // Mixing efficiency

// Grid spacing
#define DX ((1.0/5)*M_PI/180)    // Grid spacing in m or radians

// Grid size
#define NPX 1                     // Number of processors in X
#define NPY 2                     // Number of processors in X

//#define NX (3600/NPX)             // Grid size, this must be an integer
//#define NY (1460/NPY)             // This must be an integer. Note: reducing the total y-grid size will eliminate arctic cells
#define NX (1800/NPX)             // Grid size, this must be an integer
#define NY (730/NPY)              // This must be an integer. Note: reducing the total y-grid size will eliminate arctic cells

#define NM  1                     // Number of modes
#define NMW 1                     // Number of modes to write
#define NC  1                     // Number of tidal frequencies

// Wind
#define NXW  576                  // Wind file grid size, this must be an integer
#define NYW  361                  // This must be an integer
#define DT_F 3600                 // Wind forcing interval

// Time steps
#define DT   1200                  // Forward model time step [sec]
                                  // Approximate stable time steps:
                                  // 10th deg = 100 steps/period (dt=447 sec)
                                  // 25th deg = 200 (224 sec)
                                  // 50th deg = 400 (112 sec)
                                  // 100th deg = 800 (56 sec)
#define DT_W (3*3600)             // Pressure write time step
#define DT_D (32*3600)            // Diagnostics write time step
#define NT   (121*24*3600/DT)     // Simulation duration (time steps)

// Dissipation (commenting these parameters removes the relevant code)
//#define CD   0.0025               // Constant quadratic bottom drag (CD=0.0025 is standard)
//#define R    (1.0/(32*24*3600))   // Constant linear "Rayleigh" damping (or minimum value)
//#define NU   1000.0                // Constant horizontal viscosity (or minimum value). Bryan (1975) uses NU=u*DX/2 
                                  // Quick reference for U=1 cm/s: 
                                  // 1/10 deg = 50
                                  // 1/25 deg = 20
                                  // 1/50 deg = 10
                                  // 1/100 deg = 5
//#define KAPPA 0.0               // Constant horizontal diffusivity (or minimum value)
//#define R_MASK (1.0/(12*3600))  // Damping time-scale in low-wave resolution regions (applied in timestep_p.c)

// Set minimum depths for dynamics, forcing, and topographic coupling
#define H_MIN        100.0         // Minimum depth to solve equations
//#define H_MIN_FORCE  100.0         // Minimum depth to force internal tides (applied in read_tides.c)
#define H_MIN_COUPLE 100.0        // Minimum depth to couple modes (applied in read_grid.c)

// Flags for dynamics (These could be written to the input file, but it's quicker to compile than run MATLAB)
#define CORIOLIS                  // Include the Coriolis force
//#define MODECOUPLE                // Include topographic coupling 
#define SPHERE                    // Use spherical coordinates (Cartesian is default)
#define PERIODICBC                // Use periodic boundary conditions in longitude
#define NO_ANTARCTIC              // No forcing or topographic coupling south of 60 S (this should also be prefiltered in the wind file)

//#define TIDE_FORCING            // Use an Internal-Tide Generating Function
#define WIND_FORCING              // Use wind forcing

// Flags to write snapshots
//#define WRITE_VELOCITY          // Controls whether snapshots of velocity will be written
#define WRITE_PRESSURE            // Controls whether snapshots of pressure will be written
//#define WRITE_WIND              // Controls whether snapshots of wind stress will be written

// Flags to write time-averaged diagnostics
//#define ENERGY                    // Compute and write energy diagnostics
//#define FLUX                      // Compute and write energy diagnostics
//#define WORK                      // Compute and write energy diagnostics
//#define WRITE_SSH               // Compute and write the amplitude and phase of SSH
//#define WRITE_TRANSPORT         // Compute and write the amplitude and phase of transport

// Constants 
#define A   6371000.0             // radius of Earth
#define RHO 1000.0                // Reference density

// NetCDF stuff
#define ERRCODE 2
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}


////////////////////////////////////////////////////////////////////////
// Structure definitions
struct type_N {
	int x, y, m, t, c;
};

struct type_BOUNDARY {
	double p_Se[NM*NY], p_Re[NM*NY];
	double p_Sw[NM*NY], p_Rw[NM*NY];

	double p_Sn[NM*NX], p_Rn[NM*NX];
	double p_Ss[NM*NX], p_Rs[NM*NX];
	
	double v_Se[NM*(NY+1)], v_Re[NM*(NY+1)];
	double v_Sw[NM*(NY+1)], v_Rw[NM*(NY+1)];

	double u_Sn[NM*(NX+1)], u_Rn[NM*(NX+1)];
	double u_Ss[NM*(NX+1)], u_Rs[NM*(NX+1)];
	
	// Only need to send u east/west and v north/south for diffusion
	#ifdef AX
		double u_Se[NM*NY], u_Re[NM*NY];
		double u_Sw[NM*NY], u_Rw[NM*NY];

		double v_Sn[NM*NX], v_Rn[NM*NX];
		double v_Ss[NM*NX], v_Rs[NM*NX];
	#endif
};

struct type_ITGF {
	double omega[NC];
	float Ur[NC][NY][NX], Ui[NC][NY][NX];
	float Vr[NC][NY][NX], Vi[NC][NY][NX];
	double complex F[NC][NY][NX];
};

struct type_WIND {
	double lon[NXW+3], lat[NYW];
	float  tau_x[NYW][NXW+3], tau_y[NYW][NXW+3]; // Just assume you need the entire wind grid for each slab
	float  low_x[NYW][NXW+3], low_y[NYW][NXW+3]; 
};
	
	
////////////////////////////////////////////////////////////////////////
// Function prototypes
void read_tides(int);
void read_grid(int);
void read_wind(int, int);

void init_output(int);

void write_output(int, double, int);
void write_diagnostics(int, int, int);

void calc_forces(int);
void calc_divergence(void);

void calc_ITGF(double);

void pass_p(int);
void pass_uv(int);

void timestep_p(void);
void timestep_uv(void);


////////////////////////////////////////////////////////////////////////
// Declare variables (we're using global variables here, because the biguns wont fit on the stack. This is a bummer because the size has to be written in the header as a macro instead of read off the grid file during runtime)

// Temporary place to writing data 
float tmp[NMW][NY][NX];

// For trading boundary data via MPI
struct type_BOUNDARY BOUNDARY;

// Forcing 
struct type_ITGF ITGF; 
struct type_WIND WIND;
double tau_x[NY+2][NX+2];
double tau_y[NY+2][NX+2];

// Grid
double lon[NX+2];
double lat[NY+2];
double f[NY+2];
double H[NY+2][NX+2];
double c[NM][NY+2][NX+2];
double phi_bott[NM][NY+2][NX+2];
double phi_surf[NM][NY+2][NX+2];

#ifdef FILE_R
	double r[NM][NY+2][NX+2];
#endif

#ifdef FILE_NU
	double nu[NM][NY+2][NX+2];
#endif

#ifdef FILE_KAPPA
	double kappa[NM][NY+2][NX+2];
#endif

// Variables and forces
double U[NM][NY+2][NX+3];
double U1[NM][NY+2][NX+3];
double Fu[NM][NY][NX+1]; 
double Fu_1[NM][NY][NX+1];
double Fu_2[NM][NY][NX+1];
double Fu_eps[NM][NY][NX+1];
double dHdx_u[NY][NX+1];

double V[NM][NY+3][NX+2];
double V1[NM][NY+3][NX+2];
double Fv[NM][NY+1][NX]; 
double Fv_1[NM][NY+1][NX]; 
double Fv_2[NM][NY+1][NX];
double Fv_eps[NM][NY+1][NX]; 
double dHdy_v[NY+1][NX];

double p[NM][NY+2][NX+2]; 
double p1[NM][NY+2][NX+2];
double Fp[NM][NY][NX];
double Fp_1[NM][NY][NX];
double Fp_2[NM][NY][NX];
double dHdx[NY][NX];
double dHdy[NY][NX];

// Energy diagnostics and temporary storage for writing output
#ifdef WRITE_SSH
	float SSH_amp[NMW][NY][NX];
	float SSH_phase[NMW][NY][NX];
#endif

#ifdef WRITE_TRANSPORT
	float U_amp[NMW][NY][NX];
	float U_phase[NMW][NY][NX];
	float V_amp[NMW][NY][NX];
	float V_phase[NMW][NY][NX];
#endif

#ifdef ENERGY
	float KE[NMW][NY][NX];
	float PE[NMW][NY][NX];
#endif

#ifdef FLUX
	float up[NMW][NY][NX];
	float vp[NMW][NY][NX];
#endif

#ifdef WORK
	float W[NMW][NY][NX];
	float C0[NMW][NY][NX];
	float Cn[NMW][NY][NX];
	float D[NMW][NY][NX];
	float divF[NMW][NY][NX];
#endif
