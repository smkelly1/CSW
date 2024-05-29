// CSW header file

////////////////////////////////////////////////////////////////////////
// User configuration settings
////////////////////////////////////////////////////////////////////////

// Resolution
#define RES 25				        // Grid cells per degree
#define NM  1                       // Number of modes
#define NMW 1                       // Number of modes to write

// Processor layout
#define NPX 2                       // Number of processors in X
#define NPY 2                       // Number of processors in X

// Input/output files
#define FILE_GRID  "../25th_deg_grid.nc" // Grid file
//#define FILE_WIND  "../ONEDAY.nc" // Wind forcing file
#define FILE_OUT   "out"			// Snapshot output file
#define FILE_DIAG  "diag"			// Diagnostic average output file
#define FILE_TIDES "../TPXO.nc"     // Tidal forcing file
//#define FILE_R     "../r.nc"      // Spatially variable linear damping
//#define FILE_NU    "../nu.nc"     // Spatially variable horizontal viscosity
//#define FILE_KAPPA "../kappa.nc"  // Spatially variable horizontal diffusivity

// Time steps
#define DT          178.848         // Model time step [sec].  
									// Try: 5th deg = 1200 sec, 10th deg = 600 sec, 25th deg = 200                                   
#define DT_W        (12.42*3600)    // Snapshot time step
#define DT_D        (12.42*3600)    // Diagnostics averaging time step
#define NT          (30*12.42*3600/DT) // Simulation duration (time steps)

#define BETA        0.281105        // time step parameter for maximum stability (5/12 reduces to AB3) CSW uses AB3-AM4, the ROMS barotropic algorithm
#define GAMMA		0.088
#define EPSILON     0.013

// Dissipation (commenting these parameters removes the relevant code)
//#define CD          0.0025        // Constant quadratic bottom drag (CD=0.0025 is standard)
#define R          (1.0/(10*24*3600)) // Constant linear "Rayleigh" damping (or minimum value).
//#define RC2        (1.0/(32*24*3600)) // Constant linear "Rayleigh" damping (or minimum value) will be divided by c^2 for each mode.
//#define R_MAX      (1.0/(12*3600)) 
//#define NU          20.0            // Constant horizontal viscosity (or minimum value). Bryan (1975) uses NU=u*DX/2 
                                    // Quick reference for U=1 cm/s (but can use x10 bigger): 
                                    // 1/5  deg = 100 (for 100 km wavelength tau = 29 days)
                                    // 1/10 deg = 50  (tau = 59 days)
                                    // 1/25 deg = 20  (tau = 147 days)
//#define NU_MAX      2000.0
//#define KAPPA      	20.0            // Constant horizontal diffusivity (or minimum value)
//#define KAPPA_MAX   2000.0

// Controlling exponential growth
//#define DAMP_GROWTH                 // Compute the triple exponential running average to identify exponential growth, also write snapshots.
//#define NUM_PERIODS       2         // Number of inertial periods for the running average
//#define GROWTH_THRESHOLD  0.002     // SSH in m of bad signals (0.2 cm is a good start)
//#define GROWTH_MAX        0.01      // SSH in m of bad signals (1 cm is a good start)
//#define AMP_THRESHOLD     0.02      // SSH in m of bad signals (2 cm is a good start)
//#define AMP_MAX           0.04      // SSH in m of bad signals (4 cm is a good start)

// Set minimum depths for dynamics, forcing, and topographic coupling
#define H_MIN        16.0           // Minimum depth to solve equations
#define H_MIN_COUPLE 100.0          // Minimum depth to couple modes (applied in read_grid.c)
#define DH_MAX       0.25           // Maximum fractional change in depth between grid points (0.25 is good)

// Flags for dynamics 
#define CORIOLIS                    // Include the Coriolis force
//#define MODECOUPLE                  // Include topographic coupling 

// Flags to write Output
//#define WRITE_VELOCITY            // Write snapshots of velocity 
//#define WRITE_ETA                 // Write snapshots of modal SSH 
//#define WRITE_WIND                // Write snapshots of wind stress
//#define WRITE_MIX					// Write snapshots of mixed layer velocity

//#define ENERGY                    // Compute and write energy 
//#define FLUX                      // Compute and write energy flux
//#define WORK                      // Compute and write wind work, tidal generation, and scattering
#define WRITE_SSH                 // Compute and write the amplitude and phase of SSH (for tides)
//#define WRITE_TRANSPORT           // Compute and write the amplitude and phase of transport (for tides)


////////////////////////////////////////////////////////////////////////
// Stuff below here is probably okay to leave alone

// Grid spacing from RES
#if RES == 25
	#define DX ((1.0/25)*M_PI/180)  // Grid spacing in radians
	#define NX (9000/NPX)           // Grid size, this must be an integer
	#define NY (3650/NPY)           // This must be an integer. Note: reducing the total y-grid size will eliminate arctic cells
#elif RES == 10
	#define DX ((1.0/10)*M_PI/180)  // Grid spacing in radians
	#define NX (3600/NPX)           // Grid size, this must be an integer
	#define NY (1460/NPY)           // This must be an integer. Note: reducing the total y-grid size will eliminate arctic cells
#elif RES == 5
	#define DX ((1.0/5)*M_PI/180)   // Grid spacing in radians
	#define NX (1800/NPX)           // Grid size, this must be an integer
	#define NY (730/NPY)            // This must be an integer. Note: reducing the total y-grid size will eliminate arctic cells
#endif

// MERRA2 wind (don't edit unless you change wind products)
#ifdef FILE_WIND
	#define WIND_FORCING            // Use wind forcing
	#define NXW  576                // Wind file grid size, this must be an integer
	#define NYW  361                // This must be an integer
	#define DT_F 3600               // Wind forcing interval
#endif

// Tidal forcing 
#ifdef FILE_TIDES
	#define TIDE_FORCING            // Use an Internal-Tide Generating Function
	#define H_MIN_FORCE  100.0      // Minimum depth to force internal tides (applied in read_tides.c)
	#define NC  1                   // Number of tidal frequencies
#endif

// Constants 
#define A   6371000.0               // radius of Earth
#define RHO 1000.0                  // Reference density

// Define generic flags
#if defined(WRITE_ETA) || defined(WRITE_VELOCITY) || defined(WRITE_WIND) || defined(DAMP_GROWTH)
	#define WRITE_OUTPUT
#endif

#if defined(ENERGY) || defined(FLUX) || defined(WORK) || defined(WRITE_SSH) || defined(WRITE_TRANSPORT)
	#define WRITE_DIAGNOSTICS
	
#endif

// NetCDF stuff
#define ERRCODE 2
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}


////////////////////////////////////////////////////////////////////////
// Non-configuration stuff below here
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
// Standard C libraries
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>

////////////////////////////////////////////////////////////////////////
// Structure definitions
struct type_N {
	int x, y, m, t, c;
};

struct type_BOUNDARY {
	double u_Se[NM*NY], u_Re[NM*NY];
	double u_Sw[NM*NY], u_Rw[NM*NY];
	double u_Sn[NM*NX], u_Rn[NM*NX];
	double u_Ss[NM*NX], u_Rs[NM*NX];

	double v_Se[NM*NY], v_Re[NM*NY];
	double v_Sw[NM*NY], v_Rw[NM*NY];	
	double v_Sn[NM*NX], v_Rn[NM*NX];
	double v_Ss[NM*NX], v_Rs[NM*NX];

	double p_Se[NM*NY], p_Re[NM*NY];
	double p_Sw[NM*NY], p_Rw[NM*NY];
	double p_Sn[NM*NX], p_Rn[NM*NX];
	double p_Ss[NM*NX], p_Rs[NM*NX];
};

#ifdef TIDE_FORCING  
	struct type_ITGF {
		double omega[NC];
		float Ur[NC][NY][NX], Ui[NC][NY][NX];
		float Vr[NC][NY][NX], Vi[NC][NY][NX];
		double complex F[NC][NY][NX];
	};
#endif

#ifdef WIND_FORCING  
	struct type_WIND {
		double lon[NXW+3], lat[NYW];
		float  tau_x[NYW][NXW+3], tau_y[NYW][NXW+3]; // Assume you need the entire wind grid for every processor
	};
#endif
	
	
////////////////////////////////////////////////////////////////////////
// Function prototypes
void read_tides(int);
void read_grid(int);
void read_wind(int, int);

void write_output(int, int);
void write_diagnostics(int, int);

void calc_divergence(void);
void calc_forces(void);
void calc_diagnostics(void);

void pass_p(int);
void pass_uv(int);

////////////////////////////////////////////////////////////////////////
// Declare variables (we're using global variables here, because the biguns wont fit on the stack. This is a bummer because the size has to be written in the header as a macro instead of read off the grid file during runtime)

// Time
double t;

#ifdef WRITE_DIAGNOSTICS
	int Na;
#endif

// Temporary place to writing data 
float tmp[NMW][NY][NX];
float tmp2D[NY][NX];

// For trading boundary data via MPI
struct type_BOUNDARY BOUNDARY;

// Forcing 
#ifdef TIDE_FORCING  
	struct type_ITGF ITGF; 
#endif

#ifdef WIND_FORCING  
	struct type_WIND WIND;
	double tau_x[NY+2][NX+2];
	double tau_y[NY+2][NX+2];
#endif

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

#if defined(FILE_NU) || defined(DAMP_GROWTH) 
	double nu[NY+2][NX+2];
#endif

#if defined(FILE_KAPPA) || defined(DAMP_GROWTH) 
	double kappa[NY+2][NX+2];
#endif

// Variables and forces
double U[NM][NY+2][NX+2];
double UE[NM][NY+2][NX+2];
double Fu[NM][NY][NX]; 
double Fu_eps[NM][NY][NX];
double dHdx_u[NY][NX+1];
double U1[NM][NY+2][NX+2];
double U2[NM][NY+2][NX+2];
double U3[NM][NY+2][NX+2];	

double V[NM][NY+2][NX+2];
double VE[NM][NY+2][NX+2];
double Fv[NM][NY][NX]; 
double Fv_eps[NM][NY][NX]; 
double dHdy_v[NY+1][NX];
double V1[NM][NY+2][NX+2];
double V2[NM][NY+2][NX+2];
double V3[NM][NY+2][NX+2];	

double p[NM][NY+2][NX+2]; 
double pE[NM][NY+2][NX+2];
double Fp[NM][NY][NX];
double Fp_eps[NM][NY][NX];
double dHdx[NY][NX];
double dHdy[NY][NX];
double p1[NM][NY+2][NX+2];
double p2[NM][NY+2][NX+2];
double p3[NM][NY+2][NX+2];	

#ifdef DAMP_GROWTH
	double eta[NY+2][NX+2];
	double eta1[NY+2][NX+2];
	double eta2[NY+2][NX+2];
	double eta3[NY+2][NX+2];
	double etaG[NY+2][NX+2];
	double etaA[NY+2][NX+2];
	double flag_growth[NY+2][NX+2];
#endif

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
