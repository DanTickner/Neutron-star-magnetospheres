/*
Initial_conditions_03.h

Specify the magnetic field strength, rotation rate, extent of radial coordinates and maximum basis function indices.
Doing this guarantees that values are in agreement between all files, and prevents lengthy-in-lines-of-code requirement to ensure maximum indices are equal.

V03: \ell_max should be different depending whether there is rotation or no rotation.
*/

const double pi = acos( -1 );


//----- Number of gridpoints and timesteps -----
int n_points_r  = 1000;
int n_points_t  = 1000;
int n_timesteps = 70;

double delta_T = 0.008;

double T_max = n_timesteps * delta_T;


//----- Parameters for ramping function of Omega -----
double Omega_max = 0.2;

int T_index_Omega_ramp_start = 0;
int T_index_Omega_ramp_stop  = 800;

double dOmega_by_dT_index = Omega_max / ( T_index_Omega_ramp_stop - T_index_Omega_ramp_start + 1.0 );
double dOmega_by_dT       = dOmega_by_dT_index / delta_T;


//----- Maximum indices -----
int ell_max_rotation_off = 1;	// Maximum order of (associated) Legendre polynomials to use when rotation is switched off.
int ell_max_rotation_on  = 20;	// Maximum order of (associated) Legendre polynomials to use when rotation is switched on.


//----- Extent of radial coordinates (in code units) -----
double r_min    = 1.0;
double r_max    = 2.0 * 1.0 / Omega_max;


//----- Extent of polar coordinates -----
double t_min = 0;
double t_max = pi;


//----- CSV parameters -----
// Choose output filename, including parent file structure. Don't include ".csv" extension because this will be modified to add _profiles and _history.
std::string output_filename  = "20231009_e_ramp_immediately";
std::string log_file_comment = "No initial non-rotating period. Does it blow up after 5-8 timesteps now? If so, rotation is the cause of the blowup.";

int csv_profiles_write_freq_r    = n_points_r/40;	// Write to profiles   CSV after this many steps in the radial direction.
int csv_profiles_write_freq_t    = n_points_t/40;	// Write to profiles   CSV after this many steps in the polar  direction.
int csv_profiles_write_freq_T    = 1;				// Write to profiles   CSV after this many timesteps.
int csv_history_write_freq_T     = 1;				// Write to history    CSV after this many timesteps.
int csv_VSH_coeffs_write_freq_T  = 10;				// Write to VSH coeffs CSV after this many timesteps.

int csv_profiles_write_i_min = 0;				// Only write to profiles CSV if i in this range. Useful for troubleshooting problematic radial intervals.
int csv_profiles_write_i_max = n_points_r;		// Only write to profiles CSV if i in this range. Useful for troubleshooting problematic radial intervals.
int csv_profiles_write_T_min = 0;				// Only write to profiles CSV if T_index in this range. Useful for troubleshooting problematic time intervals.
int csv_profiles_write_T_max = n_timesteps;		// Only write to profiles CSV if T_index in this range. Useful for troubleshooting problematic time intervals.

int csv_precision = 15;							// Desired precision of CSV output.


//----- Screen output (cout) parameters -----
int cout_i      = 998;
int cout_j      = 5;
int cout_freq_T = 1;




//----- Outer sponge layer -----
/*
See Parfrey 2012, PhD thesis, \S3.9, Eqs. (3.39) and (3.40).
r_sponge is the radius at which we want the sponge layer to start.
Higher beta  = steeper rise, so focuses more on the outer boundary and less on the rest of the sponge layer.
Higher gamma = "flatter peak" so more of the outermost points experience more friction.
The trends in beta, gamma are guessed from playing with my Python code to plot the graph of friction( r ). PhD\20230911 Sponge layer\Plot_sponge_layer_exponential_01.py
*/

bool use_outer_sponge_layer = false;

double r_sponge         = 0.9 * r_max;
double friction_sigma_0 = 1.0;
double friction_gamma   = 6.0;
double friction_beta    = 4.0;


//----- Properties of the star -----
double R_star = 1.0e4;			// Radius of the star in metres.
double m_star = 1.0e32;			// Magnetic dipole moment of the star in A m^2 (10^3 emu).


//----- Constants -----
double c    = 299792458;			// Speed of light in m s^-1.
double mu_0 = 1.25663706212e-6;		// Permeability of free space in kg m s^-2 A^-2.


//----- Factors to convert between SI units and code units -----
double T_factor = R_star / c;
double B_factor = mu_0/(4.0*pi) * m_star * pow( R_star, -3 );
double E_factor = B_factor / c;


//----- Miscellaneous -----
int w = 14;			// fixed width for text output