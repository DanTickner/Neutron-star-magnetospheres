/*
Header_Initial_Conditions.h

Specify the magnetic field strength, rotation rate, extent of radial coordinates and maximum basis function indices.
Doing this guarantees that values are in agreement between all files, and prevents lengthy-in-lines-of-code requirement to ensure maximum indices are equal.

V05: Moving calculated quantities to Header_Time_Evolution_xx.h (from _32), making it clearer what the user can and cannot alter.
*/

const double pi = acos( -1 );


//----- Number of gridpoints and timesteps -----
int    n_points_r  = 300;
int    n_points_t  = 1000;
int    n_timesteps = 300;//6000;
double delta_T     = 0.002;


//----- Parameters for ramping function of Omega -----
double Omega_max = 0.095;

int T_index_Omega_ramp_start = 50;
int T_index_Omega_ramp_stop  = 800;


//----- Maximum indices -----
int ell_max = 10;	// Maximum order of (associated) Legendre polynomials to use with VSH decomposition.


//----- Extent of radial coordinates (in code units) -----
double r_min = 1.0;
double r_max = 10.0;//2.0 * R_LC_min;


//----- CSV parameters -----
// Choose output filename, including parent file structure. Don't include ".csv" extension because this will be modified to add _profiles and _history.
std::string output_filename  = "20240229_check_compiles_with_moved_parameters";
std::string log_file_comment = "Moving calculated quantities from Header_Initial_Conditions_xx.h (_05) to Header_Time_Evolution_xx.h (_32), making it clearer what the user can and cannot alter. Checking that the code still compiles when the ordering of variable declaration is changed.";

int csv_profiles_write_freq_r    = 10;//n_points_r/30;	// Write to profiles   CSV after this many steps in the radial direction.
int csv_profiles_write_freq_t    = 10;//n_points_t/30;	// Write to profiles   CSV after this many steps in the polar  direction.

int csv_profiles_write_freq_T    = 500;				// Write to profiles   CSV after this many timesteps.
int csv_history_write_freq_T     = 1;				// Write to history    CSV after this many timesteps.
int csv_VSH_coeffs_write_freq_T  = 1000;				// Write to VSH coeffs CSV after this many timesteps.
int csv_BCs_write_freq_T         = 100;				// Write to BCs        CSV after this many timesteps.

int csv_profiles_write_i_min = 0;				// Only write to profiles CSV if i in this range. Useful for troubleshooting problematic radial intervals.
int csv_profiles_write_i_max = n_points_r;		// Only write to profiles CSV if i in this range. Useful for troubleshooting problematic radial intervals.
int csv_profiles_write_T_min = 2000;				// Only write to profiles CSV if T_index in this range. Useful for troubleshooting problematic time intervals.
int csv_profiles_write_T_max = n_timesteps;		// Only write to profiles CSV if T_index in this range. Useful for troubleshooting problematic time intervals.

int csv_precision = 15;							// Desired precision of CSV output.


//----- Screen output (cout) parameters -----
int cout_i      = 290;
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

bool use_outer_sponge_layer = true;

double r_sponge         = 0.9 * r_max;
double friction_sigma_0 = 1.0;
double friction_gamma   = 6.0;
double friction_beta    = 4.0;


//----- Properties of the star -----
double R_star = 1.0e4;			// Radius of the star in metres.
double m_star = 1.0e32;			// Magnetic dipole moment of the star in A m^2 (10^3 emu).


//----- Miscellaneous -----
int w = 14;			// fixed width for text output