/*
Header_Initial_Conditions_06.h

Specify the magnetic field strength, rotation rate, extent of radial coordinates and maximum basis function indices.
Doing this guarantees that values are in agreement between all files, and prevents lengthy-in-lines-of-code requirement to ensure maximum indices are equal.

V06: Moving initial magnetic field from dedicated header file to here, so that all user-controlled values are in one place.
*/

const double pi = acos( -1 );




//----- Functions for initial magnetic field -----
// Non-rotating magnetic dipole.
// B = mu_0 m / ( 4pi ) ( 2cos(theta)/r^3 e_r + sin(theta)/r^3 e_theta ).
// We are working in code units, so the leading factors are set to 1, and r is normalised by the stellar radius.

double B_r_function( double r, double t ){
	return pow( r, -3 ) * 2.0 * cos( t );
}

double B_t_function( double r, double t ){
	return pow( r, -3 ) * sin( t );
}

double B_p_function( double r, double t ){
	//return 0;
	double Psi = pow( sin(t), 2 ) / r;
	return -2.0 * Psi;
}


//----- Number of gridpoints and timesteps -----
int    n_points_r  = 300;
int    n_points_t  = 1000;
int    n_timesteps = 10;//6000;
double delta_T     = 0.002;




//----- Parameters for ramping function of Omega -----
double Omega_max = 0.095;

int T_index_Omega_ramp_start = 50;
int T_index_Omega_ramp_stop  = 800;




//----- Maximum indices -----
int n_max   = 20;	// Maximum order of Chebyshev polynomials to use in approximation of radial functions and their derivatives.
int ell_max = 10;	// Maximum order of (associated) Legendre polynomials to use with VSH decomposition.




//----- Extent of radial coordinates (in code units) -----
double r_min = 1.0;
double r_max = 10.0;//2.0 * R_LC_min;




//----- CSV parameters -----
// Choose output filename, including parent file structure. Don't include ".csv" extension because this will be modified to add _profiles and _history.
std::string output_filename  = "20240311_check_compiles_with_new_VSH_decomposition_quick";
std::string log_file_comment = "Renamed the VSH decomposition functions and now pass arguments to them. Check that the code still proceeds as normal. Faster version with less timesteps to check small tweaks.";

int csv_profiles_write_freq_r    = 10;//n_points_r/30;	// Write to profiles   CSV after this many steps in the radial direction.
int csv_profiles_write_freq_t    = 10;//n_points_t/30;	// Write to profiles   CSV after this many steps in the polar  direction.

int csv_profiles_write_freq_T    = 500;				// Write to profiles   CSV after this many timesteps.
int csv_history_write_freq_T     = 1;				// Write to history    CSV after this many timesteps.
int csv_VSH_coeffs_write_freq_T  = 1000;			// Write to VSH coeffs CSV after this many timesteps.
int csv_BCs_write_freq_T         = 100;				// Write to BCs        CSV after this many timesteps.

int csv_profiles_write_i_min = 0;				// Only write to profiles CSV if i in this range. Useful for troubleshooting problematic radial intervals.
int csv_profiles_write_i_max = n_points_r;		// Only write to profiles CSV if i in this range. Useful for troubleshooting problematic radial intervals.
int csv_profiles_write_T_min = 2000;			// Only write to profiles CSV if T_index in this range. Useful for troubleshooting problematic time intervals.
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