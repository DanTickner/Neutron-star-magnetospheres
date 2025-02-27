/*
Header_Initial_Conditions_07.h

Specify the magnetic field strength, rotation rate, extent of radial coordinates and maximum basis function indices.
Doing this guarantees that values are in agreement between all files, and prevents lengthy-in-lines-of-code requirement to ensure maximum indices are equal.

V07: B_p_function for the case A(Psi) = -2 Psi was not accurate.
     Rotational velocity needed proper conversion to code units. Instead, let the user set the rotation period in seconds.
	 Rename Omega_ramp to rotation_ramp to move away from talking about Omega directly.
	 Rename max -> final when talking about rotation, since we are now decreasing the rotation period.
	 Rename n_timesteps to n_T for consistency with the coordinates.
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
	return 0;				// Case with A(Psi) = 0.
	//return -2.0 * sin( t );	// Case with A(Psi) = -2 Psi.
}


//----- Number of gridpoints and timesteps -----
int    n_points_r  = 30;
int    n_points_t  = 100;
int    n_T         = 150;
double delta_T     = 0.002;




//----- Parameters for rotation ramping -----
double P_SI_final = 2.0;	// Final rotion period in SI units (seconds).

int T_index_rotation_ramp_start = 50;
int T_index_rotation_ramp_stop  = 100;




//----- Maximum indices -----
int n_max   = 20;	// Maximum order of Chebyshev polynomials to use in approximation of radial functions and their derivatives.
int ell_max = 10;	// Maximum order of (associated) Legendre polynomials to use with VSH decomposition.




//----- Extent of radial coordinates (in code units) -----
double r_min = 1.0;
double r_max = 10.0;

bool set_r_max_to_2_R_LC = true;	// Automatically set r_max to twice the final value of R_LC.




//----- CSV parameters -----
// Choose output filename, including parent file structure. Don't include ".csv" extension because this will be modified to add _profiles and _history.
std::string output_filename  = "20240325_b_check_history_csv";
std::string log_file_comment = "Seems to work now. Only difference is that I'm running on my personal laptop. Run again to confirm. Also added options to not output each of the four CSVs but it won't affect the issue of _1_history being empty.";

bool output_csv_profiles   = true;	// Option to output at certain timesteps to a CSV containing radial and angular variation of properties of the NS.
bool output_csv_history    = true;	// Option to output at certain timesteps to a CSV containing global parameters of the NS.
bool output_csv_VSH_coeffs = false;	// Option to output at certain timesteps to a CSV containing values of the VSH coefficients and their radial derivatives.
bool output_csv_BCs        = false;	// Option to output at certain timesteps to a CSV containing the boundary values of the fields before and after application of the boundary conditions.

int csv_profiles_write_freq_r    = n_points_r/30;	// Write to profiles   CSV after this many steps in the radial direction.
int csv_profiles_write_freq_t    = n_points_t/30;	// Write to profiles   CSV after this many steps in the polar  direction.

int csv_profiles_write_freq_T    = 1;//500;				// Write to profiles   CSV after this many timesteps.
int csv_history_write_freq_T     = 1;				// Write to history    CSV after this many timesteps.
int csv_VSH_coeffs_write_freq_T  = 1000;			// Write to VSH coeffs CSV after this many timesteps.
int csv_BCs_write_freq_T         = 100;				// Write to BCs        CSV after this many timesteps.

int csv_profiles_write_i_min = 0;				// Only write to profiles CSV if i in this range. Useful for troubleshooting problematic radial intervals.
int csv_profiles_write_i_max = n_points_r;		// Only write to profiles CSV if i in this range. Useful for troubleshooting problematic radial intervals.
int csv_profiles_write_T_min = 0;			// Only write to profiles CSV if T_index in this range. Useful for troubleshooting problematic time intervals.
int csv_profiles_write_T_max = n_T;		// Only write to profiles CSV if T_index in this range. Useful for troubleshooting problematic time intervals.

int csv_precision = 15;							// Desired precision of CSV output.




//----- Screen output (cout) parameters -----
int cout_i      = 2;
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