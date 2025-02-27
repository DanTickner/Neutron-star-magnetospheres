/*
Header_Initial_Conditions_08.h

Specify the magnetic field strength, rotation rate, extent of radial coordinates and maximum basis function indices.
Doing this guarantees that values are in agreement between all files, and prevents lengthy-in-lines-of-code requirement to ensure maximum indices are equal.

V08: Rename n_points_r to n_r and n_points_t to n_t for brevity.
     Move pi to Header_Time_Evolution_38.h
	 Replace n_T by T_index_max.
*/




//----- Functions for initial magnetic field -----
// Non-rotating magnetic dipole.
// B = mu_0 m / ( 4pi ) ( 2cos(theta)/r^3 e_r + sin(theta)/r^3 e_theta ).
// We are working in code units, so the leading factors are set to 1, and r is normalised by the stellar radius.

double B_r_function( double r, double t ){
	return pow( r, -3.0 ) * 2.0 * cos( t );
}

double B_t_function( double r, double t ){
	return pow( r, -3.0 ) * sin( t );
}

double B_p_function( double r, double t ){
	return 0;				// Case with A(Psi) = 0.
	//return -2.0 * sin( t );	// Case with A(Psi) = -2 Psi.
}


//----- Number of gridpoints and timesteps -----
int         n_r                            = 500;//1000;
int         n_t                            = 500;//2000;
int         T_index_max                    = 1000;//5000;
double      delta_T                        = 0.0015;




//----- Option to start from a saved point in a previous run -----
std::string input_file_previous_run_base   = "20240509_a_test_csv_fields_from_scratch";
bool        load_fields_from_previous_run  = false;
int         T_index_initial                = 60;




//----- Parameters for rotation ramping -----
double      P_SI_final                     = 0.0015;		// Final rotion period in SI units (seconds).

int         T_index_rotation_ramp_start    = 50;
int         T_index_rotation_ramp_stop     = 2050;




//----- Maximum indices -----
int         n_max                           = 10;			// Maximum order of Chebyshev polynomials to use in approximation of radial functions and their derivatives.
int         ell_max                         = 10;			// Maximum order of (associated) Legendre polynomials to use with VSH decomposition.




//----- Extent of radial coordinates (in code units) -----
double      r_min                           = 1.0;
double      r_max                           = 10.0;
bool        set_r_max_to_2_R_LC             = false;			// Automatically set r_max to twice the final value of R_LC, overriding the choice of r_max in the line above.




//----- CSV parameters -----
// Choose output filename, including parent file structure. Don't include ".csv" extension because this will be modified to add _profiles and _history.
std::string output_filename                = "20240509_c_benchmark";
std::string log_file_comment               = "The system to start from a previous timestep is accurate to floating-point precision. Now make a benchmark file evolving to a certain point at which the code remains stable.";

bool        output_csv_fields              = true;			// Option to output at certain timesteps to a CSV containing the magnetic field values at each gridpoint.
bool        output_csv_profiles            = true;			// Option to output at certain timesteps to a CSV containing radial and angular variation of properties of the NS.
bool        output_csv_history             = true;			// Option to output at certain timesteps to a CSV containing global parameters of the NS.
bool        output_csv_VSH_coeffs          = false;			// Option to output at certain timesteps to a CSV containing values of the VSH coefficients and their radial derivatives.
bool        output_csv_BCs                 = false;			// Option to output at certain timesteps to a CSV containing the boundary values of the fields before and after application of the boundary conditions.

int         csv_profiles_write_freq_r      = n_r/30;		// Write to profiles CSV after this many steps in the radial direction.
int         csv_profiles_write_freq_t      = n_t/30;		// Write to profiles CSV after this many steps in the polar  direction.

int         csv_fields_write_freq_T        = 10;			// Write after this many timesteps to CSV containing magnetic and electric field values.
int         csv_profiles_write_freq_T      = 10;			// Write after this many timesteps to CSV containing profiles.
int         csv_history_write_freq_T       = 1;				// Write after this many timesteps to CSV containing history.
int         csv_VSH_coeffs_write_freq_T    = 1000;			// Write after this many timesteps to CSV containing VSH coeffs.
int         csv_BCs_write_freq_T           = 1000;			// Write after this many timesteps to CSV containing boundary values.

int         csv_fields_write_T_min         = 100;			// Only write to fields   CSV if T_index above this value. Saves filespace.
int         csv_profiles_write_i_min       = 0;				// Only write to profiles CSV if i in this range. Useful for troubleshooting problematic radial intervals.
int         csv_profiles_write_i_max       = n_r;			// Only write to profiles CSV if i in this range. Useful for troubleshooting problematic radial intervals.
int         csv_profiles_write_T_min       = 0;				// Only write to profiles CSV if T_index in this range. Useful for troubleshooting problematic time intervals.
int         csv_profiles_write_T_max       = T_index_max;	// Only write to profiles CSV if T_index in this range. Useful for troubleshooting problematic time intervals.




//----- Screen output (cout) parameters -----
int         cout_i                         = 2;
int         cout_j                         = 5;
int         cout_freq_T                    = 1;




//----- Outer sponge layer -----
/*
See Parfrey 2012, PhD thesis, \S3.9, Eqs. (3.39) and (3.40).
r_sponge is the radius at which we want the sponge layer to start.
Higher beta  = steeper rise, so focuses more on the outer boundary and less on the rest of the sponge layer.
Higher gamma = "flatter peak" so more of the outermost points experience more friction.
The trends in beta, gamma are guessed from playing with my Python code to plot the graph of friction( r ). PhD\20230911 Sponge layer\Plot_sponge_layer_exponential_01.py
*/

bool        use_outer_sponge_layer         = true;

double      r_sponge                       = 0.9 * r_max;
double      friction_sigma_0               = 1.0;
double      friction_gamma                 = 6.0;
double      friction_beta                  = 4.0;




//----- Properties of the star -----
double      R_star                         = 1.0e4;			// Radius of the star in metres.
double      m_star                         = 1.0e32;		// Magnetic dipole moment of the star in A m^2 (10^3 emu).




//----- Miscellaneous -----
int         w                              = 14;			// Fixed width for text output as the code runs.
int         w0                             = 20;			// Fixed width for text output of initial variables.

int        precision_csv                   = 15;			// Desired precision of CSV output. Tradeoff between precision and filesize.
int        precision_cout                  = 6;				// Desired precision when outputting evolved values to the screen.