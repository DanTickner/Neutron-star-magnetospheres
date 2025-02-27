/*
Header_Initial_Conditions_11.h

Specify the magnetic field strength, rotation rate, extent of radial coordinates and maximum basis function indices.
Doing this guarantees that values are in agreement between all files, and prevents lengthy-in-lines-of-code requirement to ensure maximum indices are equal.

V11: Add functionality for adding an angular twist.
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




//----- Extent of radial coordinates (in code units) -----
double      r_min                           = 1.0;
double      r_max                           = 10.0;




//----- Number of gridpoints and timesteps -----
int         n_r                            = 50;
int         n_t                            = 101;
int         T_index_max                    = 3000;
double      delta_T                        = 0.003;//0.0003;




//----- Parameters for rotation ramping and imparting of angular twist -----
double      P_SI_final                     = 0.0015;		// Final rotion period in SI units (seconds).

int         T_index_rotation_ramp_start    = 1e5+1;
int         T_index_rotation_ramp_stop     = 1e5+150;

int         T_index_twist_ramp_start       = 1;//1e5 + 151;
int         T_index_twist_ramp_stop        = 150;//1e5 + 200;

double      twist_final_Omega              = 0.15;	// Final angular velocity of twist relative to Omega_final.

double      twist_r_min                    = r_min;                  	// Radius r at which the twist starts, in NS radii.
double      twist_r_min_99pct              = r_min + 0.1*(r_max-r_min); // Radius r at which the twist reaches 99% of its maximum, in NS radii.
double      twist_r_max_99pct              = r_max - 0.1*(r_max-r_min); // Radius r at which the twist falls to 99% of its maximum, in NS radii.
double      twist_r_max                    = r_max;	                    // Radius r at which the twist ends, in NS radii.

double      twist_t_min                    = 30.0 * acos(-1.0) / 180.0;	// Angle theta at which the twist starts, in radians.
double      twist_t_min_99pct              = 34.0 * acos(-1.0) / 180.0; // Angle theta at which the twist reaches 99% of its maximum, in radians.
double      twist_t_max_99pct              = acos(-1.0) - twist_t_min_99pct;//40.0 * acos(-1.0) / 180.0; // Angle theta at which the twist falls to 99% of its maximum, in radians.
double      twist_t_max                    = acos(-1.0) - twist_t_min;//44.0 * acos(-1.0) / 180.0;	// Angle theta at which the twist ends, in radians.




//----- Maximum indices -----
int         ell_max                        = 20;			// Maximum order of (associated) Legendre polynomials to use with VSH decomposition.




//----- CSV parameters -----
// Choose output filename, including parent file structure. Don't include ".csv" extension because this will be modified to add _profiles and _history.
// Recommend using absolute filepaths, not relative filepaths. That is, "C:/Users/Dan/Documents/NS_code_offline_results/" and not "../../../../Dogyukfycuments/NS code offline results/";.
// This is because, if the filepath does not exist, no files will be saved but the code will still run as normal. It will be difficult to determine that the cause of no files is an incorrect filepath.
// Absolute filepaths are easier to check for errors than relative filepaths, especially if folders get moved or subfolders created after the fact. Further, in Windows they can be copied directly from the File Explorer.
std::string output_folder                  = "C:/Users/Dan/Documents/NS_code_offline_results/";			//"";		// Host folder for output CSV files. Include "/" at end. Set as empty string "" to use current folder.
std::string output_filename                = "20250220_a_twisted_axisymmetric";						  // Base filename to identify output CSVs.
std::string log_file_comment               = "Implement an axisymmetric twist becuase it ties in more nicely with what people have studied in the literature, allowing for direct comparisons to papers.";

bool        output_csv_profiles            = true;			// Option to output at certain timesteps to a CSV containing radial and angular variation of properties of the NS.
bool        output_csv_history             = true;			// Option to output at certain timesteps to a CSV containing global parameters of the NS.
bool        output_csv_VSH_coeffs          = false;			// Option to output at certain timesteps to a CSV containing values of the VSH coefficients and their radial derivatives.
bool        output_csv_BCs                 = false;			// Option to output at certain timesteps to a CSV containing the boundary values of the fields before and after application of the boundary conditions.

int         csv_profiles_write_freq_r      = 1;//n_r/50;		// Write to profiles CSV after this many steps in the radial direction.
int         csv_profiles_write_freq_t      = 1;//n_t/50;		// Write to profiles CSV after this many steps in the polar  direction.

int         csv_profiles_write_freq_T      = 50;			// Write after this many timesteps to CSV containing profiles.
int         csv_history_write_freq_T       = 1;				// Write after this many timesteps to CSV containing history.
int         csv_VSH_coeffs_write_freq_T    = 1000;			// Write after this many timesteps to CSV containing VSH coeffs.
int         csv_BCs_write_freq_T           = 50;			// Write after this many timesteps to CSV containing boundary values.

int         csv_profiles_write_i_min       = 0;				// Only write to profiles CSV if i in this range. Useful for troubleshooting problematic radial intervals.
int         csv_profiles_write_i_max       = n_r;			// Only write to profiles CSV if i in this range. Useful for troubleshooting problematic radial intervals.
int         csv_profiles_write_T_min       = 0;				// Only write to profiles CSV if T_index in this range. Useful for troubleshooting problematic time intervals.
int         csv_profiles_write_T_max       = T_index_max;	// Only write to profiles CSV if T_index in this range. Useful for troubleshooting problematic time intervals.




//----- Screen output (cout) parameters -----
int         cout_i                         = 20;
int         cout_j                         = 20;
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
double      friction_sigma_0               = 0.8;
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