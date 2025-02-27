/*
Header_Time_Evolution_36.h

Code snippets for the Time_Evolution_xx.cpp codes. Keeps them shorter and allows for two separate for-loops with Euler and Adams-Bashforth integration
without copying code.

V36: Rotational velocity needed proper conversion to code units. Instead, let the user set the rotation period in seconds.
     Move calculation of time and rotation values to a function void calculate_time_and_rotation_values(), instead of incrementing at each timestep, so that they can be checked.
	 Rename Omega_ramp to rotation_ramp to move away from talking about Omega directly.
	 Rename T_seconds to T_SI.
	 Rename n_timesteps to n_T for consistency with the coordinates.
	 In void output_parameters_to_screen() and void output_parameters_to_log_file(), change setw(16) to setw(w) so that width is customisable like the tabular output.
	 Move CFL max timestep calculation to function void calculate_CFL_max_timestep().
	 Re-introduce the factor 1/2 to magnetic_energy_density from my previous non-dimensionalisation.
	 Set output filenames to variables, e.g. std::string output_filename_profiles, for consistency if the naming convention changes.
	 Rename double T_factor = R_star / c to tau, for consistency with the writeup.
	 In void calculate_n_gridpoints_within_FF_region(), replace 1.0 / Omega[T_index] by R_LC[T_index] now that I know R_LC is calculated correctly.
*/

//----- Constants -----
double c    = 299792458.0;			// Speed of light in m s^-1.
double mu_0 = 1.25663706212e-6;		// Permeability of free space in kg m s^-2 A^-2.

//----- Factors to convert between SI units and code units -----
double tau      = R_star / c;
double B_factor = mu_0/(4.0*pi) * m_star * pow( R_star, -3 );
double E_factor = B_factor / c;

//----- Global variables -----
double T_max          = n_T * delta_T;
int    n_points_total = n_points_r * n_points_t;


double P_final            = P_SI_final / tau;			// Final rotation period  in code units.
double Omega_final        = 2.0 * pi / P_final;			// Final angular velocity in code units.
double dOmega_by_dT_index = Omega_final / ( T_index_rotation_ramp_stop - T_index_rotation_ramp_start + 1.0 );
double dOmega_by_dT       = dOmega_by_dT_index / delta_T;

double t_min       = 0;														// Minimum theta. Don't expect this to change from zero but keep as variable.
double t_max       = pi;													// Maximum theta. Don't expect this to change from pi   but keep as variable.
double delta_r     = ( r_max - r_min ) / ( (double) n_points_r - 1.0 );		// Linear grid spacing in radial direction.
double delta_t     = pi / ( (double) n_points_t - 1.0 );					// Linear grid spacing in polar  direction.
double delta_T_CFL = 0;														// Maximum timestep allowed by the CFL condition. Calculated by void calculate_CFL_max_timestep().


//--- Output files -----
std::ofstream output_file_profiles;
std::ofstream output_file_history;
std::ofstream output_file_VSH_coeffs;
std::ofstream output_file_BCs;
std::ofstream output_file_log;

std::string output_filename_profiles   = "CSV/"   + output_filename + "_1_profiles.csv";
std::string output_filename_history    = "CSV/"   + output_filename + "_2_history.csv";
std::string output_filename_VSH_coeffs = "CSV/"   + output_filename + "_3_VSH_coeffs.csv";
std::string output_filename_BCs        = "CSV/"   + output_filename + "_4_BCs.csv";
std::string output_filename_log        = "Logs/"  + output_filename + ".txt";


//----- Coordinates and arrays -----
std::vector<double> r  ( n_points_r );	// Radial coordinate r.
std::vector<double> R  ( n_points_r );	// Chebyshev integration variable R = cos^-1( lambda^-1( r ) ).
std::vector<double> t  ( n_points_t );	// Polar coordinate theta.
std::vector<double> st ( n_points_t );	// sin(theta).
std::vector<double> ct ( n_points_t );	// cos(theta).

std::vector<double> T                   ( n_T );	// Time in code units.
std::vector<double> T_SI                ( n_T );	// Time in SI   units.
std::vector<double> Omega               ( n_T );	// Angular velocity in code units.
std::vector<double> Omega_SI            ( n_T );	// Angular velocity in SI   units.
std::vector<double> P                   ( n_T );	// Rotation period in code units.
std::vector<double> P_SI                ( n_T );	// Rotation period in SI   units.
std::vector<double> R_LC                ( n_T );	// Light cylinder radius in code units.
std::vector<double> R_LC_SI             ( n_T );	// Light cylinder radius in SI units.
std::vector<double> star_rotation_angle ( n_T );	// Total angle through which the star has rotated, in radians.

std::vector< std::vector<double> > x ( n_points_r, std::vector<double> ( n_points_t ) );
std::vector< std::vector<double> > z ( n_points_r, std::vector<double> ( n_points_t ) );

std::vector< std::vector<double> > Tn ( n_points_r, std::vector<double> ( n_max + 1 ) );
std::vector< std::vector<double> > Un ( n_points_r, std::vector<double> ( n_max + 1 ) );
std::vector< std::vector<double> > P0 ( n_points_t, std::vector<double> ( ell_max+1 ) );
std::vector< std::vector<double> > P1 ( n_points_t, std::vector<double> ( ell_max+1 ) );

std::vector<double> sqrt_2Lplus1_over_4pi        ( ell_max+1 );
std::vector<double> sqrt_2Lplus1_pi              ( ell_max+1 );
std::vector<double> sqrt_2Lplus1_pi_over_L_Lplus1( ell_max+1 );

std::vector< std::vector< std::vector<double> > > B_VSH_coeffs     ( 3, std::vector< std::vector<double> > ( ell_max+1, std::vector<double> ( n_points_r ) ) );
std::vector< std::vector< std::vector<double> > > B_VSH_coeffs_dr  ( 3, std::vector< std::vector<double> > ( ell_max+1, std::vector<double> ( n_points_r ) ) );
std::vector< std::vector< std::vector<double> > > B_VSH_coeffs_dr2 ( 3, std::vector< std::vector<double> > ( ell_max+1, std::vector<double> ( n_points_r ) ) );
std::vector< std::vector< std::vector<double> > > E_VSH_coeffs     ( 3, std::vector< std::vector<double> > ( ell_max+1, std::vector<double> ( n_points_r ) ) );
std::vector< std::vector< std::vector<double> > > E_VSH_coeffs_dr  ( 3, std::vector< std::vector<double> > ( ell_max+1, std::vector<double> ( n_points_r ) ) );
std::vector< std::vector< std::vector<double> > > E_VSH_coeffs_dr2 ( 3, std::vector< std::vector<double> > ( ell_max+1, std::vector<double> ( n_points_r ) ) );

std::vector< std::vector< std::vector<double> > > B                 ( 3, std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) ) );
std::vector< std::vector< std::vector<double> > > B_minus1          ( 3, std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) ) );
std::vector< std::vector< std::vector<double> > > B_minus2          ( 3, std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) ) );
std::vector< std::vector< std::vector<double> > > B_minus3          ( 3, std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) ) );
std::vector< std::vector< std::vector<double> > > B_dT              ( 3, std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) ) );
std::vector< std::vector< std::vector<double> > > B_dT_minus1       ( 3, std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) ) );
std::vector< std::vector< std::vector<double> > > B_dT_minus2       ( 3, std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) ) );
std::vector< std::vector< std::vector<double> > > B_dT_minus3       ( 3, std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) ) );
std::vector< std::vector< std::vector<double> > > B_curl            ( 3, std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) ) );
std::vector< std::vector< std::vector<double> > > B_Cartesian       ( 3, std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) ) );
std::vector< std::vector< std::vector<double> > > B_dr              ( 3, std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) ) );
std::vector< std::vector< std::vector<double> > > B_dr2             ( 3, std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) ) );
std::vector< std::vector< std::vector<double> > > B_series          ( 3, std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) ) );	// Values of the VSH decomposition of the field components.
std::vector< std::vector< std::vector<double> > > B_abs_rel_dev     ( 3, std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) ) );
std::vector< std::vector< std::vector<double> > > B_from_VSH_coeffs ( 3, std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) ) );

std::vector< std::vector< std::vector<double> > > E                 ( 3, std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) ) );
std::vector< std::vector< std::vector<double> > > E_minus1          ( 3, std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) ) );
std::vector< std::vector< std::vector<double> > > E_minus2          ( 3, std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) ) );
std::vector< std::vector< std::vector<double> > > E_minus3          ( 3, std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) ) );
std::vector< std::vector< std::vector<double> > > E_dT              ( 3, std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) ) );
std::vector< std::vector< std::vector<double> > > E_dT_minus1       ( 3, std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) ) );
std::vector< std::vector< std::vector<double> > > E_dT_minus2       ( 3, std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) ) );
std::vector< std::vector< std::vector<double> > > E_dT_minus3       ( 3, std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) ) );
std::vector< std::vector< std::vector<double> > > E_curl            ( 3, std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) ) );
std::vector< std::vector< std::vector<double> > > E_Cartesian       ( 3, std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) ) );
std::vector< std::vector< std::vector<double> > > E_dr              ( 3, std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) ) );
std::vector< std::vector< std::vector<double> > > E_dr2             ( 3, std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) ) );
std::vector< std::vector< std::vector<double> > > E_series          ( 3, std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) ) );	// Values of the VSH decomposition of the field components.
std::vector< std::vector< std::vector<double> > > E_abs_rel_dev     ( 3, std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) ) );
std::vector< std::vector< std::vector<double> > > E_from_VSH_coeffs ( 3, std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) ) );

std::vector< std::vector< std::vector<double> > > J                 ( 3, std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) ) );	// Current density vector.
std::vector< std::vector< std::vector<double> > > J_Cartesian       ( 3, std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) ) );


std::vector< std::vector<double> > B_squared                     ( n_points_r, std::vector<double> ( n_points_t ) );
std::vector< std::vector<double> > E_squared                     ( n_points_r, std::vector<double> ( n_points_t ) );
std::vector< std::vector<double> > B_squared_minus_E_squared     ( n_points_r, std::vector<double> ( n_points_t ) );
std::vector< std::vector<double> > E_dot_B                       ( n_points_r, std::vector<double> ( n_points_t ) );
std::vector< std::vector<double> > E_div                         ( n_points_r, std::vector<double> ( n_points_t ) );
std::vector< std::vector<double> > magnetic_energy_density       ( n_points_r, std::vector<double> ( n_points_t ) );		// Energy density
std::vector< std::vector<double> > friction                      ( n_points_r, std::vector<double> ( n_points_t ) );		// Sponge layer at outer boundary to absorb outgoing waves (see Parfrey, 2012, PhD thesis, S3.9).
std::vector< std::vector<int   > > is_gridpoint_within_FF_region ( n_points_r, std::vector<int   > ( n_points_t ) );
std::vector< std::vector<int   > > were_FF_conditions_applied    ( n_points_r, std::vector<int   > ( n_points_t ) );		// +1: FF conditions applied here. -1: FF conditions not applied here. 0: Outside FF region.

std::vector<double> finite_line_element_dR ( n_points_r );	// Finite line elements dR for integrals to calculate Chebyshev series coefficients.
std::vector<double> finite_line_element_dt ( n_points_t );	// Finite line elements dtheta for integral over theta.
std::vector< std::vector<double> > finite_area_element_for_r_t ( n_points_r, std::vector<double> ( n_points_t ) );	// Finite area elements r dr dtheta for 2D integral over r and theta.

std::vector<double> B_stdev          ( 3 );	// Standard devation of the VSH decomposition of each vector component.
std::vector<double> B_max_abs_rel_dev( 3 );	// Maximum absolute relative deviation between the VSH decomposition and the exact value of each vector component.
std::vector<double> E_stdev          ( 3 );	// Standard devation of the VSH decomposition of each vector component.
std::vector<double> E_max_abs_rel_dev( 3 );	// Maximum absolute relative deviation between the VSH decomposition and the exact value of each vector component.

std::vector<double> B_t_dt_inner (n_points_t); // 20231214 Specify d/dT B at inner boundary; requires expression for d/dt B_t.

double total_magnetic_energy = 0;

std::vector<int> B_nans ( 6 );	// Count the number of "nan" or "inf" values of each field component at a given timestep.
std::vector<int> E_nans ( 6 );	// Count the number of "nan" or "inf" values of each field component at a given timestep.
int nans_tot = 0;

std::vector<double> VSH_integrand_r ( n_points_t );
std::vector<double> VSH_integrand_1 ( n_points_t );
std::vector<double> VSH_integrand_2 ( n_points_t );

double time_start_seconds = std::chrono::high_resolution_clock::now().time_since_epoch().count() * 1e-9;
time_t time_start_integer = time( 0 );
char*  time_start_string  = ctime( &time_start_integer );

int    T_index = 0;

int    n_gridpoints_within_FF_region                              = 0;
int    num_E_points_changed_for_second_FF_condition               = 0;	// Output to CSV to keep track of how active the FF condition application is.
double percentage_E_points_points_changed_for_second_FF_condition = 0;

bool   is_rotation_on            = false;	// Used by functions with different output depending whether rotation is on or off.

std::vector< std::vector<double> > fields_inner_before_BCs ( 6, std::vector<double> ( n_points_t ) );
std::vector< std::vector<double> > fields_inner_after_BCs  ( 6, std::vector<double> ( n_points_t ) );
std::vector< std::vector<double> > fields_outer_before_BCs ( 6, std::vector<double> ( n_points_t ) );
std::vector< std::vector<double> > fields_outer_after_BCs  ( 6, std::vector<double> ( n_points_t ) );

std::vector< std::vector<double> > B_inner_before_BCs ( 3, std::vector<double> ( n_points_t ) );
std::vector< std::vector<double> > B_inner_after_BCs  ( 3, std::vector<double> ( n_points_t ) );
std::vector< std::vector<double> > B_outer_before_BCs ( 3, std::vector<double> ( n_points_t ) );
std::vector< std::vector<double> > B_outer_after_BCs  ( 3, std::vector<double> ( n_points_t ) );
std::vector< std::vector<double> > E_inner_before_BCs ( 3, std::vector<double> ( n_points_t ) );
std::vector< std::vector<double> > E_inner_after_BCs  ( 3, std::vector<double> ( n_points_t ) );
std::vector< std::vector<double> > E_outer_before_BCs ( 3, std::vector<double> ( n_points_t ) );
std::vector< std::vector<double> > E_outer_after_BCs  ( 3, std::vector<double> ( n_points_t ) );




//----- Functions -----

double lambda( double x, double A, double B ){
	// Linear transformation from x in [-1,1] to X in [A,B] such that X = lambda( x, A, B ).
	return 0.5 * ( B - A ) * x + 0.5 * ( A + B );
}




double lambda_inverse( double X, double A, double B ){
	// Linear transformation from X in [A,B] to x in [-1,1] such that X = lambda_inverse( x, A, B ).
	return ( 2.0 * X - B - A ) / ( B - A );
}




void calculate_gridpoints(){
	
	//--- Radial and polar coordinates, plus Chebyshev integration variable R, sin(theta) and cos(theta) ---
	for( int i=0; i<n_points_r; i++ ){
		r[i] = r_min + i * delta_r;
		R[i] = acos( lambda_inverse( r[i], r_min, r_max ) );
	}
	for( int j=0; j<n_points_t; j++ ){
		t [j] = t_min + j * delta_t;
		st[j] = sin( t[j] );
		ct[j] = cos( t[j] );
	}
	
	//--- Cartesian coordinates ---
	for( int i=0; i<n_points_r; i++ ){
		for( int j=0; j<n_points_t; j++ ){
			x[i][j] = r[i] * st[j];
			z[i][j] = r[i] * ct[j];
		}
	}
	
	//--- Finite line and area elements for numerical integration ---
	for( int i=0; i<=n_points_r-2; i++ ){
		finite_line_element_dR[i] = R[i+1] - R[i];
	}
	
	for( int j=0; j<=n_points_t-2; j++ ){
		finite_line_element_dt[j] = ct[j] - ct[j+1];
	}
	
	for( int i=0; i<=n_points_r-2; i++ ){
		for( int j=0; j<=n_points_t-2; j++ ){
			finite_area_element_for_r_t[i][j] = ( r[i+1] - r[i] ) * ( t[j+1] - t[j] ) * ( r[i] + 0.5 * delta_r );
		}
	}
	
}




void calculate_remapped_chebyshev_polynomials(){
	
	for( int i=0; i<n_points_r; i++ ){
		
		double lambda_inverse_r_i = lambda_inverse( r[i], r_min, r_max ); 	// Saves repeating this calculation for all n.
		
		Tn[i][0] = 1.0;
		Tn[i][1] = lambda_inverse_r_i;
		Un[i][0] = 1.0;
		Un[i][1] = 2.0 * lambda_inverse_r_i;
		
		for( int n=2; n<=n_max; n++ ){
			Tn[i][n] = 2.0 * lambda_inverse_r_i * Tn[i][n-1] - Tn[i][n-2];
			Un[i][n] = 2.0 * lambda_inverse_r_i * Un[i][n-1] - Un[i][n-2];
		}
		
	}
	
}




void calculate_associated_legendre_functions(){
	
	for( int j=0; j<n_points_t; j++ ){
		
		//--- m = 0 ---
		P0[j][0] = 1.0;
		P0[j][1] = ct[j];
		
		for( int ell=2; ell<=ell_max; ell++ ){
			P0[j][ell] = ( (2.0*ell-1.0) * ct[j] * P0[j][ell-1] - (ell-1.0) * P0[j][ell-2] ) / ( (double) ell );
		}
		
		//--- m = 1 ---
		P1[j][1] = - st[j];
		P1[j][2] = - 3.0 * ct[j] * st[j];
		
		for( int ell=3; ell<=ell_max; ell++ ){
			P1[j][ell] = ( (2.0*ell-1.0) * ct[j] * P1[j][ell-1] - (double) ell * P1[j][ell-2] ) / ( (double) ell - 1.0 );
		}
		
	}
	
}




void calculate_sqrt_2Lplus1_over_4pi(){
	for( int ell=0; ell<=ell_max; ell++ ){
		sqrt_2Lplus1_over_4pi        [ell] = sqrt( (2.0*ell+1.0) / ( 4.0*pi ) );
		sqrt_2Lplus1_pi              [ell] = sqrt( (2.0*ell+1.0) * pi );
		sqrt_2Lplus1_pi_over_L_Lplus1[ell] = sqrt( (2.0*ell+1.0) * pi ) / ( (double)ell*(ell+1.0) );
	}
}




void calculate_time_and_rotation_values(){
	
	//--- Time in code units ---
	for( int T_index=0; T_index<n_T; T_index++ ){	
		T[T_index] = (double) T_index * delta_T;
	}
	
	//--- Angular velocity in code units ---
	Omega[ T_index_rotation_ramp_start ] = dOmega_by_dT_index;
	
	for( int T_index = T_index_rotation_ramp_start + 1; T_index <= T_index_rotation_ramp_stop; T_index++ ){
		Omega[T_index] = Omega[ T_index - 1 ] + dOmega_by_dT_index;
	}
	
	for( int T_index = T_index_rotation_ramp_stop + 1; T_index < n_T; T_index++ ){
		Omega[T_index] = Omega[ T_index - 1 ];
	}
	
	//--- Rotation period and light cylinder radius ---
	for( int T_index=0; T_index<n_T; T_index++ ){
		if( Omega[T_index] != 0 ){
			P   [T_index] = 2.0 * pi / Omega[T_index];
			R_LC[T_index] = 1.0      / Omega[T_index];
		}
	}
	
	//--- Total angle through which star has rotated in radians ---
	for( int T_index=T_index_rotation_ramp_start; T_index<n_T; T_index++ ){
		star_rotation_angle[T_index] = star_rotation_angle[ T_index - 1 ] + Omega[T_index] * delta_T;
	}
	
	//--- SI units ---
	for( int T_index=0; T_index<n_T; T_index++ ){
		T_SI    [T_index] = T    [T_index] * tau;
		Omega_SI[T_index] = Omega[T_index] / tau;
		P_SI    [T_index] = P    [T_index] * tau;
		R_LC_SI [T_index] = R_LC [T_index] * R_star;
	}
	
}




void calculate_CFL_max_timestep(){
	// This is the approximation for constant grid spacing in spherical coordinates described in the writeup.
	// Replace by an exact expression if needed, possibly using functions from Test codes/Test_CFL_condition.cpp.
	delta_T_CFL = std::min( r_min * delta_t, delta_r );
}




void calculate_friction_for_sponge_layer(){
	// Parfrey, 2012, PhD thesis, \S3.9.2.
	
	if( use_outer_sponge_layer ){
		for( int i=0; i<n_points_r; i++ ){
			for( int j=0; j<n_points_t; j++ ){
				friction[i][j] = pow( ( r[i] - r_sponge ) / ( r_max - r_sponge ), friction_beta );
				friction[i][j] = friction_sigma_0 * ( 1.0 - exp( - friction_gamma * friction[i][j] ) );
			}
		}
	}
	
}




double integral_trapezium_1D( std::vector<double> f, std::vector<double> delta_u ){
	// 1D integral of f(theta) in polar coordinates, given pre-calculated list of function values.
	
	double ret = 0;
	
	for( int j=0; j<=f.size()-2; j++ ){
		ret += ( f[j+1] + f[j] ) * delta_u[j];
	}
	
	return 0.5 * ret;

}




double integral_trapezium_2D( std::vector< std::vector<double> > f, std::vector< std::vector<double> > delta_A ){
	// 2D integral of f(r,theta) in circular/spherical polar coordinates, given pre-calculated 2D array of function values.
	
	double ret = 0;
	
	for( int i=0; i<=f.size()-2; i++ ){
		for( int j=0; j<=f[0].size()-2; j++ ){
			ret += ( f[i][j] + f[i][j+1] + f[i+1][j] + f[i+1][j+1] ) * delta_A[i][j];
		}
	}
	
	return 0.25 * ret;
	
}




void apply_initial_field_values_for_B(){
	for( int i=0; i<n_points_r; i++ ){
		for( int j=0; j<n_points_t; j++ ){
			B[0][i][j] = B_r_function( r[i], t[j] );
			B[1][i][j] = B_t_function( r[i], t[j] );
			B[2][i][j] = B_p_function( r[i], t[j] );
		}
	}
}




void apply_inner_boundary_conditions(){
	// Petri 2012, \S3.3.
	
	//--- Record values before inner BCs applied ---
	for( int j=0; j<n_points_t; j++ ){
		for( int f=0; f<3; f++ ){
			B_inner_before_BCs[f][j] = B[f][0][j];
			E_inner_before_BCs[f][j] = E[f][0][j];
		}
	}
	
	//--- Apply inner BCs ---
	for( int j=0; j<n_points_t; j++ ){
		B[0][0][j] = B_r_function( r_min, t[j] );
		//B[2][0][j] = 2.0 * delta_T * Omega * st[j] * cos(t[j]);					// SPECIFIC TO A DIPOLE AND NOT SURE IT'S THEORETICALLY GROUNDED; JUST SEEMS TO FIT THE DATA.
		E[1][0][j] = B_r_function( r_min, t[j] ) * - Omega[T_index] * r_min * st[j];
		E[2][0][j] = 0;
	}
	
	//B[2][0][ (int) n_points_t / 2 - 1 ] = 0;	// Set B_phi ( R_NS, pi/2 ) = 0 to avoid evolution of noise. Suggestion from Sam 20231208.
	
	//--- Record values after inner BCs applied ---
	for( int j=0; j<n_points_t; j++ ){
		for( int f=0; f<3; f++ ){
			B_inner_after_BCs[f][j] = B[f][0][j];
			E_inner_after_BCs[f][j] = E[f][0][j];
		}
	}

}




void apply_outer_boundary_conditions(){
	// Petri 2012, \S3.3.
	
	//--- Record values before outer BCs applied ---
	for( int j=0; j<n_points_t; j++ ){
		for( int f=0; f<3; f++ ){
			B_outer_before_BCs[f][j] = B[f].back()[j];
			E_outer_before_BCs[f][j] = E[f].back()[j];
		}
	}
	
	//--- Apply outer BCs ---
	if( T_index < T_index_rotation_ramp_start ){
		for( int j=0; j<n_points_t; j++ ){
			B[1].back()[j] = B_t_function( r_max, t[j] );
			B[2].back()[j] = B_p_function( r_max, t[j] );
			E[1].back()[j] = 0;
			E[2].back()[j] = 0;
		}
	}
	else {
		for( int j=0; j<n_points_t; j++ ){
			B[1].back()[j] = 0.5 * ( B[1].back()[j] - E[2].back()[j] );
			B[2].back()[j] = 0.5 * ( B[2].back()[j] + E[1].back()[j] );
			E[1].back()[j] =   B[2].back()[j];
			E[2].back()[j] = - B[1].back()[j];
		}
	}
	
	//--- Record values after outer BCs applied ---
	for( int j=0; j<n_points_t; j++ ){
		for( int f=0; f<3; f++ ){
			B_outer_after_BCs[f][j] = B[f].back()[j];
			E_outer_after_BCs[f][j] = E[f].back()[j];
		}
	}
	
}




std::vector<double> chebyshev_series_coeffs( std::vector<double> f ){
	// Calculate the coefficients of the Chebyshev series of a discrete list of function values f on the interval [r_min,r_max]. NOT WORKING 20240304.
	
	std::vector<double> coeffs    ( n_max + 1  );
	std::vector<double> integrand ( n_points_r );
	
	
	//--- Loop through indices n and calculate coefficients ---
	for( int n=0; n<=n_max; n++ ){
		
		//--- Build list of integrand values ---
		for( int i=0; i<n_points_r; i++ ){
			integrand[i] = f[i] * cos( (double) n * R[i] );
		}
		
		//--- Calculate coefficient by trapezium rule and multiply by constants ---
		coeffs[n] = integral_trapezium_1D( integrand, finite_line_element_dR );
		
		double h_n = 1.0;
		if( n == 0 ){
			h_n = 0.5;
		}
		
		coeffs[n] *= - h_n * 2.0 / pi;
		
	}
	
	return coeffs;
	
}




double evaluate_chebyshev_series_at_gridpoint( std::vector<double> coeffs, int i ){
	// Evaluate the Chebyshev series of a function on the interval [r_min,r_max] for the case where we are already on a gridpoint i.
	// This is computationally faster than for a general point because we have already calculated the Chebyshev polynomials at the gridpoints.
	
	double ret = 0;
	
	for( int n=0; n<=n_max; n++ ){
		ret += coeffs[n] * Tn[i][n];
	}
	
	return ret;

}




double evaluate_chebyshev_series_first_derivative_at_gridpoint( std::vector<double> coeffs, int i ){
	// Evaluate the Chebyshev series of the first derivative of a function on the interval [r_min,r_max] for the case where we are already on a gridpoint i.
	// This is computationally faster than for a general point because we have already calculated the Chebyshev polynomials at the gridpoints.
	
	double ret = 0;
	
	for( int n=1; n<=n_max; n++ ){
		ret += coeffs[n] * n * Un[i][n-1] * 2.0 / ( r_max - r_min );
	}
	
	return ret;

}




void VSH_decomposition( std::vector< std::vector< std::vector<double> > > &v, std::vector< std::vector< std::vector<double> > > &v_VSH_coeffs ){
	// Decomposition for vectors which in general are NOT divergence-free.
	// First, built lists of integrand values. Then, numerically integrate and put the values into an array.
	// Two separate loops due to two separate ranges on ell.
	
	for( int i=0; i<n_points_r; i++ ){
		
		//--- r-coefficients go from ell=0 ---
		for( int ell=0; ell<=ell_max; ell++ ){
			
			for( int j=0; j<n_points_t; j++ ){
				VSH_integrand_r[j] = v[0][i][j] * P0[j][ell];
			}
			
			v_VSH_coeffs[0][ell][i] = integral_trapezium_1D( VSH_integrand_r, finite_line_element_dt ) * sqrt_2Lplus1_pi[ell];
			
		}
		
		
		//--- (1)- and (2)-coefficients go from ell=1 ---
		for( int ell=1; ell<=ell_max; ell++ ){
		
			for( int j=0; j<n_points_t; j++ ){
				VSH_integrand_1[j] = v[1][i][j] * P1[j][ell];
				VSH_integrand_2[j] = v[2][i][j] * P1[j][ell];
			}
			
			v_VSH_coeffs[1][ell][i] = integral_trapezium_1D( VSH_integrand_1, finite_line_element_dt ) * sqrt_2Lplus1_pi_over_L_Lplus1[ell];
			v_VSH_coeffs[2][ell][i] = integral_trapezium_1D( VSH_integrand_2, finite_line_element_dt ) * sqrt_2Lplus1_pi_over_L_Lplus1[ell];
			
		}
		
	}
}




void VSH_decomposition_divergenceless( std::vector< std::vector< std::vector<double> > > &v, std::vector< std::vector< std::vector<double> > > &v_VSH_coeffs ){
	// Specific decomposition for divergence-free vectors, where the (1)-coeffs need not be calculated by integration.
	// Instead, they are given in terms of the r-coeffs, but this is only calculable once the radial derivatives are calculated.
	// THIS DOES NOT CALCULATE THE (1)-COEFFICIENTS!
	// First, build lists of integrand values. Then, numerically integrate and put the values into an array.
	// Two separate loops due to two separate ranges on ell.
	
	for( int i=0; i<n_points_r; i++ ){
		
		//--- r-coefficients go from ell=0 ---
		for( int ell=0; ell<=ell_max; ell++ ){
			
			for( int j=0; j<n_points_t; j++ ){
				VSH_integrand_r[j] = v[0][i][j] * P0[j][ell];
			}
			
			v_VSH_coeffs[0][ell][i] = integral_trapezium_1D( VSH_integrand_r, finite_line_element_dt ) * sqrt_2Lplus1_pi[ell];
			
		}
		
		
		//--- (2)-coefficients go from ell=1 ---
		for( int ell=1; ell<=ell_max; ell++ ){
		
			for( int j=0; j<n_points_t; j++ ){
				VSH_integrand_2[j] = v[2][i][j] * P1[j][ell];
			}
			
			v_VSH_coeffs[2][ell][i] = integral_trapezium_1D( VSH_integrand_2, finite_line_element_dt ) * sqrt_2Lplus1_pi_over_L_Lplus1[ell];
			
		}
		
	}
}




void calculate_radial_derivatives_of_VSH_coeffs_FD(){
	
	for( int L=0; L<=ell_max; L++ ){
		
		// Right derivatives at inner boundary.
		B_VSH_coeffs_dr [0][L][0] = ( B_VSH_coeffs[0][L][1] - B_VSH_coeffs[0][L][0] ) / delta_r;
		B_VSH_coeffs_dr [2][L][0] = ( B_VSH_coeffs[2][L][1] - B_VSH_coeffs[2][L][0] ) / delta_r;
		E_VSH_coeffs_dr [0][L][0] = ( E_VSH_coeffs[0][L][1] - E_VSH_coeffs[0][L][0] ) / delta_r;
		E_VSH_coeffs_dr [1][L][0] = ( E_VSH_coeffs[1][L][1] - E_VSH_coeffs[1][L][0] ) / delta_r;
		E_VSH_coeffs_dr [2][L][0] = ( E_VSH_coeffs[2][L][1] - E_VSH_coeffs[2][L][0] ) / delta_r;
		
		B_VSH_coeffs_dr2[0][L][0] = ( B_VSH_coeffs[0][L][2] - 2.0*B_VSH_coeffs[0][L][1] + B_VSH_coeffs[0][L][0] ) / ( delta_r*delta_r );
		
		// Symmetric derivatives for all intermediate points.
		for( int i=1; i<n_points_r-1; i++ ){
			
			B_VSH_coeffs_dr [0][L][i] = ( B_VSH_coeffs[0][L][i+1] - B_VSH_coeffs[0][L][i-1] ) / ( 2.0*delta_r );
			B_VSH_coeffs_dr [2][L][i] = ( B_VSH_coeffs[2][L][i+1] - B_VSH_coeffs[2][L][i-1] ) / ( 2.0*delta_r );
			E_VSH_coeffs_dr [0][L][i] = ( E_VSH_coeffs[0][L][i+1] - E_VSH_coeffs[0][L][i-1] ) / ( 2.0*delta_r );
			E_VSH_coeffs_dr [1][L][i] = ( E_VSH_coeffs[1][L][i+1] - E_VSH_coeffs[1][L][i-1] ) / ( 2.0*delta_r );
			E_VSH_coeffs_dr [2][L][i] = ( E_VSH_coeffs[2][L][i+1] - E_VSH_coeffs[2][L][i-1] ) / ( 2.0*delta_r );
					
			B_VSH_coeffs_dr2[0][L][i] = ( B_VSH_coeffs[0][L][i+1] - 2.0*B_VSH_coeffs[0][L][i] + B_VSH_coeffs[0][L][i-1] ) / ( delta_r*delta_r );
		}
		
		// Left derivatives at outer boundary.
		B_VSH_coeffs_dr[0][L].back() = ( B_VSH_coeffs[0][L].back() - B_VSH_coeffs[0][L][n_points_r-2] ) / delta_r;
		B_VSH_coeffs_dr[2][L].back() = ( B_VSH_coeffs[2][L].back() - B_VSH_coeffs[2][L][n_points_r-2] ) / delta_r;
		E_VSH_coeffs_dr[0][L].back() = ( E_VSH_coeffs[0][L].back() - E_VSH_coeffs[0][L][n_points_r-2] ) / delta_r;
		E_VSH_coeffs_dr[1][L].back() = ( E_VSH_coeffs[1][L].back() - E_VSH_coeffs[1][L][n_points_r-2] ) / delta_r;
		E_VSH_coeffs_dr[2][L].back() = ( E_VSH_coeffs[2][L].back() - E_VSH_coeffs[2][L][n_points_r-2] ) / delta_r;
		
		B_VSH_coeffs_dr2[0][L].back() = ( B_VSH_coeffs[0][L].back()  - 2.0*B_VSH_coeffs[0][L][n_points_r-2] + B_VSH_coeffs[0][L][n_points_r-3] ) / ( delta_r*delta_r );
		
	}

}




void calculate_radial_derivatives_of_VSH_coeffs_Chebyshev(){
	// DOES NOT YET CALCULATE SECOND DERIVATIVES 20240308
	
	for( int L=0; L<=ell_max; L++ ){
		
		for( int c=0; c<3; c++ ){
			
			std::vector<double> chebyshev_coeffs_B = chebyshev_series_coeffs( B_VSH_coeffs[c][L] );
			std::vector<double> chebyshev_coeffs_E = chebyshev_series_coeffs( E_VSH_coeffs[c][L] );
		
			for( int i=0; i<n_points_r; i++ ){
				B_VSH_coeffs_dr[c][L][i] = evaluate_chebyshev_series_first_derivative_at_gridpoint( chebyshev_coeffs_B, i );
				E_VSH_coeffs_dr[c][L][i] = evaluate_chebyshev_series_first_derivative_at_gridpoint( chebyshev_coeffs_E, i );
			}
		}
	}
}
		
		
		
		





void calculate_B_1_L_and_B_1_L_dr(){
	// This version of the function is for T>0, at which time we want to update B^{(1),L}.
	
	for( int L=1; L<=ell_max; L++ ){
		for( int i=0; i<n_points_r; i++ ){
			B_VSH_coeffs   [1][L][i] = ( r[i] * B_VSH_coeffs_dr [0][L][i] + 2.0 * B_VSH_coeffs   [0][L][i] ) / ( (double) L*(L+1.0) );
			B_VSH_coeffs_dr[1][L][i] = ( r[i] * B_VSH_coeffs_dr2[0][L][i] + 3.0 * B_VSH_coeffs_dr[0][L][i] ) / ( (double) L*(L+1.0) );
		}
	}
	
}




void recalculate_B_from_VSH_coeffs(){
	
	//--- Reset values ---
	B_from_VSH_coeffs = std::vector< std::vector< std::vector<double> > > ( 3, std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) ) );
	
	//--- Build values at each gridpoint ---
	for( int i=0; i<n_points_r; i++ ){
		for( int j=0; j<n_points_t; j++ ){
			
			//--- The r-coefficients start from ell=0 ---
			for( int ell=0; ell<=ell_max; ell++ ){
				B_from_VSH_coeffs[0][i][j] += B_VSH_coeffs[0][ell][i] * sqrt_2Lplus1_over_4pi[ell] * P0[j][ell];
			}
			
			//--- The (1)- and (2)-coefficients start from ell=1 ---
			for( int ell=1; ell<=ell_max; ell++ ){
				B_from_VSH_coeffs[1][i][j] += B_VSH_coeffs[1][ell][i] * sqrt_2Lplus1_over_4pi[ell] * P1[j][ell];
				B_from_VSH_coeffs[2][i][j] += B_VSH_coeffs[2][ell][i] * sqrt_2Lplus1_over_4pi[ell] * P1[j][ell];
			}
			
		}
	}
}




void recalculate_E_from_VSH_coeffs(){
	
	//--- Reset values ---
	E_from_VSH_coeffs = std::vector< std::vector< std::vector<double> > > ( 3, std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) ) );
	
	//--- Build values at each gridpoint ---
	for( int i=0; i<n_points_r; i++ ){
		for( int j=0; j<n_points_t; j++ ){
			
			//--- The r-coefficients start from ell=0 ---
			for( int ell=0; ell<=ell_max; ell++ ){
				E_from_VSH_coeffs[0][i][j] += E_VSH_coeffs[0][ell][i] * sqrt_2Lplus1_over_4pi[ell] * P0[j][ell];
			}
			
			//--- The (1)- and (2)-coefficients start from ell=1 ---
			for( int ell=1; ell<=ell_max; ell++ ){
				E_from_VSH_coeffs[1][i][j] += E_VSH_coeffs[1][ell][i] * sqrt_2Lplus1_over_4pi[ell] * P1[j][ell];
				E_from_VSH_coeffs[2][i][j] += E_VSH_coeffs[2][ell][i] * sqrt_2Lplus1_over_4pi[ell] * P1[j][ell];
			}
			
		}
	}
}




double stdev_of_array( std::vector<double> exact, std::vector<double> calc ){
	// Calculate standard deviation of an array of calculated results, compared to their known exact values.
	
	double ret = 0;
	
	for( int i=0; i<exact.size(); i++ ){
		ret += pow( exact[i] - calc[i], 2 );
	}
	
	return sqrt( ret / (double) exact.size() );
	
}




void calculate_stdev_of_VSH_decomposition(){
	
	//--- Reset values to zero ---
	B_abs_rel_dev     = std::vector< std::vector< std::vector<double> > > ( 3, std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) ) );
	E_abs_rel_dev     = std::vector< std::vector< std::vector<double> > > ( 3, std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) ) );
	B_stdev           = std::vector<double> ( 3 );
	E_stdev           = std::vector<double> ( 3 );
	B_max_abs_rel_dev = std::vector<double> ( 3 );
	E_max_abs_rel_dev = std::vector<double> ( 3 );
	
	//--- Recalculate vectors from VSH coefficients ---
	recalculate_B_from_VSH_coeffs();
	recalculate_E_from_VSH_coeffs();
	
	//--- Calculate stdev ---
	for( int f=0; f<3; f++ ){
		
		for( int i=0; i<n_points_r; i++ ){
			for( int j=0; j<n_points_t; j++ ){
				
				B_stdev[f] += pow( B[f][i][j] - B_from_VSH_coeffs[f][i][j], 2 );
				E_stdev[f] += pow( E[f][i][j] - E_from_VSH_coeffs[f][i][j], 2 );
				
				if( B[f][i][j] == 0 ){
					B_abs_rel_dev[f][i][j] = 0;
				} else {
					B_abs_rel_dev[f][i][j] = abs( 1.0 - B_from_VSH_coeffs[f][i][j] / B[f][i][j] );
				}
				
				if( E[f][i][j] == 0 ){
					E_abs_rel_dev[f][i][j] = 0;
				} else {
					E_abs_rel_dev[f][i][j] = abs( 1.0 - E_from_VSH_coeffs[f][i][j] / E[f][i][j] );
				}
				
				if( B_abs_rel_dev[f][i][j] > B_max_abs_rel_dev[f] ){
					B_max_abs_rel_dev[f] = B_abs_rel_dev[f][i][j];
				}
				
				if( E_abs_rel_dev[f][i][j] > E_max_abs_rel_dev[f] ){
					E_max_abs_rel_dev[f] = E_abs_rel_dev[f][i][j];
				}
				
			}
		}
		
		B_stdev[f] = sqrt( B_stdev[f] / ( (double) n_points_total ) );
		E_stdev[f] = sqrt( E_stdev[f] / ( (double) n_points_total ) );
		
	}
	
}




void calculate_radial_derivatives_of_vector_components(){
	// Not used by the code, but useful for CSV output to check for strange behaviour in dB/dr and dE/dr.
	// To do <20230824: Second derivatives too.
	
	//--- Reset values to zero ---
	B_dr  = std::vector< std::vector< std::vector<double> > > ( 3, std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) ) );
	E_dr  = std::vector< std::vector< std::vector<double> > > ( 3, std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) ) );
	B_dr2 = std::vector< std::vector< std::vector<double> > > ( 3, std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) ) );
	
	//--- Sum over ell for each gridpoint ---
	for( int i=0; i<n_points_r; i++ ){
		for( int j=0; j<n_points_t; j++ ){
			
			for( int ell=0; ell<=ell_max; ell++ ){
				B_dr [0][i][j] += sqrt_2Lplus1_over_4pi[ell] * P0[j][ell] * B_VSH_coeffs_dr [0][ell][i];
				E_dr [0][i][j] += sqrt_2Lplus1_over_4pi[ell] * P0[j][ell] * E_VSH_coeffs_dr [0][ell][i];
				B_dr2[0][i][j] += sqrt_2Lplus1_over_4pi[ell] * P0[j][ell] * B_VSH_coeffs_dr2[0][ell][i];
			}
			
			for( int ell=1; ell<=ell_max; ell++ ){
				B_dr [1][i][j] += sqrt_2Lplus1_over_4pi[ell] * P1[j][ell] * B_VSH_coeffs_dr [1][ell][i];
				B_dr [2][i][j] += sqrt_2Lplus1_over_4pi[ell] * P1[j][ell] * B_VSH_coeffs_dr [2][ell][i];
				E_dr [1][i][j] += sqrt_2Lplus1_over_4pi[ell] * P1[j][ell] * E_VSH_coeffs_dr [1][ell][i];
				E_dr [2][i][j] += sqrt_2Lplus1_over_4pi[ell] * P1[j][ell] * E_VSH_coeffs_dr [2][ell][i];
			}
			
		}
	}
}




void calculate_Cartesian_vector_components(){
	
	for( int i=0; i<n_points_r; i++ ){
		for( int j=0; j<n_points_t; j++ ){
			
			B_Cartesian[0][i][j] = B[0][i][j] * st[j] + B[1][i][j] * ct[j];
			B_Cartesian[1][i][j] = B[2][i][j];									// Extra variable assignment not necessary, but use for now to avoid confusion.
			B_Cartesian[2][i][j] = B[0][i][j] * ct[j] - B[1][i][j] * st[j];
			
			E_Cartesian[0][i][j] = E[0][i][j] * st[j] + E[1][i][j] * ct[j];
			E_Cartesian[1][i][j] = E[2][i][j];									// Extra variable assignment not necessary, but use for now to avoid confusion.
			E_Cartesian[2][i][j] = E[0][i][j] * ct[j] - E[1][i][j] * st[j];
			
			J_Cartesian[0][i][j] = J[0][i][j] * st[j] + J[1][i][j] * ct[j];
			J_Cartesian[1][i][j] = J[2][i][j];									// Extra variable assignment not necessary, but use for now to avoid confusion.
			J_Cartesian[2][i][j] = J[0][i][j] * ct[j] - J[1][i][j] * st[j];
			
		}
	}
}




void calculate_dot_products(){
	
	for( int i=0; i<n_points_r; i++ ){
		for( int j=0; j<n_points_t; j++ ){
			B_squared                [i][j] = pow( B[0][i][j], 2 ) + pow( B[1][i][j], 2 ) + pow( B[2][i][j], 2 );
			E_squared                [i][j] = pow( E[0][i][j], 2 ) + pow( E[1][i][j], 2 ) + pow( E[2][i][j], 2 );
			B_squared_minus_E_squared[i][j] = B_squared[i][j] - E_squared[i][j];
			E_dot_B                  [i][j] = E[0][i][j]*B[0][i][j] + E[1][i][j]*B[1][i][j] + E[2][i][j]*B[2][i][j];
		}
	}
	
}




void ramp_up_electric_fields_for_rotation(){
	if( ( T_index >= T_index_rotation_ramp_start ) and ( T_index <= T_index_rotation_ramp_stop ) ){
		
		//--- V1 (Header_Time_Evolution_24.h and OneNote notes 20231020Fri) ---
		/*
		for( int i=0; i<n_points_r; i++ ){
			for( int j=1; j<n_points_t-1; j++ ){
				E[0][i][j] += dOmega_by_dT * delta_T * r[i] * st[j] *  B[1][i][j];
				E[1][i][j] += dOmega_by_dT * delta_T * r[i] * st[j] * -B[0][i][j];
			}
		}
		*/
		
		//--- V2 (Header_Time_Evolution_25.h and OneNote Notes 20231115 and PhD/20231115 Expression for electric field ramp-up/Expression for electric field ramp-up.pdf and log file 20231208_c_run_to_T_equals_4_and_rewrite_rampup.txt) ---
		
		for( int i=0; i<n_points_r; i++ ){
			for( int j=1; j<n_points_t-1; j++ ){
				//E[0][i][j] += ( Omega - Omega_prev ) * pow( r[i], -2 ) * pow( st[j], 2 );
				//E[1][i][j] += ( Omega - Omega_prev ) * pow( r[i], -2 ) * -2.0 * st[j] * cos(t[j]);
			}
		}
		
		
		apply_inner_boundary_conditions();
		apply_outer_boundary_conditions();
	}
}




void calculate_n_gridpoints_within_FF_region(){
	
	is_gridpoint_within_FF_region = std::vector< std::vector<int> > ( n_points_r, std::vector<int> ( n_points_t ) );
	n_gridpoints_within_FF_region = 0;
	
	if( Omega[T_index] > 0 ){
		for( int i=0; i<n_points_r; i++ ){
			for( int j=1; j<n_points_t-1; j++ ){
				if( r[i] * pow( st[j], -2 ) < R_LC[T_index] ){
					is_gridpoint_within_FF_region[i][j]  = 1;
					n_gridpoints_within_FF_region       += 1;
				}
			}
		}
	}
}



void apply_force_free_conditions_within_force_free_region(){
	// Update the electric field at all points within the last closed field line of a rotating dipole, except the poles,
	// everywhere the first or second force-free condition is violated, with an arbitrary tolerance for the value of E dot B.
	// Petri 2012, Section 3.2.
	
	calculate_n_gridpoints_within_FF_region();
	
	num_E_points_changed_for_second_FF_condition = 0;
	were_FF_conditions_applied = std::vector< std::vector<int> > ( n_points_r, std::vector<int> ( n_points_t ) );
	
	if( Omega[T_index] > 0 ){
		for( int i=0; i<n_points_r; i++ ){
			for( int j=1; j<n_points_t-1; j++ ){
				
				if( is_gridpoint_within_FF_region[i][j] ){
					
					double E_dot_B_tol = 1e-3;	// Arbitrarily chosen 20230921.
					
					if( ( E_dot_B[i][j] > E_dot_B_tol ) or ( B_squared_minus_E_squared[i][j] < 0 ) ){
						
						num_E_points_changed_for_second_FF_condition += 1;
						
						double E_dot_B_over_Bsq = E_dot_B[i][j] / B_squared[i][j];
						double E_prime_r = E[0][i][j] - E_dot_B_over_Bsq * B[0][i][j];
						double E_prime_t = E[1][i][j] - E_dot_B_over_Bsq * B[1][i][j];
						double E_prime_p = E[2][i][j] - E_dot_B_over_Bsq * B[2][i][j];
						double sqrt_Bsq_over_Eprimesq = sqrt( B_squared[i][j] / ( E_prime_r*E_prime_r + E_prime_t*E_prime_t + E_prime_p*E_prime_p ) );
						
						E[0][i][j] = E_prime_r * sqrt_Bsq_over_Eprimesq;
						E[1][i][j] = E_prime_t * sqrt_Bsq_over_Eprimesq;
						E[2][i][j] = E_prime_p * sqrt_Bsq_over_Eprimesq;
						
						were_FF_conditions_applied[i][j] = 1;
						
					}
				}
					
				else {
					were_FF_conditions_applied[i][j] = -1;
				}

			}
		}
	}
	
	// Calculate percentage of gridpoints at which E changed. Below format makes it exactly one decimal place.
	percentage_E_points_points_changed_for_second_FF_condition = round( (double)1000 * num_E_points_changed_for_second_FF_condition / n_gridpoints_within_FF_region ) / 10;
}




void calculate_time_derivatives_of_vector_components(){
	
	B_curl = std::vector< std::vector< std::vector<double> > > ( 3, std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) ) );
	E_curl = std::vector< std::vector< std::vector<double> > > ( 3, std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) ) );
	J      = std::vector< std::vector< std::vector<double> > > ( 3, std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) ) );
	E_div = std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) );
									 
	for( int i=0; i<n_points_r; i++ ){
		for( int j=1; j<n_points_t-1; j++ ){
			
			//--- Spatial derivatives ---
			//- L=0 -
			E_div[i][j] += sqrt_2Lplus1_over_4pi[0] * ( E_VSH_coeffs_dr[0][0][i] + 2.0*pow(r[i],-1) * E_VSH_coeffs[0][0][i] );
			
			//- L>0 -
			for( int L=1; L<=ell_max; L++ ){
				
				E_div[i][j] += sqrt_2Lplus1_over_4pi[L] * P0[j][L] * ( E_VSH_coeffs_dr[0][L][i] + 2.0*pow(r[i],-1) * E_VSH_coeffs[0][L][i] - L*(L+1.0)*pow(r[i],-1) * E_VSH_coeffs[1][L][i] );
				
				B_curl[0][i][j] += sqrt_2Lplus1_over_4pi[L] * P0[j][L] * -L*(L+1.0) * B_VSH_coeffs[2][L][i] / r[i];
				B_curl[1][i][j] += sqrt_2Lplus1_over_4pi[L] * P1[j][L] * ( -B_VSH_coeffs_dr[2][L][i] - B_VSH_coeffs[2][L][i] / r[i] );
				B_curl[2][i][j] += sqrt_2Lplus1_over_4pi[L] * P1[j][L] * ( r[i] * B_VSH_coeffs_dr2[0][L][i] + 4.0 * B_VSH_coeffs_dr[0][L][i] - (L+2.0)*(L-1.0) * B_VSH_coeffs[0][L][i] / r[i] ) / ( (double) L*(L+1.0) );
				
				E_curl[0][i][j] += sqrt_2Lplus1_over_4pi[L] * P0[j][L] * -L*(L+1.0) * E_VSH_coeffs[2][L][i] / r[i];
				E_curl[1][i][j] += sqrt_2Lplus1_over_4pi[L] * P1[j][L] * ( -E_VSH_coeffs_dr[2][L][i] - E_VSH_coeffs[2][L][i] / r[i] );
				E_curl[2][i][j] += sqrt_2Lplus1_over_4pi[L] * P1[j][L] * ( E_VSH_coeffs_dr[1][L][i] + ( - E_VSH_coeffs[0][L][i] + E_VSH_coeffs[1][L][i] ) / r[i] );
				
			}
			
			//--- E cross B ---
			std::vector<double> E_cross_B (3);
			E_cross_B[0] = E[1][i][j] * B[2][i][j] - E[2][i][j] * B[1][i][j];
			E_cross_B[1] = E[2][i][j] * B[0][i][j] - E[0][i][j] * B[2][i][j];
			E_cross_B[2] = E[0][i][j] * B[1][i][j] - E[1][i][j] * B[0][i][j];
			
			//--- J ---
			double B_dot_curl_B = B[0][i][j]*B_curl[0][i][j] + B[1][i][j]*B_curl[1][i][j] + B[2][i][j]*B_curl[2][i][j];
			double E_dot_curl_E = E[0][i][j]*E_curl[0][i][j] + E[1][i][j]*E_curl[1][i][j] + E[2][i][j]*E_curl[2][i][j];
			
			for( int f=0; f<3; f++ ){
				J[f][i][j] = ( E_div[i][j] * E_cross_B[f] + ( B_dot_curl_B - E_dot_curl_E ) * B[f][i][j] ) / B_squared[i][j];
			}
			
			//--- Time derivatives ---
			for( int f=0; f<3; f++ ){
				B_dT[f][i][j] = - E_curl[f][i][j];
				E_dT[f][i][j] =   B_curl[f][i][j] - J[f][i][j];
			}
			
			//--- Frictional term, applied to all components except B_r_dT ---
			if( use_outer_sponge_layer and ( r[i] >= r_sponge ) ){
				B_dT[1][i][j] -= friction[i][j] * B[1][i][j];
				B_dT[2][i][j] -= friction[i][j] * B[2][i][j];
				E_dT[0][i][j] -= friction[i][j] * E[0][i][j];
				E_dT[1][i][j] -= friction[i][j] * E[1][i][j];
				E_dT[2][i][j] -= friction[i][j] * E[2][i][j];
				
			}	
		}		
	}

	//--- 20231214 Specific expressions for inner boundary ---

	B_t_dt_inner[0]     = ( B[1][0][1]     - B[1][0][0]            ) / delta_t;
	B_t_dt_inner.back() = ( B[1][0].back() - B[1][0][n_points_t-2] ) / delta_t;
        for( int j=1; j<n_points_t-1; j++ ){
	  B_t_dt_inner[j] = ( B[1][0][j+1] - B[1][0][j-1] ) / ( 2.0*delta_t );
	}

	for( int j=0; j<n_points_t; j++ ){
	  B_dT[0][0][j] = 0;
	  B_dT[1][0][j] = 0;
	  B_dT[2][0][j] = Omega[T_index] * ( 2.0*st[j]*B[0][0][j] + ct[j]*B[1][0][j] + st[j]*B_dr[0][0][j] + st[j]*B_t_dt_inner[j] );
	}
	
	//--- Set the time-derivatives of the components with boundary conditions to zero, to avoid unwanted evolution ---
	for( int j=0; j<n_points_t; j++ ){
	  //B_dT[0][0]    [j] = 0;
	  //B_dT[2][0]    [j] = 0;
		E_dT[1][0]    [j] = 0;
		E_dT[2][0]    [j] = 0;
		//B_dT[1].back()[j] = 0;
		//B_dT[2].back()[j] = 0;
		E_dT[1].back()[j] = 0;
		E_dT[2].back()[j] = 0;
	}
	//B_dT[2][0][ (int) n_points_t / 2 - 1 ] = 0;	// Set B_phi ( R_NS, pi/2 ) = 0 to avoid evolution of noise. Suggestion from Sam 20231208.
}




void integrate_vectors_wrt_time_Adams_Bashforth_Order_1(){
	// For the Adams-Bashforth method. When we increase the timestep, all of the n=0 values should become n=1 values, and all of the n=1 values should become n=2 values.
	// This is the first-order AB method, equivalent to the Euler method. It is used at the first timestep only, so we don't need to consider two steps ago.
	
	for( int i=0; i<n_points_r; i++ ){
		for( int j=1; j<n_points_t-1; j++ ){
			
			//--- Setup the values at the previous timestep ---
			for( int f=0; f<3; f++ ){
				B_minus1   [f][i][j] = B   [f][i][j];
				E_minus1   [f][i][j] = E   [f][i][j];
				B_dT_minus1[f][i][j] = B_dT[f][i][j];
				E_dT_minus1[f][i][j] = E_dT[f][i][j];
			}
			
			//--- Calculate the values at the current timestep by Adams-Bashforth integration ---
			for( int f=0; f<3; f++ ){
				B[f][i][j] = B_minus1[f][i][j] + delta_T * B_dT_minus1[f][i][j];
				E[f][i][j] = E_minus1[f][i][j] + delta_T * E_dT_minus1[f][i][j];
			}
			
		}
	}
	
}




void integrate_vectors_wrt_time_Adams_Bashforth_Order_2(){
	// For the Adams-Bashforth method. When we increase the timestep, all of the n=0 values should become n=1 values, and all of the n=1 values should become n=2 values.
	
	for( int i=0; i<n_points_r; i++ ){
		for( int j=1; j<n_points_t-1; j++ ){
			
			//--- Setup the values two steps ago ---
			for( int f=0; f<3; f++ ){
				B_minus2   [f][i][j] = B_minus1   [f][i][j];
				E_minus2   [f][i][j] = E_minus1   [f][i][j];
				B_dT_minus2[f][i][j] = B_dT_minus1[f][i][j];
				E_dT_minus2[f][i][j] = E_dT_minus1[f][i][j];
			}
			
			//--- Setup the values at the previous timestep ---
			for( int f=0; f<3; f++ ){
				B_minus1   [f][i][j] = B   [f][i][j];
				E_minus1   [f][i][j] = E   [f][i][j];
				B_dT_minus1[f][i][j] = B_dT[f][i][j];
				E_dT_minus1[f][i][j] = E_dT[f][i][j];
			}
			
			//--- Calculate the values at the current timestep by Adams-Bashforth integration ---
			for( int f=0; f<3; f++ ){
				B[f][i][j] = B_minus1[f][i][j] + delta_T * ( 1.5 * B_dT_minus1[f][i][j] - 0.5 * B_dT_minus2[f][i][j] );
				E[f][i][j] = E_minus1[f][i][j] + delta_T * ( 1.5 * E_dT_minus1[f][i][j] - 0.5 * E_dT_minus2[f][i][j] );
			}
			
		}
	}
	
}




void integrate_vectors_wrt_time_Adams_Bashforth_Order_3(){
	// For the Adams-Bashforth method. When we increase the timestep, all of the n=0 values should become n=1 values, and all of the n=1 values should become n=2 values.
	
	for( int i=0; i<n_points_r; i++ ){
		for( int j=1; j<n_points_t-1; j++ ){
			
			//--- Setup the values three steps ago ---
			for( int f=0; f<3; f++ ){
				B_minus3   [f][i][j] = B_minus2   [f][i][j];
				E_minus3   [f][i][j] = E_minus2   [f][i][j];
				B_dT_minus3[f][i][j] = B_dT_minus2[f][i][j];
				E_dT_minus3[f][i][j] = E_dT_minus2[f][i][j];
			}
			
			//--- Setup the values two steps ago ---
			for( int f=0; f<3; f++ ){
				B_minus2   [f][i][j] = B_minus1   [f][i][j];
				E_minus2   [f][i][j] = E_minus1   [f][i][j];
				B_dT_minus2[f][i][j] = B_dT_minus1[f][i][j];
				E_dT_minus2[f][i][j] = E_dT_minus1[f][i][j];
			}
			
			//--- Setup the values at the previous timestep ---
			for( int f=0; f<3; f++ ){
				B_minus1   [f][i][j] = B   [f][i][j];
				E_minus1   [f][i][j] = E   [f][i][j];
				B_dT_minus1[f][i][j] = B_dT[f][i][j];
				E_dT_minus1[f][i][j] = E_dT[f][i][j];
			}
	
			//--- Calculate the values at the current timestep by Adams-Bashforth integration ---
			for( int f=0; f<3; f++ ){
				B[f][i][j] = B_minus1[f][i][j] + delta_T * ( 23.0 * B_dT_minus1[f][i][j] - 16.0 * B_dT_minus2[f][i][j] + 5.0 * B_dT_minus3[f][i][j] ) / 12.0;
				E[f][i][j] = E_minus1[f][i][j] + delta_T * ( 23.0 * E_dT_minus1[f][i][j] - 16.0 * E_dT_minus2[f][i][j] + 5.0 * E_dT_minus3[f][i][j] ) / 12.0;
			}
			
		}
	}
}




void calculate_magnetic_energy(){
	
	for( int i=0; i<n_points_r; i++ ){
		for( int j=0; j<n_points_t; j++ ){
			magnetic_energy_density[i][j] = 0.5 * ( pow( B[0][i][j], 2 ) + pow( B[1][i][j], 2 ) + pow( B[2][i][j], 2 ) );
		}
	}
	
	total_magnetic_energy = 2.0 * pi * integral_trapezium_2D( magnetic_energy_density, finite_area_element_for_r_t );
	
}




double estimate_csv_profiles_filesize(){
	
	double fudge_factor = 17.0;		// Free parameter to tune final result.
	int n_columns = 19 + 8 * 6;		// You are responsible for keeping this value updated based on how many values are specified in void output_headers_to_csv_profiles().
	int n_rows = 0;					// Calculated by the for-loop below.
	
	for( int T_index_here=0; T_index_here<=n_T; T_index_here++ ){
		bool condition_T_1 = ( T_index_here >= csv_profiles_write_T_min ) and ( T_index_here <= csv_profiles_write_T_max );
		bool condition_T_2 = ( T_index_here % csv_profiles_write_freq_T == 0 ) or ( T_index_here == n_T );
		if( condition_T_1 and condition_T_2 ){
			for( int i=csv_profiles_write_i_min; i<n_points_r; i++ ){
				for( int j=0; j<n_points_t; j++ ){
					bool condition_r = ( i <= csv_profiles_write_i_max ) and ( ( i % csv_profiles_write_freq_r == 0 ) or ( i == n_points_r - 1 ) );
					bool condition_t = ( j % csv_profiles_write_freq_t == 0 ) or ( j == n_points_t - 1 );
					if( condition_r and condition_t ){
						n_rows += 1;
					}
				}
			}
		}
	}
	
	return n_rows * n_columns * pow( 2, -20 ) * fudge_factor;	// 2^-20 goes from bytes to megabytes.
}
		




void output_headers_to_csv_profiles(){
	// Split over separate lines for readability and ease of checking consistency with output_to_csv_profiles().
	
	// Remember to update int n_columns in void estimate_csv_filesize().
	
	output_file_profiles << "T_index,T,T_SI,i,j,r,t,x,z";
	output_file_profiles << ",B,E,J,E_dot_B,B_squared_minus_E_squared,B_r_dr2,magnetic_energy_density,friction,is_gridpoint_within_FF_region,were_FF_conditions_applied";
	
	output_file_profiles << ",B_r,E_r,B_t,E_t,B_p,E_p";
	output_file_profiles << ",B_x,E_x,B_y,E_y,B_z,E_z";
	output_file_profiles << ",J_r,J_x,J_t,J_y,J_p,J_z";
	output_file_profiles << ",B_r_dr,E_r_dr,B_t_dr,E_t_dr,B_p_dr,E_p_dr";
	output_file_profiles << ",abs_rel_dev_B_r,abs_rel_dev_E_r,abs_rel_dev_B_t,abs_rel_dev_E_t,abs_rel_dev_B_p,abs_rel_dev_E_p";
	output_file_profiles << ",curl_B_r,curl_E_r,curl_B_t,curl_E_t,curl_B_p,curl_E_p";
	output_file_profiles << ",B_r_dT,E_r_dT,B_t_dT,E_t_dT,B_p_dT,E_p_dT";
	output_file_profiles << ",B_r_dT_over_B_r,E_r_dT_over_E_r,B_t_dT_over_B_t,E_t_dT_over_E_t,B_p_dT_over_B_p,E_p_dT_over_E_p";
	output_file_profiles << ",B_t_dt_inner";
	
	output_file_profiles << "\n";
}




void output_to_csv_profiles( bool condition_override = false ){
	
	bool condition_T_1 = ( T_index >= csv_profiles_write_T_min ) and ( T_index <= csv_profiles_write_T_max );
	bool condition_T_2 = ( T_index % csv_profiles_write_freq_T == 0 ) or ( T_index == n_T );
	
	if( condition_override ){
		condition_T_1 = true;
		condition_T_2 = true;
	}
	
	if( condition_T_1 and condition_T_2 ){
	
		for( int i=csv_profiles_write_i_min; i<n_points_r; i++ ){
			for( int j=0; j<n_points_t; j++ ){
				
				bool condition_r = ( i <= csv_profiles_write_i_max ) and ( ( i % csv_profiles_write_freq_r == 0 ) or ( i == n_points_r - 1 ) );
				bool condition_t = ( j % csv_profiles_write_freq_t == 0 ) or ( j == n_points_t - 1 );
				
				if( condition_override ){
					condition_r = true;
					condition_t = true;
				}
		
				if( condition_r and condition_t ){
					
					double tolerance_for_zero = 1e-13;
					
					std::vector<double> B_dT_over_B ( 3 );
					std::vector<double> E_dT_over_E ( 3 );
					
					for( int f=0; f<3; f++ ){
						if( abs( B[f][i][j] ) > tolerance_for_zero ){
							B_dT_over_B[f] = B_dT[f][i][j] / B[f][i][j];
						}
						if( abs( E[f][i][j] ) > tolerance_for_zero ){
							E_dT_over_E[f] = E_dT[f][i][j] / E[f][i][j];
						}
					}
					
					
					output_file_profiles <<std::setprecision(csv_precision)<< T_index <<","<< T[T_index] <<","<< T_SI[T_index] <<","<< i <<","<< j <<","<< r[i] <<","<< t[j] <<","<< x[i][j] <<","<< z[i][j]
																		   <<","<< sqrt(B_squared[i][j]) <<","<< sqrt(E_squared[i][j]) <<","<< sqrt( pow(J[0][i][j],2) + pow(J[1][i][j],2) + pow(J[2][i][j],2) )
																		   <<","<< E_dot_B[i][j] <<","<< B_squared_minus_E_squared[i][j]
																		   <<","<< B_dr2[0][i][j]
																		   <<","<< magnetic_energy_density[i][j]
																		   <<","<< friction[i][j]
																		   <<","<< is_gridpoint_within_FF_region[i][j] <<","<< were_FF_conditions_applied[i][j];
					
					for( int f=0; f<3; f++ ){ output_file_profiles <<std::setprecision(csv_precision)<<","<< B            [f][i][j] <<std::setprecision(csv_precision)<<","<< E            [f][i][j]; }
					for( int f=0; f<3; f++ ){ output_file_profiles <<std::setprecision(csv_precision)<<","<< B_Cartesian  [f][i][j] <<std::setprecision(csv_precision)<<","<< E_Cartesian  [f][i][j]; }
					for( int f=0; f<3; f++ ){ output_file_profiles <<std::setprecision(csv_precision)<<","<< J            [f][i][j] <<std::setprecision(csv_precision)<<","<< J_Cartesian  [f][i][j]; }
					for( int f=0; f<3; f++ ){ output_file_profiles <<std::setprecision(csv_precision)<<","<< B_dr         [f][i][j] <<std::setprecision(csv_precision)<<","<< E_dr         [f][i][j]; }
					for( int f=0; f<3; f++ ){ output_file_profiles <<std::setprecision(csv_precision)<<","<< B_abs_rel_dev[f][i][j] <<std::setprecision(csv_precision)<<","<< E_abs_rel_dev[f][i][j]; }
					for( int f=0; f<3; f++ ){ output_file_profiles <<std::setprecision(csv_precision)<<","<< B_curl       [f][i][j] <<std::setprecision(csv_precision)<<","<< E_curl       [f][i][j]; }
					for( int f=0; f<3; f++ ){ output_file_profiles <<std::setprecision(csv_precision)<<","<< B_dT         [f][i][j] <<std::setprecision(csv_precision)<<","<< E_dT         [f][i][j]; }
					for( int f=0; f<3; f++ ){ output_file_profiles <<std::setprecision(csv_precision)<<","<< B_dT_over_B  [f]       <<std::setprecision(csv_precision)<<","<< E_dT_over_E  [f]      ; }

					output_file_profiles << B_t_dt_inner[j];
					
					output_file_profiles << "\n";
					
				}
			}
		}
	}
}




void count_nans(){
	
	B_nans = std::vector<int> ( 3 );
	E_nans = std::vector<int> ( 3 );
	
	for( int i=0; i<n_points_r; i++ ){
		for( int j=0; j<n_points_t; j++ ){
			for( int f=0; f<3; f++ ){
				if( ( std::isnan( B[f][i][j] ) ) or ( std::isinf( B[f][i][j] ) ) ){
					B_nans[f] += 1;
				}
				if( ( std::isnan( E[f][i][j] ) ) or ( std::isinf( E[f][i][j] ) ) ){
					E_nans[f] += 1;
				}
			}
		}
	}
	
	nans_tot = B_nans[0] + B_nans[1] + B_nans[2] + E_nans[0] + E_nans[1] + E_nans[2];
}




void output_headers_to_csv_history(){
	// Split over separate lines for readability and ease of checking consistency with output_to_csv_history().
	
	output_file_history << "T_index,T,T_SI";
	output_file_history << ",Omega,R_LC,star_rotation_angle";
	output_file_history << ",E_pts_chngd,pct_E_pts_chngd";
	output_file_history << ",stdev_B_r,stdev_E_r,stdev_B_t,stdev_E_t,stdev_B_p,stdev_E_p";
	output_file_history << ",max_abs_rel_dev_B_r,max_abs_rel_dev_E_r,max_abs_rel_dev_B_t,max_abs_rel_dev_E_t,max_abs_rel_dev_B_p,max_abs_rel_dev_E_p";
	output_file_history << ",nans_B_r,nans_E_r,nans_B_t,nans_E_t,nans_B_p,nans_E_p";
	output_file_history << ",total_magnetic_energy";
	output_file_history << "\n";
}




void output_to_csv_history(){
	
	bool condition_T = ( T_index % csv_history_write_freq_T == 0 ) or ( T_index == n_T );
	
	if( condition_T ){
		
		double R_LC = 0;
		if( Omega[T_index] > 0 ){
			R_LC = 1.0 / Omega[T_index];
		}
		
		output_file_history <<std::setprecision(csv_precision)<< T_index <<","<< T[T_index] <<","<< T_SI[T_index]
		                                                      <<","<< Omega[T_index] <<","<< R_LC <<","<< star_rotation_angle[T_index]
															  <<","<< num_E_points_changed_for_second_FF_condition <<","<< percentage_E_points_points_changed_for_second_FF_condition;
		
		for( int f=0; f<3; f++ ){ output_file_history <<std::setprecision(csv_precision)<<","<< B_stdev          [f] <<std::setprecision(csv_precision)<<","<< E_stdev          [f]; }
		for( int f=0; f<3; f++ ){ output_file_history <<std::setprecision(csv_precision)<<","<< B_max_abs_rel_dev[f] <<std::setprecision(csv_precision)<<","<< E_max_abs_rel_dev[f]; }
		for( int f=0; f<3; f++ ){ output_file_history <<std::setprecision(csv_precision)<<","<< B_nans           [f] <<std::setprecision(csv_precision)<<","<< E_nans           [f]; }
		
		output_file_history <<std::setprecision(csv_precision)<<","<< nans_tot <<","<< total_magnetic_energy << "\n";
							
	}
}




void output_headers_to_csv_VSH_coeffs(){
	// Split over separate lines for readability and ease of checking consistency with output_to_csv_history().
	
	output_file_VSH_coeffs << "T_index,i,r,ell,B_r_L,E_r_L,B_1_L,E_1_L,B_2_L,E_2_L,B_r_L_dr,E_r_L_dr,B_1_L_dr,E_1_L_dr,B_2_L_dr,E_2_L_dr,B_r_L_dr2\n";
	
}




void output_to_csv_VSH_coeffs(){
	
	bool condition_T = ( T_index % csv_VSH_coeffs_write_freq_T == 0 ) or ( T_index == n_T );
	
	if( condition_T ){
		for( int i=0; i<n_points_r; i++ ){
			for( int ell=0; ell<=ell_max; ell++ ){
				
				output_file_VSH_coeffs <<std::setprecision(csv_precision)<< T_index <<","<< i <<","<< r[i] <<","<< ell;
				
				for( int f=0; f<3; f++ ){ output_file_VSH_coeffs <<std::setprecision(csv_precision)<<","<< B_VSH_coeffs   [f][ell][i] <<std::setprecision(csv_precision)<<","<< E_VSH_coeffs   [f][ell][i]; }
				for( int f=0; f<3; f++ ){ output_file_VSH_coeffs <<std::setprecision(csv_precision)<<","<< B_VSH_coeffs_dr[f][ell][i] <<std::setprecision(csv_precision)<<","<< E_VSH_coeffs_dr[f][ell][i]; }
				
				output_file_VSH_coeffs <<std::setprecision(csv_precision)<< ","<< B_VSH_coeffs_dr2[0][ell][i] << "\n";
				
			}
		}
	}
	
}




void output_headers_to_csv_BCs(){
	// Split over separate lines for readability and ease of checking consistency with output_to_csv_profiles().
	
	// Remember to update int n_columns in void estimate_csv_filesize().
	
	output_file_BCs << "T_index,T,T_SI,j,t";
	output_file_BCs << ",B_r_inner_before_BCs,E_r_inner_before_BCs,B_t_inner_before_BCs,E_t_inner_before_BCs,B_p_inner_before_BCs,E_p_inner_before_BCs";
	output_file_BCs << ",B_r_inner_after_BCs,E_r_inner_after_BCs,B_t_inner_after_BCs,E_t_inner_after_BCs,B_p_inner_after_BCs,E_p_inner_after_BCs";
	output_file_BCs << ",rel_diff_B_r_inner,rel_diff_E_r_inner,rel_diff_B_t_inner,rel_diff_E_t_inner,rel_diff_B_p_inner,rel_diff_E_p_inner";
	output_file_BCs << ",B_r_outer_before_BCs,E_r_outer_before_BCs,B_t_outer_before_BCs,E_t_outer_before_BCs,B_p_outer_before_BCs,E_p_outer_before_BCs";
	output_file_BCs << ",B_r_outer_after_BCs,E_r_outer_after_BCs,B_t_outer_after_BCs,E_t_outer_after_BCs,B_p_outer_after_BCs,E_p_outer_after_BCs";
	output_file_BCs << ",rel_diff_B_r_outer,rel_diff_E_r_outer,rel_diff_B_t_outer,rel_diff_E_t_outer,rel_diff_B_p_outer,rel_diff_E_p_outer";
	output_file_BCs << "\n";
	
}




void output_to_csv_BCs(){
	
	bool condition_T = ( T_index % csv_BCs_write_freq_T == 0 ) or ( T_index == n_T );
	
	if( condition_T ){
		
		for( int j=0; j<n_points_t; j++ ){
			
			std::vector<double> B_rel_diff_inner( 3 );
			std::vector<double> E_rel_diff_inner( 3 );
			std::vector<double> B_rel_diff_outer( 3 );
			std::vector<double> E_rel_diff_outer( 3 );
			
			for( int f=0; f<3; f++ ){
				B_rel_diff_inner[f] = 1.0 - B_inner_after_BCs[f][j] / B_inner_before_BCs[f][j];
				E_rel_diff_inner[f] = 1.0 - E_inner_after_BCs[f][j] / E_inner_before_BCs[f][j];
				B_rel_diff_outer[f] = 1.0 - B_outer_after_BCs[f][j] / B_outer_before_BCs[f][j];
				E_rel_diff_outer[f] = 1.0 - E_outer_after_BCs[f][j] / E_outer_before_BCs[f][j];
			}
		
			output_file_BCs <<std::setprecision(csv_precision)<< T_index <<","<< T[T_index] <<","<< T_SI[T_index] <<","<< j <<","<< t[j];
			
			for( int f=0; f<3; f++ ){ output_file_BCs <<std::setprecision(csv_precision)<<","<< B_inner_before_BCs[f][j] <<std::setprecision(csv_precision)<<","<< E_inner_before_BCs[f][j]; }
			for( int f=0; f<3; f++ ){ output_file_BCs <<std::setprecision(csv_precision)<<","<< B_inner_after_BCs [f][j] <<std::setprecision(csv_precision)<<","<< E_inner_after_BCs [f][j]; }
			for( int f=0; f<3; f++ ){ output_file_BCs <<std::setprecision(csv_precision)<<","<< B_rel_diff_inner  [f]    <<std::setprecision(csv_precision)<<","<< E_rel_diff_inner  [f]   ; }
			for( int f=0; f<3; f++ ){ output_file_BCs <<std::setprecision(csv_precision)<<","<< B_outer_before_BCs[f][j] <<std::setprecision(csv_precision)<<","<< E_outer_before_BCs[f][j]; }
			for( int f=0; f<3; f++ ){ output_file_BCs <<std::setprecision(csv_precision)<<","<< B_outer_after_BCs [f][j] <<std::setprecision(csv_precision)<<","<< E_outer_after_BCs [f][j]; }
			for( int f=0; f<3; f++ ){ output_file_BCs <<std::setprecision(csv_precision)<<","<< B_rel_diff_outer  [f]    <<std::setprecision(csv_precision)<<","<< E_rel_diff_outer  [f]   ; }
			
			output_file_BCs << "\n";
		}
	}
	
}




void output_headers_to_screen(){
	std::cout <<std::left<<std::setw(w)<< "T_index" <<std::left<<std::setw(w)<< "T" <<"|\t"
	          <<std::left<<std::setw(w)<< "B_r" <<std::left<<std::setw(w)<< "B_t" <<std::left<<std::setw(w)<< "B_p" <<"|\t"
			  <<std::left<<std::setw(w)<< "E_r" <<std::left<<std::setw(w)<< "E_t" <<std::left<<std::setw(w)<< "E_p" <<"|\t"
			  <<std::left<<std::setw(w)<< "% E pts chngd" <<"|\t"
			  <<std::left<<std::setw(w)<< "nans_tot" <<std::left<<std::setw(w)<< "% done" <<std::left<<std::setw(w)<< "Est. time left"
			  << std::endl;
}




void output_to_screen(){
	
	if( ( T_index % cout_freq_T == 0 ) or ( T_index == n_T ) ){
	
		double percent_done = round( (double)1000*T_index/n_T )/10; // This format makes it exactly one decimal place.
		
		std::cout <<std::left<<std::setw(w)<< T_index <<std::left<<std::setw(w)<< T[T_index] <<"|\t"
				  <<std::left<<std::setw(w)<< B[0][cout_i][cout_j] <<std::left<<std::setw(w)<< B[1][cout_i][cout_j] <<std::left<<std::setw(w)<< B[2][cout_i][cout_j] <<"|\t"
				  <<std::left<<std::setw(w)<< E[0][cout_i][cout_j] <<std::left<<std::setw(w)<< E[1][cout_i][cout_j] <<std::left<<std::setw(w)<< E[2][cout_i][cout_j] <<"|\t"
				  <<std::left<<std::setw(w)<< percentage_E_points_points_changed_for_second_FF_condition <<"|\t"
				  <<std::left<<std::setw(w)<< nans_tot <<std::left<<std::setw(w)<< percent_done;
		
		if( T_index == 0 ){
			std::cout << std::endl;
		} else {
			double time_now_seconds = std::chrono::high_resolution_clock::now().time_since_epoch().count() * 1e-9;
			double runtime = time_now_seconds - time_start_seconds;
			double time_left = ( (double) n_T / T_index - 1.0 ) * runtime;
			int time_left_h = time_left / 3600;
			int time_left_m = time_left / 60 - time_left_h * 60;
			int time_left_s = time_left - time_left_h * 3600 - time_left_m * 60;
			
			if( time_left_m < 10 ){
				std::cout << time_left_h << ":0" << time_left_m;
			} else {
				std::cout << time_left_h << ":"  << time_left_m;
			}
			if( time_left_s < 10 ){
				std::cout << ":0" << time_left_s << std::endl;
			} else {
				std::cout << ":"  << time_left_s << std::endl;
			}
		}
	}
}




void output_parameters_to_screen(){
	// No need to estimate the size of the history CSV since it's far smaller than the profiles CSV.
	// For the future, perhaps split into multiple profiles CSVs. Then, add them all together to get an estimated total space required. A nice feature would be to fail if that much space isn't available, or if it exceeds some specified threshold.
	// NOTE: Number of rotation periods needs updating to include the ramp. Calculate on paper by integrating dOmega_by_dT.
	// Or add a variable rotation_angle which just gets Omega * delta_T added to it at each timestep, and can be recorded in the _history.csv file.
	
	std::cout << log_file_comment << "\n" << std::endl;
	
	std::cout << "Time start:\t" << time_start_string << std::endl;
	
	std::cout << "Parameter values (given in code units, then in SI units)" << std::endl;
	std::cout << "P_final                :\t" <<std::left<<std::setw(w)<< P_final      <<std::left<<std::setw(w)<< P_SI_final     << "\n";
	std::cout << "Omega_final            :\t" <<std::left<<std::setw(w)<< Omega.back() <<std::left<<std::setw(w)<< Omega_SI.back() << std::endl;
	std::cout << "R_LC_max               :\t" <<std::left<<std::setw(w)<< R_LC .back() <<std::left<<std::setw(w)<< R_LC_SI .back() << std::endl;
	std::cout << "r_min, r_max           :\t" <<std::left<<std::setw(w)<< r_min   <<std::left<<std::setw(w)<< r_max              << std::endl;
	std::cout << "n_points_r, delta_r    :\t" <<std::left<<std::setw(w)<< n_points_r <<std::left<<std::setw(w)<< delta_r << std::endl;
	std::cout << "n_points_t, delta_t    :\t" <<std::left<<std::setw(w)<< n_points_t <<std::left<<std::setw(w)<< delta_t << std::endl;
	std::cout << "Timestep               :\t" <<std::left<<std::setw(w)<< delta_T <<std::left<<std::setw(w)<< delta_T * tau << std::endl;
	std::cout << "CFL max timestep       :\t" <<std::left<<std::setw(w)<< delta_T_CFL <<std::left<<std::setw(w)<< delta_T_CFL * tau << std::endl;
	std::cout << "Length of simulation   :\t" <<std::left<<std::setw(w)<< T_max << T_max * tau <<"\tsteps: " << n_T << std::endl;
	std::cout << "Est. CSV size (MB)     :\t" << estimate_csv_profiles_filesize() << std::endl;
	std::cout << "ell_max                :\t" << ell_max << std::endl;
	std::cout << "Printing values at ( r[" << cout_i <<"], theta[" << cout_j <<"] ) = ( " << r[cout_i] <<", " << t[cout_j] << " ) = ( " << r[cout_i] <<", " << t[cout_j]/pi << " pi )." << std::endl;
	std::cout << "dOmega_by_dT_index     :\t" << dOmega_by_dT_index <<"\tdOmega_by_dT:\t"<< dOmega_by_dT << std::endl;
	std::cout << "Ramp start             :\t" << T_index_rotation_ramp_start <<"\t"<< T_index_rotation_ramp_start * delta_T << std::endl;
	std::cout << "Ramp stop              :\t" << T_index_rotation_ramp_stop  <<"\t"<< T_index_rotation_ramp_stop  * delta_T << std::endl;
	
	std::cout << "\nuse_outer_sponge_layer:\t" << use_outer_sponge_layer << std::endl;
	std::cout << "sigma_0, gamma, beta    :\t" << friction_sigma_0 <<"\t"<< friction_gamma <<"\t"<< friction_beta << std::endl;
	
	std::cout << "\nCSV file saved:\t" << output_filename_profiles   << std::endl;
	std::cout <<   "CSV file saved:\t" << output_filename_history    << std::endl;
	std::cout <<   "CSV file saved:\t" << output_filename_VSH_coeffs << std::endl;
	std::cout <<   "CSV file saved:\t" << output_filename_BCs        << std::endl;
	std::cout <<   "Log file saved:\t" << output_filename_log        << std::endl;
	
	std::cout << "\nTo cut evolution short, press CTRL+C at any time. This is perfectly safe and the CSV file will still be saved up to that point.\n" << std::endl;
}




void output_parameters_to_log_file(){
	
	output_file_log << log_file_comment << "\n\n";
	
	output_file_log << "Time start:\t" << time_start_string << "\n";
	
	output_file_log << "Parameter values (given in code units, then in SI units)\n";
	output_file_log << "P_final                :\t" <<std::left<<std::setw(w)<< P_final     <<std::left<<std::setw(w)<< P_SI_final     << "\n";
	output_file_log << "Omega_final            :\t" <<std::left<<std::setw(w)<< Omega.back() <<std::left<<std::setw(w)<< Omega_SI.back() << "\n";
	output_file_log << "R_LC_max               :\t" <<std::left<<std::setw(w)<< R_LC .back() <<std::left<<std::setw(w)<< R_LC_SI .back() << "\n";
	output_file_log << "r_min, r_max           :\t" <<std::left<<std::setw(w)<< r_min   <<std::left<<std::setw(w)<< r_max << "\n";
	output_file_log << "n_points_r, delta_r    :\t" <<std::left<<std::setw(w)<< n_points_r <<std::left<<std::setw(w)<< delta_r << std::endl;
	output_file_log << "n_points_t, delta_t    :\t" <<std::left<<std::setw(w)<< n_points_t <<std::left<<std::setw(w)<< delta_t << std::endl;
	output_file_log << "Timestep               :\t" <<std::left<<std::setw(w)<< delta_T <<std::left<<std::setw(w)<< delta_T * tau << "\n";
	output_file_log << "CFL max timestep       :\t" <<std::left<<std::setw(w)<< delta_T_CFL <<std::left<<std::setw(w)<< delta_T_CFL * tau << "\n";
	output_file_log << "Length of simulation   :\t" <<std::left<<std::setw(w)<< T_max << T_max * tau <<"\tsteps: " << n_T << "\n";
	output_file_log << "Est. CSV size (MB)     :\t" << estimate_csv_profiles_filesize() << "\n";
	output_file_log << "ell_max                :\t" << ell_max << "\n";
	output_file_log << "Printing values at ( r[" << cout_i <<"], theta[" << cout_j <<"] ) = ( " << r[cout_i] <<", " << t[cout_j] << " ) = ( " << r[cout_i] <<", " << t[cout_j]/pi << " pi )." << "\n";
	output_file_log << "dOmega_by_dT_index     :\t" << dOmega_by_dT_index <<"\tdOmega_by_dT:\t"<< dOmega_by_dT << "\n";
	output_file_log << "Ramp start             :\t" << T_index_rotation_ramp_start <<"\t"<< T_index_rotation_ramp_start * delta_T << "\n";
	output_file_log << "Ramp stop              :\t" << T_index_rotation_ramp_stop  <<"\t"<< T_index_rotation_ramp_stop  * delta_T << "\n";
	
	output_file_log << "\nuse_outer_sponge_layer:\t" << use_outer_sponge_layer << "\n";
	output_file_log << "sigma_0, gamma, beta    :\t" << friction_sigma_0 <<"\t"<< friction_gamma <<"\t"<< friction_beta << "\n";
	
	output_file_log << "\nCSV file saved:\t" << output_filename_profiles   << "\n";
	output_file_log <<   "CSV file saved:\t" << output_filename_history    << "\n";
	output_file_log <<   "CSV file saved:\t" << output_filename_VSH_coeffs << "\n";
	output_file_log <<   "CSV file saved:\t" << output_filename_BCs        << "\n";
	output_file_log <<   "Log file saved:\t" << output_filename_log        << "\n";
}




void output_execution_time(){
	
	double time_stop_seconds = std::chrono::high_resolution_clock::now().time_since_epoch().count() * 1e-9;
	time_t time_stop_integer = time( 0 );
	char*  time_stop_string  = ctime( &time_stop_integer );

	double exec_time   = time_stop_seconds - time_start_seconds;
	int    exec_time_h = exec_time / 3600;
	int    exec_time_m = exec_time / 60 - exec_time_h * 60;
	int    exec_time_s = exec_time - exec_time_h * 3600 - exec_time_m * 60;
	
	std::cout       << "\nCode finished successfully.\nTime stop:\t" << time_stop_string << "Execution time (h:mm:ss):\t" << exec_time_h;
	output_file_log << "\nCode finished successfully.\nTime stop:\t" << time_stop_string << "Execution time (h:mm:ss):\t" << exec_time_h;
	
	if( exec_time_m < 10 ){
		std::cout       << ":0" << exec_time_m;
		output_file_log << ":0" << exec_time_m;
	} else {
		std::cout       << ":"  << exec_time_m;
		output_file_log << ":"  << exec_time_m;
	}
	if( exec_time_s < 10 ){
		std::cout       << ":0" << exec_time_s << std::endl;
		output_file_log << ":0" << exec_time_s << std::endl;
	} else {
		std::cout       << ":"  << exec_time_s << std::endl;
		output_file_log << ":"  << exec_time_s << std::endl;
	}
	
}
