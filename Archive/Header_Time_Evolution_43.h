/*
Header_Time_Evolution_43.h

Code snippets for the Time_Evolution_xx.cpp codes. Keeps them shorter and allows for two separate for-loops with Euler and Adams-Bashforth integration
without copying code.

V43: Remove old versions of radial derivative and numerical integration functions.
     Remove fields CSV read and write because they take up too much memory.
	 Remove Chebyshev decomposition functions.
*/

//----- Constants -----
const double pi = acos( -1 );
const double c    = 299792458.0;			// Speed of light in m s^-1.
const double mu_0 = 1.25663706212e-6;		// Permeability of free space in kg m s^-2 A^-2.

//----- Factors to convert between SI units and code units -----
const double tau      = R_star / c;
const double B_factor = mu_0/(4.0*pi) * m_star * pow( R_star, -3.0 );
const double E_factor = B_factor / c;

//----- Global variables -----
const double T_max              = T_index_max * delta_T;
const double P_final            = P_SI_final / tau;				// Final rotation period  in code units.
const double Omega_final        = 2.0 * pi / P_final;			// Final angular velocity in code units.
const double dOmega_by_dT_index = Omega_final / ( T_index_rotation_ramp_stop - T_index_rotation_ramp_start + 1.0 );
double       dOmega_by_dT       = dOmega_by_dT_index / delta_T;

double       delta_r     = ( r_max - r_min ) / ( (double) n_r - 1.0 );			// Linear grid spacing in radial direction.
double       delta_t     = pi / ( (double) n_t - 1.0 );							// Linear grid spacing in polar  direction.
double       delta_T_CFL = 0;													// Maximum timestep allowed by the CFL condition. Calculated by void calculate_CFL_max_timestep().
const double t_min       = 0.5 * delta_t;										// Minimum theta. Don't expect this to change from zero but keep as variable. Offset by half a gridpoint to avoid evolving at the poles.
const double t_max       = 0.5 * delta_t + pi;													// Maximum theta. Don't expect this to change from pi   but keep as variable.


//--- Output files -----
std::ofstream output_file_time_values;
std::ofstream output_file_gridpoints;
std::ofstream output_file_profiles;
std::ofstream output_file_history;
std::ofstream output_file_VSH_coeffs;
std::ofstream output_file_BCs;
std::ofstream output_file_log;

std::string output_filename_time_values = "CSV/"   + output_filename + "_1_time_values.csv";
std::string output_filename_gridpoints  = "CSV/"   + output_filename + "_2_gridpoints.csv";
std::string output_filename_profiles    = "CSV/"   + output_filename + "_4_profiles.csv";
std::string output_filename_history     = "CSV/"   + output_filename + "_5_history.csv";
std::string output_filename_VSH_coeffs  = "CSV/"   + output_filename + "_6_VSH_coeffs.csv";
std::string output_filename_BCs         = "CSV/"   + output_filename + "_7_BCs.csv";
std::string output_filename_log         = "Logs/"  + output_filename + ".txt";


//----- Coordinates and arrays -----
std::vector<double> r  ( n_r );	// Radial coordinate r.
std::vector<double> t  ( n_t );	// Polar coordinate theta.
std::vector<double> st ( n_t );	// sin(theta).
std::vector<double> ct ( n_t );	// cos(theta).

std::vector<double> T                   ( T_index_max );	// Time in code units.
std::vector<double> T_SI                ( T_index_max );	// Time in SI   units.
std::vector<double> Omega               ( T_index_max );	// Angular velocity in code units.
std::vector<double> Omega_SI            ( T_index_max );	// Angular velocity in SI   units.
std::vector<double> P                   ( T_index_max );	// Rotation period in code units.
std::vector<double> P_SI                ( T_index_max );	// Rotation period in SI   units.
std::vector<double> R_LC                ( T_index_max );	// Light cylinder radius in code units.
std::vector<double> R_LC_SI             ( T_index_max );	// Light cylinder radius in SI units.
std::vector<double> star_rotation_angle ( T_index_max );	// Total angle through which the star has rotated, in radians.

std::vector< std::vector<double> > x ( n_r, std::vector<double> ( n_t ) );
std::vector< std::vector<double> > z ( n_r, std::vector<double> ( n_t ) );

std::vector< std::vector<double> > P0 ( n_t, std::vector<double> ( ell_max+1 ) );
std::vector< std::vector<double> > P1 ( n_t, std::vector<double> ( ell_max+1 ) );

std::vector<double> sqrt_2Lplus1_over_4pi        ( ell_max+1 );
std::vector<double> sqrt_2Lplus1_pi              ( ell_max+1 );
std::vector<double> sqrt_2Lplus1_pi_over_L_Lplus1( ell_max+1 );

std::vector< std::vector< std::vector<double> > > B_VSH_coeffs     ( 3, std::vector< std::vector<double> > ( ell_max+1, std::vector<double> ( n_r ) ) );
std::vector< std::vector< std::vector<double> > > B_VSH_coeffs_dr  ( 3, std::vector< std::vector<double> > ( ell_max+1, std::vector<double> ( n_r ) ) );
std::vector< std::vector< std::vector<double> > > B_VSH_coeffs_dr2 ( 3, std::vector< std::vector<double> > ( ell_max+1, std::vector<double> ( n_r ) ) );
std::vector< std::vector< std::vector<double> > > E_VSH_coeffs     ( 3, std::vector< std::vector<double> > ( ell_max+1, std::vector<double> ( n_r ) ) );
std::vector< std::vector< std::vector<double> > > E_VSH_coeffs_dr  ( 3, std::vector< std::vector<double> > ( ell_max+1, std::vector<double> ( n_r ) ) );
std::vector< std::vector< std::vector<double> > > E_VSH_coeffs_dr2 ( 3, std::vector< std::vector<double> > ( ell_max+1, std::vector<double> ( n_r ) ) );

std::vector< std::vector< std::vector<double> > > B                    ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );
std::vector< std::vector< std::vector<double> > > B_minus1             ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );
std::vector< std::vector< std::vector<double> > > B_minus2             ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );
std::vector< std::vector< std::vector<double> > > B_minus3             ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );
std::vector< std::vector< std::vector<double> > > B_dT                 ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );
std::vector< std::vector< std::vector<double> > > B_dT_minus1          ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );
std::vector< std::vector< std::vector<double> > > B_dT_minus2          ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );
std::vector< std::vector< std::vector<double> > > B_dT_minus3          ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );
std::vector< std::vector< std::vector<double> > > B_curl               ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );
std::vector< std::vector< std::vector<double> > > B_Cartesian          ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );
std::vector< std::vector< std::vector<double> > > B_dr                 ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );
std::vector< std::vector< std::vector<double> > > B_dr2                ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );
std::vector< std::vector< std::vector<double> > > B_series             ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );	// Values of the VSH decomposition of the field components.
std::vector< std::vector< std::vector<double> > > B_abs_rel_dev        ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );
std::vector< std::vector< std::vector<double> > > B_from_VSH_coeffs    ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );
std::vector< std::vector< std::vector<double> > > B_dT_before_sponge   ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );
std::vector< std::vector< std::vector<double> > > B_dT_after_sponge    ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );
std::vector< std::vector< std::vector<double> > > B_dT_sponge_rel_diff ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );

std::vector< std::vector< std::vector<double> > > E                    ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );
std::vector< std::vector< std::vector<double> > > E_minus1             ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );
std::vector< std::vector< std::vector<double> > > E_minus2             ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );
std::vector< std::vector< std::vector<double> > > E_minus3             ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );
std::vector< std::vector< std::vector<double> > > E_dT                 ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );
std::vector< std::vector< std::vector<double> > > E_dT_minus1          ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );
std::vector< std::vector< std::vector<double> > > E_dT_minus2          ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );
std::vector< std::vector< std::vector<double> > > E_dT_minus3          ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );
std::vector< std::vector< std::vector<double> > > E_curl               ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );
std::vector< std::vector< std::vector<double> > > E_Cartesian          ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );
std::vector< std::vector< std::vector<double> > > E_dr                 ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );
std::vector< std::vector< std::vector<double> > > E_dr2                ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );
std::vector< std::vector< std::vector<double> > > E_series             ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );	// Values of the VSH decomposition of the field components.
std::vector< std::vector< std::vector<double> > > E_abs_rel_dev        ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );
std::vector< std::vector< std::vector<double> > > E_from_VSH_coeffs    ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );
std::vector< std::vector< std::vector<double> > > E_dT_before_sponge   ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );
std::vector< std::vector< std::vector<double> > > E_dT_after_sponge    ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );
std::vector< std::vector< std::vector<double> > > E_dT_sponge_rel_diff ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );

std::vector< std::vector< std::vector<double> > > J                 ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );	// Current density vector.
std::vector< std::vector< std::vector<double> > > J_Cartesian       ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );


std::vector< std::vector<double> > B_squared                     ( n_r, std::vector<double> ( n_t ) );
std::vector< std::vector<double> > E_squared                     ( n_r, std::vector<double> ( n_t ) );
std::vector< std::vector<double> > B_squared_minus_E_squared     ( n_r, std::vector<double> ( n_t ) );
std::vector< std::vector<double> > E_dot_B                       ( n_r, std::vector<double> ( n_t ) );
std::vector< std::vector<double> > B_div                         ( n_r, std::vector<double> ( n_t ) );
std::vector< std::vector<double> > E_div                         ( n_r, std::vector<double> ( n_t ) );
std::vector< std::vector<double> > normalised_B_div              ( n_r, std::vector<double> ( n_t ) );		// r * div(B) / B, for a volume integral to characterise the size of div(B) over the domain.
std::vector< std::vector<bool  > > force_free_condition_1        ( n_r, std::vector<bool  > ( n_t ) );		// Is it true at each gridpoint that: E dot B = 0
std::vector< std::vector<bool  > > force_free_condition_2        ( n_r, std::vector<bool  > ( n_t ) );		// Is it true at each gridpoint that: B^2 - E^2 > 0
std::vector< std::vector<bool  > > force_free_conditions         ( n_r, std::vector<bool  > ( n_t ) );		// Combine the first two conditions to give one overall condition for easy plotting or other diagnostics.
std::vector< std::vector<int   > > is_gridpoint_within_FF_region ( n_r, std::vector<int   > ( n_t ) );
std::vector< std::vector<int   > > were_FF_conditions_applied    ( n_r, std::vector<int   > ( n_t ) );		// +1: FF conditions applied here. -1: FF conditions not applied here. 0: Outside FF region.

std::vector<double> finite_line_element_dt ( n_t );	// Finite line elements dtheta for integral over theta.
std::vector< std::vector<double> > finite_area_element_for_r_t ( n_r, std::vector<double> ( n_t ) );	// Finite area elements r dr dtheta for 2D integral over r and theta.

std::vector<double> B_stdev          ( 3 );	// Standard devation of the VSH decomposition of each vector component.
std::vector<double> B_max_abs_rel_dev( 3 );	// Maximum absolute relative deviation between the VSH decomposition and the exact value of each vector component.
std::vector<double> E_stdev          ( 3 );	// Standard devation of the VSH decomposition of each vector component.
std::vector<double> E_max_abs_rel_dev( 3 );	// Maximum absolute relative deviation between the VSH decomposition and the exact value of each vector component.

std::vector<double> B_t_dt_inner ( n_t ); // 20231214 Specify d/dT B at inner boundary; requires expression for d/dt B_t.

double total_electric_energy            = 0;
double total_magnetic_energy            = 0;
double volume_of_domain                 = 0;
double normalised_B_div_volume_integral = 0;	// Measure of the strength of div(B) across the domain, to track how div(B) evolves with time.

std::vector<int> B_nans ( 6 );	// Count the number of "nan" or "inf" values of each field component at a given timestep.
std::vector<int> E_nans ( 6 );	// Count the number of "nan" or "inf" values of each field component at a given timestep.
int nans_tot = 0;

std::vector<double> friction ( n_r );		// Sponge layer at outer boundary to absorb outgoing waves (see Parfrey, 2012, PhD thesis, S3.9).

std::vector<double> VSH_integrand_r ( n_t );
std::vector<double> VSH_integrand_1 ( n_t );
std::vector<double> VSH_integrand_2 ( n_t );

double time_start_seconds = std::chrono::high_resolution_clock::now().time_since_epoch().count() * 1e-9;
time_t time_start_integer = time( 0 );
char*  time_start_string  = ctime( &time_start_integer );

int    T_index = 0;

int    n_gridpoints_within_FF_region                              = 0;
int    num_E_points_changed_for_second_FF_condition               = 0;	// Output to CSV to keep track of how active the FF condition application is.
double percentage_E_points_points_changed_for_second_FF_condition = 0;

std::vector< std::vector<double> > B_inner_before_BCs     ( 3, std::vector<double> ( n_t ) );
std::vector< std::vector<double> > B_inner_after_BCs      ( 3, std::vector<double> ( n_t ) );
std::vector< std::vector<double> > B_outer_before_BCs     ( 3, std::vector<double> ( n_t ) );
std::vector< std::vector<double> > B_outer_after_BCs      ( 3, std::vector<double> ( n_t ) );
std::vector< std::vector<double> > E_inner_before_BCs     ( 3, std::vector<double> ( n_t ) );
std::vector< std::vector<double> > E_inner_after_BCs      ( 3, std::vector<double> ( n_t ) );
std::vector< std::vector<double> > E_outer_before_BCs     ( 3, std::vector<double> ( n_t ) );
std::vector< std::vector<double> > E_outer_after_BCs      ( 3, std::vector<double> ( n_t ) );
std::vector< std::vector<double> > B_inner_BCs_difference ( 3, std::vector<double> ( n_t ) );
std::vector< std::vector<double> > B_outer_BCs_difference ( 3, std::vector<double> ( n_t ) );
std::vector< std::vector<double> > E_inner_BCs_difference ( 3, std::vector<double> ( n_t ) );
std::vector< std::vector<double> > E_outer_BCs_difference ( 3, std::vector<double> ( n_t ) );




//----- Functions -----

double next_sum_of_reciprocal_powers_of_2( double x ){
	// Finds the first number higher than x that can be represented as a sum of reciprocal powers of 2, and therefore can be expressed exactly.
	// Saves floating point errors when x is a number repeatedly added.
	double delta = pow( 2, -14.0 );
	double ret = 0;
	while( ret <= x ){
		ret += delta;
	}
	return ret;
}




void calculate_time_and_rotation_values(){
	
	//--- Set delta_T to a sum of reciprocal powers of 2 ---
	double delta_T_old = delta_T;
	delta_T = next_sum_of_reciprocal_powers_of_2( delta_T );
	std::cout       << "Setting delta_T to a sum of reciprocal powers of 2 to avoid floating point errors:\t" <<std::setprecision(precision_csv)<< delta_T <<" from "<< delta_T_old << std::endl;
	output_file_log << "Setting delta_T to a sum of reciprocal powers of 2 to avoid floating point errors:\t" <<std::setprecision(precision_csv)<< delta_T <<" from "<< delta_T_old << "\n";
	
	//--- Time in code units ---
	for( int T_index=0; T_index<T_index_max; T_index++ ){	
		T[T_index] = (double) T_index * delta_T;
	}
	
	//--- Angular velocity in code units ---
	for( int T_index=0; T_index<T_index_max; T_index++ ){
		if( T_index < T_index_rotation_ramp_start ){
			Omega[T_index] = 0;
		}
		if( ( T_index >= T_index_rotation_ramp_start ) and ( T_index <= T_index_rotation_ramp_stop ) ){
			Omega[T_index] = Omega[ T_index - 1 ] + dOmega_by_dT_index;
		}
		if( T_index > T_index_rotation_ramp_stop ){
			Omega[T_index] = Omega[ T_index - 1 ];
		}
	}
	
	//--- Rotation period and light cylinder radius ---
	for( int T_index=0; T_index<T_index_max; T_index++ ){
		if( Omega[T_index] != 0 ){
			P   [T_index] = 2.0 * pi / Omega[T_index];
			R_LC[T_index] = 1.0      / Omega[T_index];
		}
	}
	
	//--- Total angle through which star has rotated in radians ---
	for( int T_index=T_index_rotation_ramp_start; T_index<T_index_max; T_index++ ){
		star_rotation_angle[T_index] = star_rotation_angle[ T_index - 1 ] + Omega[T_index] * delta_T;
	}
	
	//--- SI units ---
	for( int T_index=0; T_index<T_index_max; T_index++ ){
		T_SI    [T_index] = T    [T_index] * tau;
		Omega_SI[T_index] = Omega[T_index] / tau;
		P_SI    [T_index] = P    [T_index] * tau;
		R_LC_SI [T_index] = R_LC [T_index] * R_star;
	}
	
}



void calculate_gridpoints(){
	
	//--- Set delta_r to a sum of reciprocal powers of 2 ---
	double delta_r_old = delta_r;
	double r_max_old   = r_max;
	delta_r = next_sum_of_reciprocal_powers_of_2( delta_r );
	r_max   = r_min + ( n_r - 1 ) * delta_r;
	std::cout       << "Setting delta_r to a sum of reciprocal powers of 2 to avoid floating point errors:\t" <<std::setprecision(precision_csv)<< delta_r <<" from "<< delta_r_old << std::endl;
	output_file_log << "Setting delta_r to a sum of reciprocal powers of 2 to avoid floating point errors:\t" <<std::setprecision(precision_csv)<< delta_r <<" from "<< delta_r_old << "\n";
	std::cout       << "New r_max = " << r_max <<" from " << r_max_old << "\n" << std::endl;
	output_file_log << "New r_max = " << r_max <<" from " << r_max_old << "\n\n";
	
	
	//--- Calculate gridpoints ---
	for( int i=0; i<n_r; i++ ){
		r[i] = r_min + i * delta_r;
	}
	
	for( int j=0; j<n_t; j++ ){
		t[j] = t_min + j * delta_t;
	}
	
}




void calculate_functions_of_gridpoints(){
	// Keep separate to void calculate_gridpoints() so that this can be called when reloading previous runs.
	
	//--- sin(theta) and cos(theta) ---
	for( int j=0; j<n_t; j++ ){
		st[j] = sin( t[j] );
		ct[j] = cos( t[j] );
	}
	
	//--- Cartesian coordinates ---
	for( int i=0; i<n_r; i++ ){
		for( int j=0; j<n_t; j++ ){
			x[i][j] = r[i] * st[j];
			z[i][j] = r[i] * ct[j];
		}
	}
	
	//--- Finite line and area elements for numerical integration ---
	for( int j=0; j<=n_t-2; j++ ){
		finite_line_element_dt[j] = ct[j] - ct[j+1];
	}
	
	for( int i=0; i<=n_r-2; i++ ){
		for( int j=0; j<=n_t-2; j++ ){
			finite_area_element_for_r_t[i][j] = delta_r * ( r[i]*r[i] + r[i]*delta_r + (1.0/3.0)*delta_r*delta_r ) * ( ct[j] - ct[j+1] );
		}
	}
	
	//--- Volume of domain ---
	volume_of_domain = (4.0*pi/3.0) * ( pow(r_max,3.0) - pow(r_min,3.0) );
	
}




void calculate_associated_legendre_functions(){
	
	for( int j=0; j<n_t; j++ ){
		
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




void calculate_CFL_max_timestep(){
	// This is the approximation for constant grid spacing in spherical coordinates described in the writeup.
	// Replace by an exact expression if needed, possibly using functions from Test codes/Test_CFL_condition.cpp.
	delta_T_CFL = std::min( r_min * delta_t, delta_r );
}




void calculate_friction_for_sponge_layer(){
	// Parfrey, 2012, PhD thesis, \S3.9.2.
	
	if( use_outer_sponge_layer ){
		for( int i=0; i<n_r; i++ ){
			friction[i] = pow( ( r[i] - r_sponge ) / ( r_max - r_sponge ), friction_beta );
			friction[i] = friction_sigma_0 * ( 1.0 - exp( - friction_gamma * friction[i] ) );
		}
	}
	
}




double integral_Booles_rule( std::vector<double> f, double h ){
	// 1D integral of f(u) in Cartesian coordinates, given pre-calculated list of function values. Only works if the number of gridpoints is one more than a multiple of 4.
	int N = f.size();
	
	if( N % 4 != 1 ){
		std::cout << "This method is most accurate when you have 4m+1 data points, where m is any natual number. Consider a different method or adding one or two extra points." << std::endl;
	}
	
	double ret = 0;
	
	for( int k=1; k<=N/4; k++ ){
		ret += 7.0 * ( f[4*k-4] + f[4*k] ) + 32.0 * ( f[4*k-3] + f[4*k-1] ) + 12.0 * f[4*k-2];
	}
	
	return ret * h * 2.0 / 45.0;
	
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
	for( int i=0; i<n_r; i++ ){
		for( int j=0; j<n_t; j++ ){
			B[0][i][j] = B_r_function( r[i], t[j] );
			B[1][i][j] = B_t_function( r[i], t[j] );
			B[2][i][j] = B_p_function( r[i], t[j] );
		}
	}
}




void apply_inner_boundary_conditions(){
	// Petri 2012, \S3.3.
	
	//--- Reset arrays ---
	B_inner_before_BCs = std::vector< std::vector<double> > ( 3, std::vector<double> ( n_t ) );
	B_inner_after_BCs  = std::vector< std::vector<double> > ( 3, std::vector<double> ( n_t ) );
	E_inner_before_BCs = std::vector< std::vector<double> > ( 3, std::vector<double> ( n_t ) );
	E_inner_after_BCs  = std::vector< std::vector<double> > ( 3, std::vector<double> ( n_t ) );
	
	//--- Record values before inner BCs applied ---
	for( int j=0; j<n_t; j++ ){
		for( int f=0; f<3; f++ ){
			B_inner_before_BCs[f][j] = B[f][0][j];
			E_inner_before_BCs[f][j] = E[f][0][j];
		}
	}
	
	//--- Apply inner BCs ---
	for( int j=0; j<n_t; j++ ){
		B[0][0][j] = B_r_function( r_min, t[j] );
		//B[2][0][j] = 2.0 * delta_T * Omega * st[j] * cos(t[j]);					// SPECIFIC TO A DIPOLE AND NOT SURE IT'S THEORETICALLY GROUNDED; JUST SEEMS TO FIT THE DATA.
		E[1][0][j] = B_r_function( r_min, t[j] ) * - Omega[T_index] * r_min * st[j];
		E[2][0][j] = 0;
	}
	
	//B[2][0][ (int) n_t / 2 - 1 ] = 0;	// Set B_phi ( R_NS, pi/2 ) = 0 to avoid evolution of noise. Suggestion from Sam 20231208.
	
	//--- Record values after inner BCs applied ---
	for( int j=0; j<n_t; j++ ){
		for( int f=0; f<3; f++ ){
			B_inner_after_BCs[f][j] = B[f][0][j];
			E_inner_after_BCs[f][j] = E[f][0][j];
			B_inner_BCs_difference[f][j] = B_inner_after_BCs[f][j] - B_inner_before_BCs[f][j];
			E_inner_BCs_difference[f][j] = E_inner_after_BCs[f][j] - E_inner_before_BCs[f][j];
		}
	}

}




void apply_outer_boundary_conditions(){
	// Petri 2012, \S3.3.
	
	//--- Reset arrays ---
	B_outer_before_BCs = std::vector< std::vector<double> > ( 3, std::vector<double> ( n_t ) );
	B_outer_after_BCs  = std::vector< std::vector<double> > ( 3, std::vector<double> ( n_t ) );
	E_outer_before_BCs = std::vector< std::vector<double> > ( 3, std::vector<double> ( n_t ) );
	E_outer_after_BCs  = std::vector< std::vector<double> > ( 3, std::vector<double> ( n_t ) );
	
	//--- Record values before outer BCs applied ---
	for( int j=0; j<n_t; j++ ){
		for( int f=0; f<3; f++ ){
			B_outer_before_BCs[f][j] = B[f].back()[j];
			E_outer_before_BCs[f][j] = E[f].back()[j];
		}
	}
	
	//--- Apply outer BCs ---
	if( T_index < T_index_rotation_ramp_start ){
		for( int j=0; j<n_t; j++ ){
			B[1].back()[j] = B_t_function( r_max, t[j] );
			B[2].back()[j] = B_p_function( r_max, t[j] );
			E[1].back()[j] = 0;
			E[2].back()[j] = 0;
		}
	}
	else {
		for( int j=0; j<n_t; j++ ){
			
			double B_t_PDE = B[1].back()[j];
			double B_p_PDE = B[2].back()[j];
			double E_t_PDE = E[1].back()[j];
			double E_p_PDE = E[2].back()[j];
			
			B[1].back()[j] = 0.5 * ( B_t_PDE - E_p_PDE );
			B[2].back()[j] = 0.5 * ( B_p_PDE + E_t_PDE );
			E[1].back()[j] = 0.5 * ( E_t_PDE + B_p_PDE );
			E[2].back()[j] = 0.5 * ( E_p_PDE - B_t_PDE );
			
		}
	}
	
	//--- Record values after outer BCs applied ---
	for( int j=0; j<n_t; j++ ){
		for( int f=0; f<3; f++ ){
			B_outer_after_BCs[f][j] = B[f].back()[j];
			E_outer_after_BCs[f][j] = E[f].back()[j];
			B_outer_BCs_difference[f][j] = B_outer_after_BCs[f][j] - B_outer_before_BCs[f][j];
			E_outer_BCs_difference[f][j] = E_outer_after_BCs[f][j] - E_outer_before_BCs[f][j];
		}
	}
	
}




std::vector< std::vector<double> > VSH_decomposition_r( std::vector< std::vector< std::vector<double> > > &v ){
	// Calculate the A^{r,ell} coefficients of the VSH series of a vector.
	
	std::vector<double> VSH_integrand( n_t );
	std::vector< std::vector<double> > v_VSH_coeff ( ell_max+1, std::vector<double> ( n_r ) );
	
	for( int i=0; i<n_r; i++ ){
		for( int ell=0; ell<=ell_max; ell++ ){
				
			for( int j=0; j<n_t; j++ ){
				VSH_integrand[j] = v[0][i][j] * P0[j][ell] * st[j];
			}
			
			v_VSH_coeff[ell][i] = integral_Booles_rule( VSH_integrand, delta_t ) * sqrt_2Lplus1_pi[ell];
			
		}
	}
	
	return v_VSH_coeff;
}




std::vector< std::vector<double> > VSH_decomposition_1( std::vector< std::vector< std::vector<double> > > &v ){
	// Calculate the A^{(1),ell} coefficients of the VSH series of a vector.
	
	std::vector<double> VSH_integrand( n_t );
	std::vector< std::vector<double> > v_VSH_coeff ( ell_max+1, std::vector<double> ( n_r ) );
	
	for( int i=0; i<n_r; i++ ){
		for( int ell=1; ell<=ell_max; ell++ ){
				
			for( int j=0; j<n_t; j++ ){
				VSH_integrand[j] = v[1][i][j] * P1[j][ell] * st[j];
			}
			
			v_VSH_coeff[ell][i] = integral_Booles_rule( VSH_integrand, delta_t ) * sqrt_2Lplus1_pi_over_L_Lplus1[ell];
			
		}
	}
	
	return v_VSH_coeff;
}




std::vector< std::vector<double> > VSH_decomposition_2( std::vector< std::vector< std::vector<double> > > &v ){
	// Calculate the A^{(2),ell} coefficients of the VSH series of a vector.
	
	std::vector<double> VSH_integrand( n_t );
	std::vector< std::vector<double> > v_VSH_coeff ( ell_max+1, std::vector<double> ( n_r ) );
	
	for( int i=0; i<n_r; i++ ){
		for( int ell=1; ell<=ell_max; ell++ ){
				
			for( int j=0; j<n_t; j++ ){
				VSH_integrand[j] = v[2][i][j] * P1[j][ell] * st[j];
			}
			
			v_VSH_coeff[ell][i] = integral_Booles_rule( VSH_integrand, delta_t ) * sqrt_2Lplus1_pi_over_L_Lplus1[ell];
			
		}
	}
	
	return v_VSH_coeff;
}




std::vector<double> radial_derivatives_1_FD( std::vector<double> &f ){
	// Calculate the first radial derivatives of a vector f representing a function evaluated at the radial gridpoints r[i] used in this code.
	// Not intended for derivatives with respect to an arbitrary coordinates.
	
	std::vector<double> f_dr ( n_r );
	
	//----- Order-h^6 symmetric finite difference method for the intermediate points; order-h^8 offset forward/backward method for near the endpoints. -----
	
	f_dr[0    ] = ( -0.125*f[8    ] + (8.0/7.0)*f[7    ] - (14.0/3.0)*f[6    ] + 11.2*f[5    ] - 17.5*f[4    ] + (56.0/3.0)*f[3    ] - 14.0*f[2    ] + 8.0*f[1    ] - (761.0/280.0)*f[0    ] ) / delta_r;
	f_dr[n_r-1] = (  0.125*f[n_r-9] - (8.0/7.0)*f[n_r-8] + (14.0/3.0)*f[n_r-7] - 11.2*f[n_r-6] + 17.5*f[n_r-5] - (56.0/3.0)*f[n_r-4] + 14.0*f[n_r-3] - 8.0*f[n_r-2] + (761.0/280.0)*f[n_r-1] ) / delta_r;
	
	f_dr[1    ] = (   (1.0/56.0)*f[8    ] - (1.0/6.0)*f[7    ] +  0.7*f[6    ] - 1.75*f[5    ] +  (35.0/12.0)*f[4    ] -  3.5*f[3    ] +  3.5*f[2    ] - (223.0/140.0)*f[1    ] -  0.125*f[0    ] ) / delta_r;
	f_dr[n_r-2] = (  -(1.0/56.0)*f[n_r-9] + (1.0/6.0)*f[n_r-8] -  0.7*f[n_r-7] + 1.75*f[n_r-6] -  (35.0/12.0)*f[n_r-5] +  3.5*f[n_r-4] -  3.5*f[n_r-3] + (223.0/140.0)*f[n_r-2] +  0.125*f[n_r-1] ) / delta_r;
	
	f_dr[2    ] = (   -(1.0/168.0)*f[8    ] +  (2.0/35.0)*f[7    ] -  0.25*f[6    ] +  (2.0/3.0)*f[5    ] -  1.25*f[4    ] +  2.0*f[3    ] -   0.95*f[2    ] -  (2.0/7.0)*f[1    ] +   (1.0/56.0)*f[0    ] ) / delta_r;
	f_dr[n_r-3] = (    (1.0/168.0)*f[n_r-9] -  (2.0/35.0)*f[n_r-8] +  0.25*f[n_r-7] -  (2.0/3.0)*f[n_r-6] +  1.25*f[n_r-5] -  2.0*f[n_r-4] +   0.95*f[n_r-3] +  (2.0/7.0)*f[n_r-2] -   (1.0/56.0)*f[n_r-1] ) / delta_r;
	
	for( int i=3; i<=n_r-4; i++ ){
		f_dr[i] = ( (1.0/60.0)*(f[i+3]-f[i-3]) - 0.15*(f[i+2]-f[i-2]) + 0.75*(f[i+1]-f[i-1]) ) / delta_r;
	}
	
	return f_dr;
	
}




std::vector<double> radial_derivatives_2_FD( std::vector<double> &f ){
	// Calculate the second radial derivatives of a vector f representing a function evaluated at the radial gridpoints r[i] used in this code.
	// Not intended for derivatives with respect to an arbitrary coordinate.
	
	std::vector<double> f_dr2 ( n_r );
	
	
	
	//----- Order-h^6 symmetric finite difference method for the intermediate points; order-h^8 forward/backward method for near the endpoints. -----
	for( int i=0; i<=3; i++ ){
		f_dr2[i] = ( 3267.0*f[i+8] - 29664.0*f[i+7] + 120008.0*f[i+6] - 284256.0*f[i+5] + 435330.0*f[i+4] - 448672.0*f[i+3] + 312984.0*f[i+2] - 138528.0*f[i+1] + 29531.0*f[i] ) / ( 5040.0*delta_r*delta_r );
	}
	
	for( int i=n_r-4; i<=n_r-1; i++ ){
		f_dr2[i] = ( 3267.0*f[i-8] - 29664.0*f[i-7] + 120008.0*f[i-6] - 284256.0*f[i-5] + 435330.0*f[i-4] - 448672.0*f[i-3] + 312984.0*f[i-2] - 138528.0*f[i-1] + 29531.0*f[i] ) / ( 5040.0*delta_r*delta_r );
	}
	
	for( int i=4; i<=n_r-5; i++ ){
		f_dr2[i] = ( 2.0*(f[i+3]+f[i-3]) - 27.0*(f[i+2]+f[i-2]) + 270.0*(f[i+1]+f[i-1]) - 490.0*f[i] ) / ( 180.0*delta_r*delta_r );
	}
	
	return f_dr2;
}







		
		
void calculate_B_1_ell_and_derivative_divergenceless(){
	// This version of the function is for T>0, at which time we want to update B^{(1),L}.
	
	for( int L=1; L<=ell_max; L++ ){
		for( int i=0; i<n_r; i++ ){
			B_VSH_coeffs   [1][L][i] = ( r[i] * B_VSH_coeffs_dr [0][L][i] + 2.0 * B_VSH_coeffs   [0][L][i] ) / ( (double) L*(L+1.0) );
			B_VSH_coeffs_dr[1][L][i] = ( r[i] * B_VSH_coeffs_dr2[0][L][i] + 3.0 * B_VSH_coeffs_dr[0][L][i] ) / ( (double) L*(L+1.0) );
		}
	}
	
}




std::vector< std::vector< std::vector<double> > > vector_values_from_VSH_coeffs( std::vector< std::vector< std::vector<double> > > &v_VSH_coeffs ){
	// Recalculate the values of a vector at all gridpoints given its VSH series coefficients.
	
	std::vector< std::vector< std::vector<double> > > v_from_VSH_coeffs ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );
	
	for( int i=0; i<n_r; i++ ){
		for( int j=0; j<n_t; j++ ){
			
			//--- The r-coefficients start from ell=0 ---
			for( int ell=0; ell<=ell_max; ell++ ){
				v_from_VSH_coeffs[0][i][j] += v_VSH_coeffs[0][ell][i] * sqrt_2Lplus1_over_4pi[ell] * P0[j][ell];
			}
			
			//--- The (1)- and (2)-coefficients start from ell=1 ---
			for( int ell=1; ell<=ell_max; ell++ ){
				v_from_VSH_coeffs[1][i][j] += v_VSH_coeffs[1][ell][i] * sqrt_2Lplus1_over_4pi[ell] * P1[j][ell];
				v_from_VSH_coeffs[2][i][j] += v_VSH_coeffs[2][ell][i] * sqrt_2Lplus1_over_4pi[ell] * P1[j][ell];
			}
			
		}
	}
	
	return v_from_VSH_coeffs;
}




double stdev_of_array_1d( std::vector<double> exact, std::vector<double> calc ){
	// Calculate standard deviation of an array of calculated results, compared to their known exact values.
	
	double ret = 0;
	
	for( int i=0; i<exact.size(); i++ ){
		ret += pow( exact[i] - calc[i], 2.0 );
	}
	
	return sqrt( ret / (double) exact.size() );
	
}

double stdev_of_array_2d( std::vector< std::vector<double> > exact, std::vector< std::vector<double> > calc ){
	// Calculate standard deviation of an array of calculated results, compared to their known exact values.
	
	double ret = 0;
	
	for( int i=0; i<exact.size(); i++ ){
		for( int j=0; j<exact[0].size(); j++ ){
			ret += pow( exact[i][j] - calc[i][j], 2.0 );
		}
	}
	
	return sqrt( ret / ( (double) exact.size() * exact[0].size() ) );
	
}

std::vector<double> max_deviation_of_array_2d( std::vector< std::vector<double> > exact, std::vector< std::vector<double> > calc ){
	// Calculate maximum absolute deviation between arrays of exact and calcualted results, and return index at which the maximum occurs.
	// Blind to the same maximum occurring more than once; in that case, returns the highest indices at which that happens.
	
	std::vector<double> ret (3);
	
	for( int i=0; i<exact.size(); i++ ){
		for( int j=0; j<exact[0].size(); j++ ){
			
			double abs_dev = abs( exact[i][j] - calc[i][j] );
			
			if( abs_dev > ret[0] ){
				ret[0] = abs_dev;
				ret[1] = i;
				ret[2] = j;
			}
		
		}
	}
	
	return ret;
}




void calculate_stdev_of_VSH_decomposition(){
	
	//--- Reset values to zero ---
	B_abs_rel_dev     = std::vector< std::vector< std::vector<double> > > ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );
	E_abs_rel_dev     = std::vector< std::vector< std::vector<double> > > ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );
	B_stdev           = std::vector<double> ( 3 );
	E_stdev           = std::vector<double> ( 3 );
	B_max_abs_rel_dev = std::vector<double> ( 3 );
	E_max_abs_rel_dev = std::vector<double> ( 3 );
	
	//--- Recalculate vectors from VSH coefficients ---
	B_from_VSH_coeffs = vector_values_from_VSH_coeffs( B_VSH_coeffs );
	E_from_VSH_coeffs = vector_values_from_VSH_coeffs( E_VSH_coeffs );
	
	//--- Calculate stdev ---
	for( int f=0; f<3; f++ ){
		
		for( int i=0; i<n_r; i++ ){
			for( int j=0; j<n_t; j++ ){
				
				B_stdev[f] += pow( B[f][i][j] - B_from_VSH_coeffs[f][i][j], 2.0 );
				E_stdev[f] += pow( E[f][i][j] - E_from_VSH_coeffs[f][i][j], 2.0 );
				
				if( B[f][i][j] != 0 ){
					B_abs_rel_dev[f][i][j] = abs( 1.0 - B_from_VSH_coeffs[f][i][j] / B[f][i][j] );
				}
				if( E[f][i][j] != 0 ){
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
		
		B_stdev[f] = sqrt( B_stdev[f] / ( (double) n_r * n_t ) );
		E_stdev[f] = sqrt( E_stdev[f] / ( (double) n_r * n_t ) );
		
	}
	
}




void calculate_radial_derivatives_of_vector_components(){
	// Not used by the code, but useful for CSV output to check for strange behaviour in dB/dr and dE/dr.
	// To do <20230824: Second derivatives too.
	
	//--- Reset values to zero ---
	B_dr  = std::vector< std::vector< std::vector<double> > > ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );
	E_dr  = std::vector< std::vector< std::vector<double> > > ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );
	B_dr2 = std::vector< std::vector< std::vector<double> > > ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );
	
	//--- Sum over ell for each gridpoint ---
	for( int i=0; i<n_r; i++ ){
		for( int j=0; j<n_t; j++ ){
			
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
	
	//--- Reset values ---
	B_Cartesian = std::vector< std::vector< std::vector<double> > > ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );
	E_Cartesian = std::vector< std::vector< std::vector<double> > > ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );
	J_Cartesian = std::vector< std::vector< std::vector<double> > > ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );
	
	//--- Perform calculation ---
	for( int i=0; i<n_r; i++ ){
		for( int j=0; j<n_t; j++ ){
			
			
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
	
	//--- Reset values ---
	B_squared                 = std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) );
	E_squared                 = std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) );
	B_squared_minus_E_squared = std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) );
	E_dot_B                   = std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) );
	
	//--- Perform calculation ---
	for( int i=0; i<n_r; i++ ){
		for( int j=0; j<n_t; j++ ){
			
			for( int f=0; f<3; f++ ){
				B_squared[i][j] += pow( B[f][i][j], 2.0 );
				E_squared[i][j] += pow( E[f][i][j], 2.0 );
				E_dot_B  [i][j] += E[f][i][j] * B[f][i][j];
			}
			
			B_squared_minus_E_squared[i][j] = B_squared[i][j] - E_squared[i][j];
			
		}
	}
	
}




void ramp_up_electric_fields_for_rotation(){
	if( ( T_index >= T_index_rotation_ramp_start ) and ( T_index <= T_index_rotation_ramp_stop ) ){
	
		int n_steps_rampup = T_index_rotation_ramp_stop - T_index_rotation_ramp_start + 1;
		
		for( int i=0; i<n_r; i++ ){
			for( int j=0; j<n_t; j++ ){
				
				double factor = Omega.back() * r[i] * st[j] / (double) n_steps_rampup;
				
				E[0][i][j] += factor * B_t_function( r[i], t[j] );
				E[1][i][j] -= factor * B_r_function( r[i], t[j] );
				
			}
		}
		
	}
}




void calculate_n_gridpoints_within_FF_region(){
	
	//--- Reset values ---
	is_gridpoint_within_FF_region = std::vector< std::vector<int> > ( n_r, std::vector<int> ( n_t ) );
	n_gridpoints_within_FF_region = 0;
	
	//--- Perform calculation ---
	if( Omega[T_index] > 0 ){
		for( int i=0; i<n_r; i++ ){
			for( int j=0; j<n_t; j++ ){
				if( r[i] * st[j] < R_LC[T_index] ){
					is_gridpoint_within_FF_region[i][j]  = 1;
					n_gridpoints_within_FF_region       += 1;
				}
			}
		}
	}
}




void check_force_free_conditions(){
	
	//--- Reset values ---
	force_free_condition_1 = std::vector< std::vector<bool> > ( n_r, std::vector<bool> ( n_t ) );
	force_free_condition_2 = std::vector< std::vector<bool> > ( n_r, std::vector<bool> ( n_t ) );
	force_free_conditions  = std::vector< std::vector<bool> > ( n_r, std::vector<bool> ( n_t ) );
	double E_dot_B_tol = 1e-3;	// Arbitrarily chosen 20230921.
	
	for( int i=0; i<n_r; i++ ){
		for( int j=0; j<n_t; j++ ){
			
			force_free_condition_1[i][j] = ( abs( E_dot_B[i][j] ) < E_dot_B_tol );
			force_free_condition_2[i][j] = ( B_squared_minus_E_squared[i][j] > 0 );
			force_free_conditions [i][j] = ( force_free_condition_1[i][j] and force_free_condition_2[i][j] );
		}
	}
	
}



void apply_force_free_conditions_within_force_free_region(){
	// Update the electric field at all points within the last closed field line of a rotating dipole, except the poles,
	// everywhere the first or second force-free condition is violated, with an arbitrary tolerance for the value of E dot B.
	// Petri 2012, Section 3.2.
	
	//--- Reset values ---
	num_E_points_changed_for_second_FF_condition = 0;
	were_FF_conditions_applied = std::vector< std::vector<int> > ( n_r, std::vector<int> ( n_t ) );
	
	//--- Perform calculation ---
	
	if( Omega[T_index] > 0 ){
		for( int i=0; i<n_r; i++ ){
			for( int j=0; j<n_t; j++ ){
				
				if( is_gridpoint_within_FF_region[i][j] and not( force_free_conditions[i][j] ) ){
						
					num_E_points_changed_for_second_FF_condition ++;
					
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
					
				else {
					were_FF_conditions_applied[i][j] = -1;
				}

			}
		}
	}
	
	// Calculate percentage of gridpoints at which E changed. Below format makes it exactly one decimal place.
	percentage_E_points_points_changed_for_second_FF_condition = round( (double)1000 * num_E_points_changed_for_second_FF_condition / n_gridpoints_within_FF_region ) / 10;
}




std::vector< std::vector<double> > divergence_from_VSH_coeffs( std::vector< std::vector< std::vector<double> > > &v_VSH_coeffs, std::vector< std::vector< std::vector<double> > > &v_VSH_coeffs_dr ){
	
	std::vector< std::vector<double> > v_div ( n_r, std::vector<double> ( n_t ) );
	
	for( int i=0; i<n_r; i++ ){
		for( int j=0; j<n_t; j++ ){
			
			//--- L = 0 ---
			v_div[i][j] += sqrt_2Lplus1_over_4pi[0] * ( v_VSH_coeffs_dr[0][0][i] + (2.0/r[i]) * v_VSH_coeffs[0][0][i] );
			
			//--- L > 0 ---
			for( int L=1; L<=ell_max; L++ ){
				v_div[i][j] += sqrt_2Lplus1_over_4pi[L] * P0[j][L] * ( v_VSH_coeffs_dr[0][L][i] + (2.0/r[i]) * v_VSH_coeffs[0][L][i] - (L*(L+1.0)/r[i]) * v_VSH_coeffs[1][L][i] );
			}
			
		}
	}
	
	return v_div;
}




std::vector< std::vector< std::vector<double> > > curl_from_VSH_coeffs( std::vector< std::vector< std::vector<double> > > &v_VSH_coeffs, std::vector< std::vector< std::vector<double> > > &v_VSH_coeffs_dr ){
	
	std::vector< std::vector< std::vector<double> > > v_curl ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );
	
	for( int i=0; i<n_r; i++ ){
		for( int j=0; j<n_t; j++ ){
			for( int L=1; L<=ell_max; L++ ){
				
				v_curl[0][i][j] += sqrt_2Lplus1_over_4pi[L] * P0[j][L] * -L*(L+1.0) * v_VSH_coeffs[2][L][i] / r[i];
				v_curl[1][i][j] += sqrt_2Lplus1_over_4pi[L] * P1[j][L] * ( -v_VSH_coeffs_dr[2][L][i] - v_VSH_coeffs[2][L][i] / r[i] );
				v_curl[2][i][j] += sqrt_2Lplus1_over_4pi[L] * P1[j][L] * ( v_VSH_coeffs_dr[1][L][i] + ( - v_VSH_coeffs[0][L][i] + v_VSH_coeffs[1][L][i] ) / r[i] );
				
			}
		}
	}
	
	return v_curl;
}




void calculate_current_density(){
	
	J = std::vector< std::vector< std::vector<double> > > ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );
	
	for( int i=0; i<n_r; i++ ){
		for( int j=0; j<n_t; j++ ){
			
			//--- E cross B ---
			std::vector<double> E_cross_B (3);
			E_cross_B[0] = E[1][i][j] * B[2][i][j] - E[2][i][j] * B[1][i][j];
			E_cross_B[1] = E[2][i][j] * B[0][i][j] - E[0][i][j] * B[2][i][j];
			E_cross_B[2] = E[0][i][j] * B[1][i][j] - E[1][i][j] * B[0][i][j];
			
			//--- B dot curl(B) and E dot curl(E) ---
			double B_dot_curl_B = 0;
			double E_dot_curl_E = 0;
			
			for( int f=0; f<3; f++ ){
				B_dot_curl_B += B[f][i][j] * B_curl[f][i][j];
				E_dot_curl_E += E[f][i][j] * E_curl[f][i][j];
			}
			
			//--- J ---
			for( int f=0; f<3; f++ ){
				J[f][i][j] = ( E_div[i][j] * E_cross_B[f] + ( B_dot_curl_B - E_dot_curl_E ) * B[f][i][j] ) / B_squared[i][j];
			}
			
		}
	}
}
	
	





void calculate_time_derivatives_of_vector_components(){
	
	//--- Reset values ---
	B_dT         = std::vector< std::vector< std::vector<double> > > ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );
	E_dT         = std::vector< std::vector< std::vector<double> > > ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );
	B_t_dt_inner = std::vector<double> ( n_t );

	
	
	//--- Perform calculation ---
									 
	for( int i=0; i<n_r; i++ ){
		for( int j=0; j<n_t; j++ ){
			
			//--- Time derivatives ---
			for( int f=0; f<3; f++ ){
				B_dT[f][i][j] = - E_curl[f][i][j];
				E_dT[f][i][j] =   B_curl[f][i][j] - J[f][i][j];
			}
			
			//--- Frictional term, applied to all components except B_r_dT ---
			if( use_outer_sponge_layer and ( r[i] >= r_sponge ) ){
				
				for( int f=0; f<3; f++ ){
					B_dT_before_sponge[f][i][j] = B_dT[f][i][j];
					E_dT_before_sponge[f][i][j] = E_dT[f][i][j];
				}
				
				
				B_dT[1][i][j] -= friction[i] * B[1][i][j];
				B_dT[2][i][j] -= friction[i] * B[2][i][j];
				E_dT[0][i][j] -= friction[i] * E[0][i][j];
				E_dT[1][i][j] -= friction[i] * E[1][i][j];
				E_dT[2][i][j] -= friction[i] * E[2][i][j];
				
				
				for( int f=0; f<3; f++ ){
					B_dT_after_sponge[f][i][j] = B_dT[f][i][j];
					E_dT_after_sponge[f][i][j] = E_dT[f][i][j];
					
					B_dT_sponge_rel_diff[f][i][j] = ( B_dT_before_sponge[f][i][j] == 0 ) ? ( 0 ) : ( 1.0 - B_dT_before_sponge[f][i][j] / B_dT_after_sponge[f][i][j] );
					E_dT_sponge_rel_diff[f][i][j] = ( E_dT_before_sponge[f][i][j] == 0 ) ? ( 0 ) : ( 1.0 - E_dT_before_sponge[f][i][j] / E_dT_after_sponge[f][i][j] );
				}
				
			}	
		}		
	}

	//--- 20231214 Specific expressions for inner boundary ---
	
	B_t_dt_inner[0]     = ( B[1][0][1]     - B[1][0][0]            ) / delta_t;
	B_t_dt_inner.back() = ( B[1][0].back() - B[1][0][n_t-2] ) / delta_t;
        for( int j=0; j<n_t; j++ ){
	  B_t_dt_inner[j] = ( B[1][0][j+1] - B[1][0][j-1] ) / ( 2.0*delta_t );
	}

	for( int j=0; j<n_t; j++ ){
	  B_dT[0][0][j] = 0;
	  B_dT[1][0][j] = 0;
	  B_dT[2][0][j] = Omega[T_index] * ( 2.0*st[j]*B[0][0][j] + ct[j]*B[1][0][j] + st[j]*B_dr[0][0][j] + st[j]*B_t_dt_inner[j] );
	}
	
	
	//--- Set the time-derivatives of the components with boundary conditions to zero, to avoid unwanted evolution ---
	for( int j=0; j<n_t; j++ ){
	    B_dT[0][0]    [j] = 0;	// B_r at inner boundary. This is required by the BCs.
	    //B_dT[2][0]    [j] = 0;
	    //E_dT[1][0]    [j] = 0;
		E_dT[2][0]    [j] = 0;	// E_p at inner boundary. This is required by the BCs.
		//B_dT[1].back()[j] = 0;
		//B_dT[2].back()[j] = 0;
		//E_dT[1].back()[j] = 0;
		//E_dT[2].back()[j] = 0;
	}
	//B_dT[2][0][ (int) n_t / 2 - 1 ] = 0;	// Set B_phi ( R_NS, pi/2 ) = 0 to avoid evolution of noise. Suggestion from Sam 20231208.
}




void integrate_vectors_wrt_time_Adams_Bashforth_Order_1(){
	// For the Adams-Bashforth method. When we increase the timestep, all of the n=0 values should become n=1 values, and all of the n=1 values should become n=2 values.
	// This is the first-order AB method, equivalent to the Euler method. It is used at the first timestep only, so we don't need to consider two steps ago.
	
	//--- Reset values ---
	B_minus1    = std::vector< std::vector< std::vector<double> > > ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );
	E_minus1    = std::vector< std::vector< std::vector<double> > > ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );
	B_dT_minus1 = std::vector< std::vector< std::vector<double> > > ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );
	E_dT_minus1 = std::vector< std::vector< std::vector<double> > > ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );
	
	//--- Perform calculation ---
	for( int i=0; i<n_r; i++ ){
		for( int j=0; j<n_t; j++ ){
			
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
	
	//--- Reset values ---
	B_minus2    = std::vector< std::vector< std::vector<double> > > ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );
	E_minus2    = std::vector< std::vector< std::vector<double> > > ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );
	B_dT_minus2 = std::vector< std::vector< std::vector<double> > > ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );
	E_dT_minus2 = std::vector< std::vector< std::vector<double> > > ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );
	
	//--- Perform calculation ---
	for( int i=0; i<n_r; i++ ){
		for( int j=0; j<n_t; j++ ){
			
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
	
	//--- Reset values ---
	B_minus3    = std::vector< std::vector< std::vector<double> > > ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );
	E_minus3    = std::vector< std::vector< std::vector<double> > > ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );
	B_dT_minus3 = std::vector< std::vector< std::vector<double> > > ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );
	E_dT_minus3 = std::vector< std::vector< std::vector<double> > > ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );
	
	//--- Perform calculation ---
	for( int i=0; i<n_r; i++ ){
		for( int j=0; j<n_t; j++ ){
			
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




void calculate_electric_and_magnetic_energy(){
	
	total_electric_energy = pi * integral_trapezium_2D( E_squared, finite_area_element_for_r_t );	// 2pi from phi-integral, and 1/2 from ( 1/2 E_squared ), makes pi.
	total_magnetic_energy = pi * integral_trapezium_2D( B_squared, finite_area_element_for_r_t );
	
}

void calculate_normalised_B_div(){
	
	normalised_B_div = std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) );
	
	for( int i=0; i<n_r; i++ ){
		for( int j=0; j<n_t; j++ ){
			normalised_B_div[i][j] = r[i] * B_div[i][j] / sqrt(B_squared[i][j]);
		}
	}
	
	normalised_B_div_volume_integral = 2.0 * pi * integral_trapezium_2D( normalised_B_div, finite_area_element_for_r_t ) / volume_of_domain;
	
}




double estimate_csv_profiles_filesize(){
	
	double fudge_factor = 1.0;		// Free parameter to tune final result.
	int n_columns = 19 + 8 * 6;		// You are responsible for keeping this value updated based on how many values are specified in void output_headers_to_csv_profiles().
	int n_rows = 0;					// Calculated by the for-loop below.
	
	for( int T_i=0; T_i<=T_index_max; T_i++ ){
		bool condition_T_index_max_1 = ( T_i >= csv_profiles_write_T_min ) and ( T_i <= csv_profiles_write_T_max );
		bool condition_T_index_max_2 = ( T_i % csv_profiles_write_freq_T == 0 ) or ( T_i == T_index_max );
		if( condition_T_index_max_1 and condition_T_index_max_2 ){
			for( int i=csv_profiles_write_i_min; i<n_r; i++ ){
				for( int j=0; j<n_t; j++ ){
					bool condition_r = ( i <= csv_profiles_write_i_max ) and ( ( i % csv_profiles_write_freq_r == 0 ) or ( i == n_r - 1 ) );
					bool condition_t = ( j % csv_profiles_write_freq_t == 0 ) or ( j == n_t - 1 );
					if( condition_r and condition_t ){
						n_rows += 1;
					}
				}
			}
		}
	}
	
	return n_rows * n_columns * pow( 2, -20.0 ) * precision_csv * fudge_factor;	// 2^-20 goes from bytes to megabytes.
}
		




void output_headers_to_csv_profiles(){
	// Split over separate lines for readability and ease of checking consistency with output_to_csv_profiles().
	
	// Remember to update int n_columns in void estimate_csv_filesize().
	
	output_file_profiles << "T_index,i,j,x,z";
	output_file_profiles << ",B,E,J,E_dot_B,B_squared_minus_E_squared,B_r_dr2,B_div,E_div,normalised_B_div";
	output_file_profiles << ",force_free_condition_1,force_free_condition_2,force_free_conditions,is_gridpoint_within_FF_region,were_FF_conditions_applied";
	output_file_profiles << ",B_r,E_r,B_t,E_t,B_p,E_p,B_x,E_x,B_y,E_y,B_z,E_z";
	output_file_profiles << ",J_r,J_x,J_t,J_y,J_p,J_z";
	output_file_profiles << ",B_r_dr,E_r_dr,B_t_dr,E_t_dr,B_p_dr,E_p_dr";
	output_file_profiles << ",abs_rel_dev_B_r,abs_rel_dev_E_r,abs_rel_dev_B_t,abs_rel_dev_E_t,abs_rel_dev_B_p,abs_rel_dev_E_p";
	output_file_profiles << ",curl_B_r,curl_E_r,curl_B_t,curl_E_t,curl_B_p,curl_E_p";
	output_file_profiles << ",B_r_dT,E_r_dT,B_t_dT,E_t_dT,B_p_dT,E_p_dT";
	output_file_profiles << ",B_r_dT_over_B_r,E_r_dT_over_E_r,B_t_dT_over_B_t,E_t_dT_over_E_t,B_p_dT_over_B_p,E_p_dT_over_E_p";
	output_file_profiles << ",B_t_dt_inner";
	//output_file_profiles << ",B_r_dT_before_sponge,E_r_dT_before_sponge,B_t_dT_before_sponge,E_t_dT_before_sponge,B_p_dT_before_sponge,E_p_dT_before_sponge";
	//output_file_profiles << ",B_r_dT_after_sponge,E_r_dT_after_sponge,B_t_dT_after_sponge,E_t_dT_after_sponge,B_p_dT_after_sponge,E_p_dT_after_sponge";
	output_file_profiles << ",B_r_dT_sponge_rel_diff,E_r_dT_sponge_rel_diff,B_t_dT_sponge_rel_diff,E_t_dT_sponge_rel_diff,B_p_dT_sponge_rel_diff,E_p_dT_sponge_rel_diff";
	
	output_file_profiles << "\n";
}




void output_to_csv_profiles( bool condition_override = false ){
	
	bool condition_T_index_max_1 = ( T_index >= csv_profiles_write_T_min ) and ( T_index <= csv_profiles_write_T_max );
	bool condition_T_index_max_2 = ( T_index % csv_profiles_write_freq_T == 0 ) or ( T_index == T_index_max );
	
	if( ( condition_T_index_max_1 and condition_T_index_max_2 ) or condition_override ){
	
		for( int i=csv_profiles_write_i_min; i<n_r; i++ ){
			for( int j=0; j<n_t; j++ ){
				
				bool condition_r = ( i <= csv_profiles_write_i_max ) and ( ( i % csv_profiles_write_freq_r == 0 ) or ( i == n_r - 1 ) );
				bool condition_t = ( j % csv_profiles_write_freq_t == 0 ) or ( j == n_t - 1 );
		
				if( ( condition_r and condition_t ) or condition_override ){
					
					double tolerance_for_zero = 1e-13;
					
					std::vector<double> B_dT_over_B ( 3 );
					std::vector<double> E_dT_over_E ( 3 );
					
					for( int f=0; f<3; f++ ){
						B_dT_over_B[f] = ( abs( B[f][i][j] ) > tolerance_for_zero ) ? ( B_dT[f][i][j] / B[f][i][j] ) : 0;
						E_dT_over_E[f] = ( abs( E[f][i][j] ) > tolerance_for_zero ) ? ( E_dT[f][i][j] / E[f][i][j] ) : 0;
					}
					
					
					output_file_profiles <<std::setprecision(precision_csv)<< T_index <<","<< i <<","<< j <<","<< x[i][j] <<","<< z[i][j]
																		   <<","<< sqrt(B_squared[i][j]) <<","<< sqrt(E_squared[i][j]) <<","<< sqrt( pow(J[0][i][j],2.0) + pow(J[1][i][j],2.0) + pow(J[2][i][j],2.0) )
																		   <<","<< E_dot_B[i][j] <<","<< B_squared_minus_E_squared[i][j]
																		   <<","<< B_dr2[0][i][j]
																		   <<","<< B_div[i][j] <<","<< E_div[i][j] <<","<< normalised_B_div[i][j]
																		   <<","<< force_free_condition_1[i][j] <<","<< force_free_condition_2[i][j] <<","<< force_free_conditions[i][j] <<","<< is_gridpoint_within_FF_region[i][j] <<","<< were_FF_conditions_applied[i][j];
					
					for( int f=0; f<3; f++ ){ output_file_profiles <<","<< B                   [f][i][j] <<","<< E                   [f][i][j]; }
					for( int f=0; f<3; f++ ){ output_file_profiles <<","<< B_Cartesian         [f][i][j] <<","<< E_Cartesian         [f][i][j]; }
					for( int f=0; f<3; f++ ){ output_file_profiles <<","<< J                   [f][i][j] <<","<< J_Cartesian         [f][i][j]; }
					for( int f=0; f<3; f++ ){ output_file_profiles <<","<< B_dr                [f][i][j] <<","<< E_dr                [f][i][j]; }
					for( int f=0; f<3; f++ ){ output_file_profiles <<","<< B_abs_rel_dev       [f][i][j] <<","<< E_abs_rel_dev       [f][i][j]; }
					for( int f=0; f<3; f++ ){ output_file_profiles <<","<< B_curl              [f][i][j] <<","<< E_curl              [f][i][j]; }
					for( int f=0; f<3; f++ ){ output_file_profiles <<","<< B_dT                [f][i][j] <<","<< E_dT                [f][i][j]; }
					for( int f=0; f<3; f++ ){ output_file_profiles <<","<< B_dT_over_B         [f]       <<","<< E_dT_over_E         [f]      ; }
					//for( int f=0; f<3; f++ ){ output_file_profiles <<","<< B_dT_before_sponge[f][i][j] <<","<< E_dT_before_sponge[f][i][j]; }
					//for( int f=0; f<3; f++ ){ output_file_profiles <<","<< B_dT_after_sponge [f][i][j] <<","<< E_dT_after_sponge [f][i][j]; }
					for( int f=0; f<3; f++ ){ output_file_profiles <<","<< B_dT_sponge_rel_diff[f][i][j] <<","<< E_dT_sponge_rel_diff[f][i][j]; }

					output_file_profiles <<","<< B_t_dt_inner[j];
					
					output_file_profiles << "\n";
					
				}
			}
		}
	}
}




void count_nans(){
	
	B_nans = std::vector<int> ( 3 );
	E_nans = std::vector<int> ( 3 );
	
	for( int i=0; i<n_r; i++ ){
		for( int j=0; j<n_t; j++ ){
			for( int f=0; f<3; f++ ){
				// Some C++ compilers require std::isnan and std::isinf, whereas some require isnan and isinf. Comment-out as appropriate.
				if( ( std::isnan( B[f][i][j] ) ) or ( std::isinf( B[f][i][j] ) ) ){
				//if( ( isnan( B[f][i][j] ) ) or ( isinf( B[f][i][j] ) ) ){
					B_nans[f] ++;
				}
				if( ( std::isnan( E[f][i][j] ) ) or ( std::isinf( E[f][i][j] ) ) ){
				//if( ( isnan( E[f][i][j] ) ) or ( isinf( E[f][i][j] ) ) ){
					E_nans[f] ++;
				}
			}
		}
	}
	
	nans_tot = B_nans[0] + B_nans[1] + B_nans[2] + E_nans[0] + E_nans[1] + E_nans[2];
}




void output_headers_to_csv_history(){
	// Split over separate lines for readability and ease of checking consistency with output_to_csv_history().
	
	output_file_history << "T_index"
	                    << ",B_r,B_t,B_p,E_r,E_t,E_p"
						<< ",N_pts_within_FF_region,E_pts_chngd,pct_E_pts_chngd"
						<< ",stdev_B_r,stdev_E_r,stdev_B_t,stdev_E_t,stdev_B_p,stdev_E_p"
						<< ",max_abs_rel_dev_B_r,max_abs_rel_dev_E_r,max_abs_rel_dev_B_t,max_abs_rel_dev_E_t,max_abs_rel_dev_B_p,max_abs_rel_dev_E_p"
						<< ",B_outer_BCs_difference_t,E_outer_BCs_difference_t,B_outer_BCs_difference_p,E_outer_BCs_difference_p"
						<< ",nans_B_r,nans_E_r,nans_B_t,nans_E_t,nans_B_p,nans_E_p"
						<< ",nans_tot,total_electric_energy,total_magnetic_energy,normalised_B_div_volume_integral"
						<< "\n";
	
}




void output_to_csv_history(){
	
	if( ( T_index % csv_history_write_freq_T == 0 ) or ( T_index == T_index_max ) ){
		
		output_file_history <<std::setprecision(precision_csv)<< T_index
		                                                      <<","<< B[0][cout_i][cout_j] <<","<< B[1][cout_i][cout_j] <<","<< B[2][cout_i][cout_j] <<","<< E[0][cout_i][cout_j] <<","<< E[1][cout_i][cout_j] <<","<< E[2][cout_i][cout_j]
															  <<","<< n_gridpoints_within_FF_region <<","<< num_E_points_changed_for_second_FF_condition <<","<< percentage_E_points_points_changed_for_second_FF_condition;
		
		for( int f=0; f<3; f++ ){ output_file_history <<","<< B_stdev[f] <<","<< E_stdev[f]; }
		for( int f=0; f<3; f++ ){ output_file_history <<","<< B_max_abs_rel_dev[f] <<","<< E_max_abs_rel_dev[f]; }
		for( int f=1; f<3; f++ ){ output_file_history <<","<< B_outer_BCs_difference[f][cout_j] <<","<< B_outer_BCs_difference[f][cout_j]; }
		for( int f=0; f<3; f++ ){ output_file_history <<","<< B_nans[f] <<","<< E_nans[f]; }
		
		output_file_history <<","<< nans_tot <<","<< total_electric_energy <<","<< total_magnetic_energy <<","<< normalised_B_div_volume_integral << "\n";
							
	}
}




void output_headers_to_csv_VSH_coeffs(){
	// Split over separate lines for readability and ease of checking consistency with output_to_csv_history().
	
	output_file_VSH_coeffs << "T_index,i,ell,B_r_L,E_r_L,B_1_L,E_1_L,B_2_L,E_2_L,B_r_L_dr,E_r_L_dr,B_1_L_dr,E_1_L_dr,B_2_L_dr,E_2_L_dr,B_r_L_dr2\n";
	
}




void output_to_csv_VSH_coeffs(){
	
	if( ( T_index % csv_VSH_coeffs_write_freq_T == 0 ) or ( T_index == T_index_max ) ){
		for( int i=0; i<n_r; i++ ){
			for( int ell=0; ell<=ell_max; ell++ ){
				
				output_file_VSH_coeffs <<std::setprecision(precision_csv)<< T_index <<","<< i <<","<< ell;
				
				for( int f=0; f<3; f++ ){ output_file_VSH_coeffs <<","<< B_VSH_coeffs   [f][ell][i] <<","<< E_VSH_coeffs   [f][ell][i]; }
				for( int f=0; f<3; f++ ){ output_file_VSH_coeffs <<","<< B_VSH_coeffs_dr[f][ell][i] <<","<< E_VSH_coeffs_dr[f][ell][i]; }
				
				output_file_VSH_coeffs << ","<< B_VSH_coeffs_dr2[0][ell][i] << "\n";
				
			}
		}
	}
	
}




void output_headers_to_csv_BCs(){
	// Split over separate lines for readability and ease of checking consistency with output_to_csv_profiles().
	
	// Remember to update int n_columns in void estimate_csv_filesize().
	
	output_file_BCs << "T_index,j"
	
	                << ",B_inner_before_BCs_r,B_inner_after_BCs_r,B_inner_BCs_difference_r"
					<< ",B_inner_before_BCs_t,B_inner_after_BCs_t,B_inner_BCs_difference_t"
					<< ",B_inner_before_BCs_p,B_inner_after_BCs_p,B_inner_BCs_difference_p"
					
					<< ",E_inner_before_BCs_r,E_inner_after_BCs_r,E_inner_BCs_difference_r"
					<< ",E_inner_before_BCs_t,E_inner_after_BCs_t,E_inner_BCs_difference_t"
					<< ",E_inner_before_BCs_p,E_inner_after_BCs_p,E_inner_BCs_difference_p"
					
					<< ",B_outer_before_BCs_r,B_outer_after_BCs_r,B_outer_BCs_difference_r"
					<< ",B_outer_before_BCs_t,B_outer_after_BCs_t,B_outer_BCs_difference_t"
					<< ",B_outer_before_BCs_p,B_outer_after_BCs_p,B_outer_BCs_difference_p"
					
					<< ",E_outer_before_BCs_r,E_outer_after_BCs_r,E_outer_BCs_difference_r"
					<< ",E_outer_before_BCs_t,E_outer_after_BCs_t,E_outer_BCs_difference_t"
					<< ",E_outer_before_BCs_p,E_outer_after_BCs_p,E_outer_BCs_difference_p"
					
					<< "\n";
	
}




void output_to_csv_BCs(){
	
	if( ( T_index % csv_BCs_write_freq_T == 0 ) or ( T_index == T_index_max ) ){
		
		for( int j=0; j<n_t; j++ ){
		
			output_file_BCs <<std::setprecision(precision_csv)<< T_index <<","<< j;
			
			for( int f=0; f<3; f++ ){ output_file_BCs <<","<< B_inner_before_BCs[f][j] <<","<< B_inner_after_BCs [f][j] <<","<< B_inner_BCs_difference[f][j]; }
			for( int f=0; f<3; f++ ){ output_file_BCs <<","<< E_inner_before_BCs[f][j] <<","<< E_inner_after_BCs [f][j] <<","<< E_inner_BCs_difference[f][j]; }
			for( int f=0; f<3; f++ ){ output_file_BCs <<","<< B_outer_before_BCs[f][j] <<","<< B_outer_after_BCs [f][j] <<","<< B_outer_BCs_difference[f][j]; }
			for( int f=0; f<3; f++ ){ output_file_BCs <<","<< E_outer_before_BCs[f][j] <<","<< E_outer_after_BCs [f][j] <<","<< E_outer_BCs_difference[f][j]; }
			
			output_file_BCs << "\n";
		}
	}
	
}




void output_to_csv_gridpoints(){
	
	output_file_gridpoints << "i_or_j,r,t\n";
	
	for( int n=0; n<std::max(n_r,n_t); n++ ){
		
		output_file_gridpoints << n <<","<<std::setprecision(precision_csv);
		
		if( n < n_r ){
			output_file_gridpoints << r[n];
		}
		
		output_file_gridpoints <<",";
		
		if( n < n_t ){
			output_file_gridpoints << t[n];
		}
		
		output_file_gridpoints << "\n";
		
	}
	
}

void output_to_csv_time_values(){
	
	output_file_time_values << "T_index,T,T_SI,Omega,Omega_SI,P,P_SI,R_LC,R_LC_SI,star_rotation_angle\n";
	
	for( int Ti=0; Ti<T_index_max; Ti++ ){
		output_file_time_values << Ti <<","<<std::setprecision(precision_csv)<< T[Ti] <<","<< T_SI[Ti] <<","<< Omega[Ti] <<","<< Omega_SI[Ti] <<","<< P[Ti] <<","<< P_SI[Ti] <<","<< R_LC[Ti] <<","<< R_LC_SI[Ti] <<","<< star_rotation_angle[Ti] <<"\n";
	}
	
}
	




void output_headers_to_screen_and_log_file(){
	std::cout       <<std::left<<std::setw(w)<< "T_index" <<std::left<<std::setw(w)<< "T" <<"|\t"
	                <<std::left<<std::setw(w)<< "B_r" <<std::left<<std::setw(w)<< "B_t" <<std::left<<std::setw(w)<< "B_p" <<"|\t"
					<<std::left<<std::setw(w)<< "E_r" <<std::left<<std::setw(w)<< "E_t" <<std::left<<std::setw(w)<< "E_p" <<"|\t"
					<<std::left<<std::setw(w)<< "% E pts chngd" <<"|\t"
					<<std::left<<std::setw(w)<< "nans_tot" <<std::left<<std::setw(w)<< "% done" <<std::left<<std::setw(w)<< "Est. time left"
					<< std::endl;
	
	output_file_log <<std::left<<std::setw(w)<< "T_index" <<std::left<<std::setw(w)<< "T" <<"|\t"
	                <<std::left<<std::setw(w)<< "B_r" <<std::left<<std::setw(w)<< "B_t" <<std::left<<std::setw(w)<< "B_p" <<"|\t"
					<<std::left<<std::setw(w)<< "E_r" <<std::left<<std::setw(w)<< "E_t" <<std::left<<std::setw(w)<< "E_p" <<"|\t"
					<<std::left<<std::setw(w)<< "% E pts chngd" <<"|\t"
					<<std::left<<std::setw(w)<< "nans_tot" <<std::left<<std::setw(w)<< "% done" <<std::left<<std::setw(w)<< "Est. time left"
					<< "\n";
}




void output_to_screen_and_log_file(){
	
	if( ( T_index % cout_freq_T == 0 ) or ( T_index == T_index_max ) ){
	
		double percent_done = round( (double)1000*(1+T_index)/T_index_max )/10; // This format makes it exactly one decimal place. T_index starts at zero so add 1.
		
		std::cout       <<std::setprecision(precision_cout)<<std::left<<std::setw(w)<< T_index <<std::left<<std::setw(w)<< T[T_index] <<"|\t"
		                <<std::left<<std::setw(w)<< B[0][cout_i][cout_j] <<std::left<<std::setw(w)<< B[1][cout_i][cout_j] <<std::left<<std::setw(w)<< B[2][cout_i][cout_j] <<"|\t"
						<<std::left<<std::setw(w)<< E[0][cout_i][cout_j] <<std::left<<std::setw(w)<< E[1][cout_i][cout_j] <<std::left<<std::setw(w)<< E[2][cout_i][cout_j] <<"|\t"
						<<std::left<<std::setw(w)<< percentage_E_points_points_changed_for_second_FF_condition <<"|\t"
						<<std::left<<std::setw(w)<< nans_tot <<std::left<<std::setw(w)<< percent_done;
				  
		output_file_log <<std::setprecision(precision_cout)<<std::left<<std::setw(w)<< T_index <<std::left<<std::setw(w)<< T[T_index] <<"|\t"
		                <<std::left<<std::setw(w)<< B[0][cout_i][cout_j] <<std::left<<std::setw(w)<< B[1][cout_i][cout_j] <<std::left<<std::setw(w)<< B[2][cout_i][cout_j] <<"|\t"
						<<std::left<<std::setw(w)<< E[0][cout_i][cout_j] <<std::left<<std::setw(w)<< E[1][cout_i][cout_j] <<std::left<<std::setw(w)<< E[2][cout_i][cout_j] <<"|\t"
						<<std::left<<std::setw(w)<< percentage_E_points_points_changed_for_second_FF_condition <<"|\t"
						<<std::left<<std::setw(w)<< nans_tot <<std::left<<std::setw(w)<< percent_done;
		
		if( T_index == 0 ){
			std::cout << std::endl;
		} else {
			double time_now_seconds = std::chrono::high_resolution_clock::now().time_since_epoch().count() * 1e-9;
			double runtime = time_now_seconds - time_start_seconds;
			double time_left = ( -1.0 + (double) ( T_index_max - 1.0 ) / T_index ) * runtime;
			int time_left_h = time_left / 3600;
			int time_left_m = time_left / 60 - time_left_h * 60;
			int time_left_s = time_left - time_left_h * 3600 - time_left_m * 60;
			
			if( time_left_m < 10 ){
				std::cout       << time_left_h << ":0" << time_left_m;
				output_file_log << time_left_h << ":0" << time_left_m;
			} else {
				std::cout       << time_left_h << ":"  << time_left_m;
				output_file_log << time_left_h << ":"  << time_left_m;
			}
			if( time_left_s < 10 ){
				std::cout       << ":0" << time_left_s << std::endl;
				output_file_log << ":0" << time_left_s << "\n";
			} else {
				std::cout       << ":"  << time_left_s << std::endl;
				output_file_log << ":"  << time_left_s << "\n";
			}
		}
	}
}




void output_parameters( std::ostream &os ){
	// No need to estimate the size of the history CSV since it's far smaller than the profiles CSV.
	// For the future, perhaps split into multiple profiles CSVs. Then, add them all together to get an estimated total space required. A nice feature would be to fail if that much space isn't available, or if it exceeds some specified threshold.
	
	os << log_file_comment << "\n\n";
	
	os << "Time start:\t" <<std::setprecision(precision_csv)<< time_start_string << "\n";
	
	os << "Parameter values (given in code units, then in SI units)\n";
	os << "P_final                :\t" <<std::left<<std::setw(w0)<< P_final     <<std::left<<std::setw(w0)<< P_SI_final     << "\n";
	os << "Omega_final            :\t" <<std::left<<std::setw(w0)<< Omega.back() <<std::left<<std::setw(w0)<< Omega_SI.back() << "\n";
	os << "R_LC_max               :\t" <<std::left<<std::setw(w0)<< R_LC .back() <<std::left<<std::setw(w0)<< R_LC_SI .back() << "\n";
	os << "r_min, r_max           :\t" <<std::left<<std::setw(w0)<< r_min   <<std::left<<std::setw(w0)<< r_max << "\n";
	os << "n_r, delta_r           :\t" <<std::left<<std::setw(w0)<< n_r <<std::left<<std::setw(w0)<< delta_r << std::endl;
	os << "n_t, delta_t           :\t" <<std::left<<std::setw(w0)<< n_t <<std::left<<std::setw(w0)<< delta_t << std::endl;
	os << "Timestep delta_T       :\t" <<std::left<<std::setw(w0)<< delta_T <<std::left<<std::setw(w0)<< delta_T * tau << "\n";
	os << "CFL max timestep       :\t" <<std::left<<std::setw(w0)<< delta_T_CFL <<std::left<<std::setw(w0)<< delta_T_CFL * tau << "\n";
	os << "Length of simulation   :\t" <<std::left<<std::setw(w0)<< T_max << T_max * tau <<"\tsteps: " << T_index_max << "\n";
	os << "Est. CSV size (MB)     :\t" << estimate_csv_profiles_filesize() << "\n";
	os << "ell_max                :\t" << ell_max << "\n";
	os << "Printing values at ( r[" << cout_i <<"], theta[" << cout_j <<"] ) = ( " << r[cout_i] <<", " << t[cout_j] << " ) = ( " << r[cout_i] <<", " << t[cout_j]/pi << " pi )." << "\n";
	os << "dOmega_by_dT_index     :\t" << dOmega_by_dT_index <<"\tdOmega_by_dT:\t"<< dOmega_by_dT << "\n";
	os << "Ramp start             :\t" << T_index_rotation_ramp_start <<"\t"<< T_index_rotation_ramp_start * delta_T << "\n";
	os << "Ramp stop              :\t" << T_index_rotation_ramp_stop  <<"\t"<< T_index_rotation_ramp_stop  * delta_T << "\n";
	
	os << "\nuse_outer_sponge_layer:\t" << use_outer_sponge_layer << "\n";
	os << "sigma_0, gamma, beta    :\t" << friction_sigma_0 <<"\t"<< friction_gamma <<"\t"<< friction_beta << "\n";
	
	if( output_csv_profiles ){
		os << "\nCSV file saved:\t" << output_filename_profiles   << "\n";
	} else {
		os << "CSV _1_profiles not saved.\n";
	}
		
	if( output_csv_history ){
		os <<   "CSV file saved:\t" << output_filename_history    << "\n";
	} else {
		os << "CSV _2_history not saved.\n";
	}
	
	if( output_csv_VSH_coeffs ){
		os <<   "CSV file saved:\t" << output_filename_VSH_coeffs << "\n";
	} else {
		os << "CSV _3_VSH_coeffs not saved.\n";
	}
	
	if( output_csv_BCs ){
		os <<   "CSV file saved:\t" << output_filename_BCs        << "\n";
	} else {
		os << "CSV _4_BCs not saved.\n";
	}
	
	os << "Log file saved:\t" << output_filename_log << "\n";
	
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
