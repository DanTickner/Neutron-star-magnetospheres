/*
Header_Time_Evolution_43.h

Code snippets for the Time_Evolution_xx.cpp codes. Keeps them shorter and allows for two separate for-loops with Euler and Adams-Bashforth integration
without copying code.

V43: Replace VSH decomposition by finite differencing in the theta direction.
     Remove all functionality for Chebyshev and VSH series.
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
const double t_min       = 0.5 * delta_t;										// Minimum theta. Use a shifted grid to avoid the poles.
const double t_max       = pi;													// Maximum theta. Don't expect this to change from pi   but keep as variable.



//--- Output files -----
std::ofstream output_file_time_values;
std::ofstream output_file_gridpoints;
std::ofstream output_file_fields;
std::ofstream output_file_profiles;
std::ofstream output_file_history;
std::ofstream output_file_BCs;
std::ofstream output_file_log;

std::string output_filename_time_values = "CSV/"   + output_filename + "_1_time_values.csv";
std::string output_filename_gridpoints  = "CSV/"   + output_filename + "_2_gridpoints.csv";
std::string output_filename_fields      = "CSV/"   + output_filename + "_3_fields.csv";
std::string output_filename_profiles    = "CSV/"   + output_filename + "_4_profiles.csv";
std::string output_filename_history     = "CSV/"   + output_filename + "_5_history.csv";
std::string output_filename_BCs         = "CSV/"   + output_filename + "_6_BCs.csv";
std::string output_filename_log         = "Logs/"  + output_filename + ".txt";


//----- Coordinates and arrays -----
std::vector<double> r          ( n_r );	// Radial coordinate r.
std::vector<double> one_over_r ( n_r );	// 1/r. Saves repeat calculations.
std::vector<double> R          ( n_r );	// Chebyshev integration variable R = cos^-1( lambda^-1( r ) ).
std::vector<double> t          ( n_t );	// Polar coordinate theta.
std::vector<double> st         ( n_t );	// sin(theta).
std::vector<double> ct         ( n_t );	// cos(theta).
std::vector<double> ct_over_st ( n_t );	// cos(theta)/sin(theta). Saves repeat calculations.

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

std::vector< std::vector< std::vector<double> > > B                 ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );
std::vector< std::vector< std::vector<double> > > B_minus1          ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );
std::vector< std::vector< std::vector<double> > > B_minus2          ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );
std::vector< std::vector< std::vector<double> > > B_minus3          ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );
std::vector< std::vector< std::vector<double> > > B_minus4          ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );
std::vector< std::vector< std::vector<double> > > B_dT              ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );
std::vector< std::vector< std::vector<double> > > B_dT_minus1       ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );
std::vector< std::vector< std::vector<double> > > B_dT_minus2       ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );
std::vector< std::vector< std::vector<double> > > B_dT_minus3       ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );
std::vector< std::vector< std::vector<double> > > B_dT_minus4       ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );
std::vector< std::vector< std::vector<double> > > B_curl            ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );
std::vector< std::vector< std::vector<double> > > B_Cartesian       ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );
std::vector< std::vector< std::vector<double> > > B_dr              ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );
std::vector< std::vector< std::vector<double> > > B_dt              ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );

std::vector< std::vector< std::vector<double> > > E                 ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );
std::vector< std::vector< std::vector<double> > > E_minus1          ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );
std::vector< std::vector< std::vector<double> > > E_minus2          ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );
std::vector< std::vector< std::vector<double> > > E_minus3          ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );
std::vector< std::vector< std::vector<double> > > E_minus4          ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );
std::vector< std::vector< std::vector<double> > > E_dT              ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );
std::vector< std::vector< std::vector<double> > > E_dT_minus1       ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );
std::vector< std::vector< std::vector<double> > > E_dT_minus2       ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );
std::vector< std::vector< std::vector<double> > > E_dT_minus3       ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );
std::vector< std::vector< std::vector<double> > > E_dT_minus4       ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );
std::vector< std::vector< std::vector<double> > > E_curl            ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );
std::vector< std::vector< std::vector<double> > > E_Cartesian       ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );
std::vector< std::vector< std::vector<double> > > E_dr              ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );
std::vector< std::vector< std::vector<double> > > E_dt              ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );

std::vector< std::vector< std::vector<double> > > J                 ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );	// Current density vector.
std::vector< std::vector< std::vector<double> > > J_Cartesian       ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );


std::vector< std::vector<double> > B_squared                     ( n_r, std::vector<double> ( n_t ) );
std::vector< std::vector<double> > E_squared                     ( n_r, std::vector<double> ( n_t ) );
std::vector< std::vector<double> > B_squared_minus_E_squared     ( n_r, std::vector<double> ( n_t ) );
std::vector< std::vector<double> > E_dot_B                       ( n_r, std::vector<double> ( n_t ) );
std::vector< std::vector<double> > B_div                         ( n_r, std::vector<double> ( n_t ) );
std::vector< std::vector<double> > E_div                         ( n_r, std::vector<double> ( n_t ) );
std::vector< std::vector<double> > normalised_B_div              ( n_r, std::vector<double> ( n_t ) );		// r * div(B) / B, for a volume integral to characterise the size of div(B) over the domain.
std::vector< std::vector<double> > friction                      ( n_r, std::vector<double> ( n_t ) );		// Sponge layer at outer boundary to absorb outgoing waves (see Parfrey, 2012, PhD thesis, S3.9).
std::vector< std::vector<int   > > is_gridpoint_within_FF_region ( n_r, std::vector<int   > ( n_t ) );
std::vector< std::vector<int   > > were_FF_conditions_applied    ( n_r, std::vector<int   > ( n_t ) );		// +1: FF conditions applied here. -1: FF conditions not applied here. 0: Outside FF region.

std::vector<double> finite_line_element_dR ( n_r );	// Finite line elements dR for integrals to calculate Chebyshev series coefficients.
std::vector<double> finite_line_element_dt ( n_t );	// Finite line elements dtheta for integral over theta.
std::vector< std::vector<double> > finite_area_element_for_r_t ( n_r, std::vector<double> ( n_t ) );	// Finite area elements r dr dtheta for 2D integral over r and theta.

std::vector<double> B_t_dt_inner ( n_t ); // 20231214 Specify d/dT B at inner boundary; requires expression for d/dt B_t.

double total_electric_energy            = 0;
double total_magnetic_energy            = 0;
double volume_of_domain                 = 0;
double normalised_B_div_volume_integral = 0;	// Measure of the strength of div(B) across the domain, to track how div(B) evolves with time.

std::vector<int> B_nans ( 6 );	// Count the number of "nan" or "inf" values of each field component at a given timestep.
std::vector<int> E_nans ( 6 );	// Count the number of "nan" or "inf" values of each field component at a given timestep.
int nans_tot = 0;

double time_start_seconds = std::chrono::high_resolution_clock::now().time_since_epoch().count() * 1e-9;
time_t time_start_integer = time( 0 );
char*  time_start_string  = ctime( &time_start_integer );

int    T_index = 0;

int    n_gridpoints_within_FF_region                              = 0;
int    num_E_points_changed_for_second_FF_condition               = 0;	// Output to CSV to keep track of how active the FF condition application is.
double percentage_E_points_points_changed_for_second_FF_condition = 0;

std::vector< std::vector<double> > B_inner_before_BCs ( 3, std::vector<double> ( n_t ) );
std::vector< std::vector<double> > B_inner_after_BCs  ( 3, std::vector<double> ( n_t ) );
std::vector< std::vector<double> > B_outer_before_BCs ( 3, std::vector<double> ( n_t ) );
std::vector< std::vector<double> > B_outer_after_BCs  ( 3, std::vector<double> ( n_t ) );
std::vector< std::vector<double> > E_inner_before_BCs ( 3, std::vector<double> ( n_t ) );
std::vector< std::vector<double> > E_inner_after_BCs  ( 3, std::vector<double> ( n_t ) );
std::vector< std::vector<double> > E_outer_before_BCs ( 3, std::vector<double> ( n_t ) );
std::vector< std::vector<double> > E_outer_after_BCs  ( 3, std::vector<double> ( n_t ) );




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
		t [j] = t_min + j * delta_t;
	}
	
}




void calculate_functions_of_gridpoints(){
	// Keep separate to void calculate_gridpoints() so that this can be called when reloading previous runs.
	
	//--- 1/r for divergence and curl ---
	for( int i=0; i<n_r; i++ ){
		one_over_r[i] = pow( r[i], -1.0 );
	}
	
	//--- sin(theta) and cos(theta) ---
	for( int j=0; j<n_t; j++ ){
		st[j] = sin( t[j] );
		ct[j] = cos( t[j] );
		//ct_over_st[j] = ( ( j == 0 ) or ( j == n_t - 1 ) ) ? 0.0 : ct[j] / st[j];
		ct_over_st[j] = ct[j] / st[j];
	}
	
	//ct_over_st[0] = ct_over_st[1];
	//ct_over_st[n_t-1] = ct_over_st[n_t-2];
	
	//--- Cartesian coordinates ---
	for( int i=0; i<n_r; i++ ){
		for( int j=0; j<n_t; j++ ){
			x[i][j] = r[i] * st[j];
			z[i][j] = r[i] * ct[j];
		}
	}
	
	//--- Finite line and area elements for numerical integration ---
	for( int i=0; i<=n_r-2; i++ ){
		finite_line_element_dR[i] = R[i+1] - R[i];
	}
	
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




void calculate_CFL_max_timestep(){
	// This is the approximation for constant grid spacing in spherical coordinates described in the writeup.
	// Replace by an exact expression if needed, possibly using functions from Test codes/Test_CFL_condition.cpp.
	delta_T_CFL = std::min( r_min * delta_t, delta_r );
}




void calculate_friction_for_sponge_layer(){
	// Parfrey, 2012, PhD thesis, \S3.9.2.
	
	if( use_outer_sponge_layer ){
		for( int i=0; i<n_r; i++ ){
			for( int j=0; j<n_t; j++ ){
				friction[i][j] = pow( ( r[i] - r_sponge ) / ( r_max - r_sponge ), friction_beta );
				friction[i][j] = friction_sigma_0 * ( 1.0 - exp( - friction_gamma * friction[i][j] ) );
			}
		}
	}
	
}




double integral_trapezium_1D( std::vector<double> f, std::vector<double> delta_u ){
	// 1D integral of f(u) in curvilinear coordinates, given pre-calculated list of function values.
	
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




std::vector< std::vector< std::string > > data_from_CSV ( std::string input_filename ){
	// Read data from a CSV and add it to a 2D array, which is returned. Keep data as "string" format, to be converted after calling this function.
	
	std::vector< std::vector< std::string > > data;
	std::vector<std::string>                  row;
	std::string                               line, word;
	
	std::fstream input_file ( input_filename, std::ios::in );
	
	std::cout << "Reading data from:\t" << input_filename << std::endl;
	
	if( input_file.is_open() ){
		
		while( std::getline( input_file, line ) ){
			
			row.clear();
			
			std::stringstream str( line );
			
			while( std::getline( str, word, ',' ) ){	// MUST use single quotes ',' and not double quotes ","
				row.push_back( word );
			}
			
			data.push_back( row  );
		}
	}
	
	else {
		std::cout << "Could not open the file " << input_filename << std::endl;
	}
	
	return data;
}




void set_time_and_rotation_values_from_CSV(){
	
	std::vector< std::vector< std::string > > data = data_from_CSV( "CSV/" + input_file_previous_run_base + "_1_time_values.csv" );
	
	T                  .clear();
	T_SI               .clear();
	Omega              .clear();
	Omega_SI           .clear();
	P                  .clear();
	P_SI               .clear();
	R_LC               .clear();
	R_LC_SI            .clear();
	star_rotation_angle.clear();
	
	for( int n=1; n<data.size(); n++ ){
		T                  .push_back( std::stod( data[n][1] ) );
		T_SI               .push_back( std::stod( data[n][2] ) );
		Omega              .push_back( std::stod( data[n][3] ) );
		Omega_SI           .push_back( std::stod( data[n][4] ) );
		P                  .push_back( std::stod( data[n][5] ) );
		P_SI               .push_back( std::stod( data[n][6] ) );
		R_LC               .push_back( std::stod( data[n][7] ) );
		R_LC_SI            .push_back( std::stod( data[n][8] ) );
		star_rotation_angle.push_back( std::stod( data[n][9] ) );
	}
	
	delta_T = T[1] - T[0];
	
}





void set_gridpoints_from_CSV(){
	
	std::vector< std::vector< std::string > > data = data_from_CSV( "CSV/" + input_file_previous_run_base + "_2_gridpoints.csv" );
	
	r.clear();
	t.clear();
	
	//--- Build lists of coordinates ---
	for( int n=1; n<data.size(); n++ ){
		
		if( data[n][1] != "" ){
			r.push_back( std::stod( data[n][1] ) );
		}
		
		if( data[n].size() == 3 ){
			t.push_back( std::stod( data[n][2] ) );
		}
		
	}
	
	//--- Recalculate n_r and n_t, and hence reset the sizes of the vectors that use them (only B, E, B_dT, E_dT; everything else is reset before it is recaclulated) ---
	n_r     = r.size();
	n_t     = t.size();
	delta_r = r[1] - r[0];
	delta_t = t[1] - t[0];
	r_min   = r[0];
	r_max   = r.back();
	
	B    = std::vector< std::vector< std::vector<double> > > ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );
	E    = std::vector< std::vector< std::vector<double> > > ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );
	B_dT = std::vector< std::vector< std::vector<double> > > ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );
	E_dT = std::vector< std::vector< std::vector<double> > > ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );
	
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

void apply_initial_field_values_from_CSV(){
	
	std::vector< std::vector< std::string > > data = data_from_CSV( "CSV/" + input_file_previous_run_base + "_3_fields.csv" );
	
	for( int n=1; n<data.size(); n++ ){	// Avoid index i.
		if( std::stoi( data[n][0] ) == T_index_initial ){
			
			int i = std::stoi( data[n][1] );
			int j = std::stoi( data[n][2] );
			
			for( int f=0; f<3; f++ ){
				B          [f][i][j] = std::stod( data[n][ 3+f] );
				E          [f][i][j] = std::stod( data[n][ 6+f] );
				B_minus1   [f][i][j] = std::stod( data[n][ 9+f] );
				E_minus1   [f][i][j] = std::stod( data[n][12+f] );
				B_minus2   [f][i][j] = std::stod( data[n][15+f] );
				E_minus2   [f][i][j] = std::stod( data[n][18+f] );
				B_minus3   [f][i][j] = std::stod( data[n][21+f] );
				E_minus3   [f][i][j] = std::stod( data[n][24+f] );
				B_dT       [f][i][j] = std::stod( data[n][27+f] );
				E_dT       [f][i][j] = std::stod( data[n][30+f] );
				B_dT_minus1[f][i][j] = std::stod( data[n][33+f] );
				E_dT_minus1[f][i][j] = std::stod( data[n][36+f] );
				B_dT_minus2[f][i][j] = std::stod( data[n][39+f] );
				E_dT_minus2[f][i][j] = std::stod( data[n][42+f] );
				B_dT_minus3[f][i][j] = std::stod( data[n][45+f] );
				E_dT_minus3[f][i][j] = std::stod( data[n][48+f] );
			}

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
		}
	}

}




void apply_outer_boundary_conditions(){
	// Petri 2012, \S3.3.
	
	//--- Reset arrays ---
	B_inner_before_BCs = std::vector< std::vector<double> > ( 3, std::vector<double> ( n_t ) );
	B_inner_after_BCs  = std::vector< std::vector<double> > ( 3, std::vector<double> ( n_t ) );
	E_inner_before_BCs = std::vector< std::vector<double> > ( 3, std::vector<double> ( n_t ) );
	E_inner_after_BCs  = std::vector< std::vector<double> > ( 3, std::vector<double> ( n_t ) );
	
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
			B[1].back()[j] = 0.5 * ( B[1].back()[j] - E[2].back()[j] );
			B[2].back()[j] = 0.5 * ( B[2].back()[j] + E[1].back()[j] );
			E[1].back()[j] =   B[2].back()[j];
			E[2].back()[j] = - B[1].back()[j];
		}
	}
	
	//--- Record values after outer BCs applied ---
	for( int j=0; j<n_t; j++ ){
		for( int f=0; f<3; f++ ){
			B_outer_after_BCs[f][j] = B[f].back()[j];
			E_outer_after_BCs[f][j] = E[f].back()[j];
		}
	}
	
}


/*
std::vector< std::vector<double> > f_dx_FD( std::vector< std::vector<double> > &f, double h_x ){
	// Partial derivative wrt x of a function f(x,y).
	
	int N_x = f   .size();
	int N_y = f[0].size();
	std::vector< std::vector<double> > f_dx ( N_x, std::vector<double> ( N_y ) );
	
	for( int j=0; j<N_y; j++ ){
		
		for( int i=0; i<N_x-1; i++ ){
			f_dx[i][j] = ( f[i+1][j] - f[i][j] ) / h_x;
		}
		
		f_dx[N_x-1][j] = ( f[N_x-1][j] - f[N_x-2][j] ) / h_x;
	}
	
	return f_dx;
	
}

std::vector< std::vector<double> > f_dy_FD( std::vector< std::vector<double> > &f, double h_y ){
	// Partial derivative wrt y of a function f(x,y).
	
	int N_x = f   .size();
	int N_y = f[0].size();
	std::vector< std::vector<double> > f_dy ( N_x, std::vector<double> ( N_y ) );
	
	for( int i=0; i<N_x; i++ ){
		
		for( int j=0; j<N_y-1; j++ ){
			f_dy[i][j] = ( f[i][j+1] - f[i][j] ) / h_y;
		}
		
		f_dy[i][N_y-1] = ( f[i][N_y-1] - f[i][N_y-2] ) / h_y;
	}
	
	return f_dy;
	
}
*/

/*
std::vector< std::vector<double> > f_dx_FD( std::vector< std::vector<double> > &f, double h_x ){
	
	int N_x = f   .size();
	int N_y = f[0].size();
	std::vector< std::vector<double> > f_dx ( N_x, std::vector<double> ( N_y ) );
	
	for( int j=0; j<N_y; j++ ){
	
		f_dx[0    ][j] = ( -105.0*f[8    ][j] + 960.0*f[7    ][j] - 3920.0*f[6    ][j] + 9408.0*f[5    ][j] - 14700.0*f[4    ][j] + 15680.0*f[3    ][j] - 11760.0*f[2    ][j] + 6720.0*f[1    ][j] - 2283.0*f[0    ][j] ) / ( 840.0*h_x );
		f_dx[N_x-1][j] = (  105.0*f[N_x-9][j] - 960.0*f[N_x-8][j] + 3920.0*f[N_x-7][j] - 9408.0*f[N_x-6][j] + 14700.0*f[N_x-5][j] - 15680.0*f[N_x-4][j] + 11760.0*f[N_x-3][j] - 6720.0*f[N_x-2][j] + 2283.0*f[N_x-1][j] ) / ( 840.0*h_x );
		
		f_dx[1    ][j] = (   15.0*f[8    ][j] - 140.0*f[7    ][j] +  588.0*f[6    ][j] - 1470.0*f[5    ][j] +  2450.0*f[4    ][j] -  2940.0*f[3    ][j] +  2940.0*f[2    ][j] - 1338.0*f[1    ][j] -  105.0*f[0    ][j] ) / ( 840.0*h_x );
		f_dx[N_x-2][j] = (  -15.0*f[N_x-9][j] + 140.0*f[N_x-8][j] -  588.0*f[N_x-7][j] + 1470.0*f[N_x-6][j] -  2450.0*f[N_x-5][j] +  2940.0*f[N_x-4][j] -  2940.0*f[N_x-3][j] + 1338.0*f[N_x-2][j] +  105.0*f[N_x-1][j] ) / ( 840.0*h_x );
		
		f_dx[2    ][j] = (   -5.0*f[8    ][j] +  48.0*f[7    ][j] -  210.0*f[6    ][j] +  560.0*f[5    ][j] -  1050.0*f[4    ][j] +  1680.0*f[3    ][j] -   798.0*f[2    ][j] -  240.0*f[1    ][j] +   15.0*f[0    ][j] ) / ( 840.0*h_x );
		f_dx[N_x-3][j] = (    5.0*f[N_x-9][j] -  48.0*f[N_x-8][j] +  210.0*f[N_x-7][j] -  560.0*f[N_x-6][j] +  1050.0*f[N_x-5][j] -  1680.0*f[N_x-4][j] +   798.0*f[N_x-3][j] +  240.0*f[N_x-2][j] -   15.0*f[N_x-1][j] ) / ( 840.0*h_x );
		
		for( int i=3; i<=N_x-4; i++){
			f_dx[i][j] = ( f[i+3][j]-f[i-3][j] - 9.0*(f[i+2][j]-f[i-2][j]) + 45.0*(f[i+1][j]-f[i-1][j]) ) / ( 60.0*h_x );
		}
	
	}
	
	return f_dx;
	
}

std::vector< std::vector<double> > f_dy_FD( std::vector< std::vector<double> > &f, double h_y ){
	
	int N_x = f   .size();
	int N_y = f[0].size();
	std::vector< std::vector<double> > f_dy ( N_x, std::vector<double> ( N_y ) );
	
	for( int i=0; i<N_x; i++ ){
	
		f_dy[i][0    ] = ( -105.0*f[i][8    ] + 960.0*f[i][7    ] - 3920.0*f[i][6    ] + 9408.0*f[i][5    ] - 14700.0*f[i][4    ] + 15680.0*f[i][3    ] - 11760.0*f[i][2    ] + 6720.0*f[i][1    ] - 2283.0*f[i][0    ] ) / ( 840.0*h_y );
		f_dy[i][N_y-1] = (  105.0*f[i][N_y-9] - 960.0*f[i][N_y-8] + 3920.0*f[i][N_y-7] - 9408.0*f[i][N_y-6] + 14700.0*f[i][N_y-5] - 15680.0*f[i][N_y-4] + 11760.0*f[i][N_y-3] - 6720.0*f[i][N_y-2] + 2283.0*f[i][N_y-1] ) / ( 840.0*h_y );
		
		f_dy[i][1    ] = (   15.0*f[i][8    ] - 140.0*f[i][7    ] +  588.0*f[i][6    ] - 1470.0*f[i][5    ] +  2450.0*f[i][4    ] -  2940.0*f[i][3    ] +  2940.0*f[i][2    ] - 1338.0*f[i][1    ] -  105.0*f[i][0    ] ) / ( 840.0*h_y );
		f_dy[i][N_y-2] = (  -15.0*f[i][N_y-9] + 140.0*f[i][N_y-8] -  588.0*f[i][N_y-7] + 1470.0*f[i][N_y-6] -  2450.0*f[i][N_y-5] +  2940.0*f[i][N_y-4] -  2940.0*f[i][N_y-3] + 1338.0*f[i][N_y-2] +  105.0*f[i][N_y-1] ) / ( 840.0*h_y );
		
		f_dy[i][2    ] = (   -5.0*f[i][8    ] +  48.0*f[i][7    ] -  210.0*f[i][6    ] +  560.0*f[i][5    ] -  1050.0*f[i][4    ] +  1680.0*f[i][3    ] -   798.0*f[i][2    ] -  240.0*f[i][1    ] +   15.0*f[i][0    ] ) / ( 840.0*h_y );
		f_dy[i][N_y-3] = (    5.0*f[i][N_y-9] -  48.0*f[i][N_y-8] +  210.0*f[i][N_y-7] -  560.0*f[i][N_y-6] +  1050.0*f[i][N_y-5] -  1680.0*f[i][N_y-4] +   798.0*f[i][N_y-3] +  240.0*f[i][N_y-2] -   15.0*f[i][N_y-1] ) / ( 840.0*h_y );
		
		for( int j=3; j<=N_y-4; j++){
			f_dy[i][j] = ( f[i][j+3]-f[i][j-3] - 9.0*(f[i][j+2]-f[i][j-2]) + 45.0*(f[i][j+1]-f[i][j-1]) ) / ( 60.0*h_y );
		}
	
	}
	
	return f_dy;
	
}
*/


std::vector< std::vector<double> > f_dx_FD( std::vector< std::vector<double> > &f, double h ){
	
	int N_x = f   .size();
	int N_y = f[0].size();
	std::vector< std::vector<double> > f_dx ( N_x, std::vector<double> ( N_y ) );
	double one_over_h = 1.0 / h;
	
	for( int j=0; j<N_y; j++ ){
	
		f_dx[0    ][j] = -0.1*f[10    ][j] + (10.0/9.0)*f[9     ][j] - 5.625*f[8    ][j] + (120.0/7.0)*f[7    ][j] - 35.0*f[6    ][j] + 50.4*f[5    ][j] - 52.5*f[4    ][j] + 40.0*f[3    ][j] - 22.5*f[2    ][j] + 10.0*f[1    ][j] - (7381.0/2520.0)*f[0    ][j];
		f_dx[N_x-1][j] =  0.1*f[N_x-11][j] - (10.0/9.0)*f[N_x-10][j] + 5.625*f[N_x-9][j] - (120.0/7.0)*f[N_x-8][j] + 35.0*f[N_x-7][j] - 50.4*f[N_x-6][j] + 52.5*f[N_x-5][j] - 40.0*f[N_x-4][j] + 22.5*f[N_x-3][j] - 10.0*f[N_x-2][j] + (7381.0/2520.0)*f[N_x-1][j];
		
		f_dx[1    ][j] =  (1.0/90.0)*f[10    ][j] -  0.125*f[9     ][j] +  (9.0/14.0)*f[8    ][j] - 2.0*f[7    ][j] + 4.2*f[6    ][j] -  6.3*f[5    ][j] +  7.0*f[4    ][j] -  6.0*f[3    ][j] + 4.5*f[2    ][j] -  (4609.0/2520.0)*f[1    ][j] - 0.1*f[0    ][j];
		f_dx[N_x-2][j] = -(1.0/90.0)*f[N_x-11][j] +  0.125*f[N_x-10][j] -  (9.0/14.0)*f[N_x-9][j] + 2.0*f[N_x-8][j] - 4.2*f[N_x-7][j] +  6.3*f[N_x-6][j] -  7.0*f[N_x-5][j] +  6.0*f[N_x-4][j] - 4.5*f[N_x-3][j] +  (4609.0/2520.0)*f[N_x-2][j] + 0.1*f[N_x-1][j];
		
		f_dx[2    ][j] = -(1.0/360.0)*f[10    ][j] + (2.0/63.0)*f[9     ][j] - (1.0/6.0)*f[8    ][j] + (8.0/15.0)*f[7    ][j] -  (7.0/6.0)*f[6    ][j] +  (28.0/15.0)*f[5    ][j] -  (7.0/3.0)*f[4    ][j] + (8.0/3.0)*f[3    ][j] - (341.0/280.0)*f[2    ][j] - (2.0/9.0)*f[1    ][j] + (1.0/90.0)*f[0    ][j];
		f_dx[N_x-3][j] =  (1.0/360.0)*f[N_x-11][j] - (2.0/63.0)*f[N_x-10][j] + (1.0/6.0)*f[N_x-9][j] - (8.0/15.0)*f[N_x-8][j] +  (7.0/6.0)*f[N_x-7][j] -  (28.0/15.0)*f[N_x-6][j] +  (7.0/3.0)*f[N_x-5][j] - (8.0/3.0)*f[N_x-4][j] + (341.0/280.0)*f[N_x-3][j] + (2.0/9.0)*f[N_x-2][j] - (1.0/90.0)*f[N_x-1][j];
		
		f_dx[3    ][j] =  (1.0/840.0)*f[10    ][j] - (1.0/72.0)*f[9     ][j] + 0.075*f[8    ][j] - 0.25*f[7    ][j] + (7.0/12.0)*f[6    ][j] - 1.05*f[5    ][j] + 1.75*f[4    ][j] - (319.0/420.0)*f[3    ][j] - 0.375*f[2    ][j] + (1.0/24.0)*f[1    ][j] - (1.0/360.0)*f[0    ][j];
		f_dx[N_x-4][j] = -(1.0/840.0)*f[N_x-11][j] + (1.0/72.0)*f[N_x-10][j] - 0.075*f[N_x-9][j] + 0.25*f[N_x-8][j] - (7.0/12.0)*f[N_x-7][j] + 1.05*f[N_x-6][j] - 1.75*f[N_x-5][j] + (319.0/420.0)*f[N_x-4][j] + 0.375*f[N_x-3][j] - (1.0/24.0)*f[N_x-2][j] + (1.0/360.0)*f[N_x-1][j];
		
		f_dx[4    ][j] = -(1.0/1260.0)*f[10    ][j] + (1.0/105.0)*f[9     ][j] - (3.0/56.0)*f[8    ][j] + (4.0/21.0)*f[7    ][j] - 0.5*f[6    ][j] + 1.2*f[5    ][j] - (11.0/30.0)*f[4    ][j] - (4.0/7.0)*f[3    ][j] + (3.0/28.0)*f[2    ][j] - (1.0/63.0)*f[1    ][j] + (1.0/840.0)*f[0    ][j];
		f_dx[N_x-5][j] =  (1.0/1260.0)*f[N_x-11][j] - (1.0/105.0)*f[N_x-10][j] + (3.0/56.0)*f[N_x-9][j] - (4.0/21.0)*f[N_x-8][j] + 0.5*f[N_x-7][j] - 1.2*f[N_x-6][j] + (11.0/30.0)*f[N_x-5][j] + (4.0/7.0)*f[N_x-4][j] - (3.0/28.0)*f[N_x-3][j] + (1.0/63.0)*f[N_x-2][j] - (1.0/840.0)*f[N_x-1][j];
		
		for( int i=5; i<=N_x-6; i++ ){
			f_dx[i][j] = (1.0/1260.0)*(f[i+5][j]-f[i-5][j]) - (5.0/504.0)*(f[i+4][j]-f[i-4][j]) + (5.0/84.0)*(f[i+3][j]-f[i-3][j]) - (5.0/21.0)*(f[i+2][j]-f[i-2][j]) + (5.0/6.0)*(f[i+1][j]-f[i-1][j]);
		}
		
		for( int i=0; i<N_x; i++ ){
			f_dx[i][j] *= one_over_h;
		}
	
	}
	
	return f_dx;
	
}

std::vector< std::vector<double> > f_dy_FD( std::vector< std::vector<double> > &f, double h ){
	
	int N_x = f   .size();
	int N_y = f[0].size();
	std::vector< std::vector<double> > f_dy ( N_x, std::vector<double> ( N_y ) );
	double one_over_h = 1.0 / h;
	
	for( int i=0; i<N_x; i++ ){
	
		f_dy[i][0    ] = -0.1*f[i][10    ] + (10.0/9.0)*f[i][9     ] - 5.625*f[i][8    ] + (120.0/7.0)*f[i][7    ] - 35.0*f[i][6    ] + 50.4*f[i][5    ] - 52.5*f[i][4    ] + 40.0*f[i][3    ] - 22.5*f[i][2    ] + 10.0*f[i][1    ] - (7381.0/2520.0)*f[i][0    ];
		f_dy[i][N_y-1] =  0.1*f[i][N_y-11] - (10.0/9.0)*f[i][N_y-10] + 5.625*f[i][N_y-9] - (120.0/7.0)*f[i][N_y-8] + 35.0*f[i][N_y-7] - 50.4*f[i][N_y-6] + 52.5*f[i][N_y-5] - 40.0*f[i][N_y-4] + 22.5*f[i][N_y-3] - 10.0*f[i][N_y-2] + (7381.0/2520.0)*f[i][N_y-1];
		
		f_dy[i][1    ] =  (1.0/90.0)*f[i][10    ] -  0.125*f[i][9     ] +  (9.0/14.0)*f[i][8    ] - 2.0*f[i][7    ] + 4.2*f[i][6    ] -  6.3*f[i][5    ] +  7.0*f[i][4    ] -  6.0*f[i][3    ] + 4.5*f[i][2    ] -  (4609.0/2520.0)*f[i][1    ] - 0.1*f[i][0    ];
		f_dy[i][N_y-2] = -(1.0/90.0)*f[i][N_y-11] +  0.125*f[i][N_y-10] -  (9.0/14.0)*f[i][N_y-9] + 2.0*f[i][N_y-8] - 4.2*f[i][N_y-7] +  6.3*f[i][N_y-6] -  7.0*f[i][N_y-5] +  6.0*f[i][N_y-4] - 4.5*f[i][N_y-3] +  (4609.0/2520.0)*f[i][N_y-2] + 0.1*f[i][N_y-1];
		
		f_dy[i][2    ] = -(1.0/360.0)*f[i][10    ] + (2.0/63.0)*f[i][9     ] - (1.0/6.0)*f[i][8    ] + (8.0/15.0)*f[i][7    ] -  (7.0/6.0)*f[i][6    ] +  (28.0/15.0)*f[i][5    ] -  (7.0/3.0)*f[i][4    ] + (8.0/3.0)*f[i][3    ] - (341.0/280.0)*f[i][2    ] - (2.0/9.0)*f[i][1    ] + (1.0/90.0)*f[i][0    ];
		f_dy[i][N_y-3] =  (1.0/360.0)*f[i][N_y-11] - (2.0/63.0)*f[i][N_y-10] + (1.0/6.0)*f[i][N_y-9] - (8.0/15.0)*f[i][N_y-8] +  (7.0/6.0)*f[i][N_y-7] -  (28.0/15.0)*f[i][N_y-6] +  (7.0/3.0)*f[i][N_y-5] - (8.0/3.0)*f[i][N_y-4] + (341.0/280.0)*f[i][N_y-3] + (2.0/9.0)*f[i][N_y-2] - (1.0/90.0)*f[i][N_y-1];
		
		f_dy[i][3    ] =  (1.0/840.0)*f[i][10    ] - (1.0/72.0)*f[i][9     ] + 0.075*f[i][8    ] - 0.25*f[i][7    ] + (7.0/12.0)*f[i][6    ] - 1.05*f[i][5    ] + 1.75*f[i][4    ] - (319.0/420.0)*f[i][3    ] - 0.375*f[i][2    ] + (1.0/24.0)*f[i][1    ] - (1.0/360.0)*f[i][0    ];
		f_dy[i][N_y-4] = -(1.0/840.0)*f[i][N_y-11] + (1.0/72.0)*f[i][N_y-10] - 0.075*f[i][N_y-9] + 0.25*f[i][N_y-8] - (7.0/12.0)*f[i][N_y-7] + 1.05*f[i][N_y-6] - 1.75*f[i][N_y-5] + (319.0/420.0)*f[i][N_y-4] + 0.375*f[i][N_y-3] - (1.0/24.0)*f[i][N_y-2] + (1.0/360.0)*f[i][N_y-1];
		
		f_dy[i][4    ] = -(1.0/1260.0)*f[i][10    ] + (1.0/105.0)*f[i][9     ] - (3.0/56.0)*f[i][8    ] + (4.0/21.0)*f[i][7    ] - 0.5*f[i][6    ] + 1.2*f[i][5    ] - (11.0/30.0)*f[i][4    ] - (4.0/7.0)*f[i][3    ] + (3.0/28.0)*f[i][2    ] - (1.0/63.0)*f[i][1    ] + (1.0/840.0)*f[i][0    ];
		f_dy[i][N_y-5] =  (1.0/1260.0)*f[i][N_y-11] - (1.0/105.0)*f[i][N_y-10] + (3.0/56.0)*f[i][N_y-9] - (4.0/21.0)*f[i][N_y-8] + 0.5*f[i][N_y-7] - 1.2*f[i][N_y-6] + (11.0/30.0)*f[i][N_y-5] + (4.0/7.0)*f[i][N_y-4] - (3.0/28.0)*f[i][N_y-3] + (1.0/63.0)*f[i][N_y-2] - (1.0/840.0)*f[i][N_y-1];
		
		for( int j=5; j<=N_y-6; j++ ){
			f_dy[i][j] = (1.0/1260.0)*(f[i][j+5]-f[i][j-5]) - (5.0/504.0)*(f[i][j+4]-f[i][j-4]) + (5.0/84.0)*(f[i][j+3]-f[i][j-3]) - (5.0/21.0)*(f[i][j+2]-f[i][j-2]) + (5.0/6.0)*(f[i][j+1]-f[i][j-1]);
		}
		
		for( int j=0; j<N_y; j++ ){
			f_dy[i][j] *= one_over_h;
		}
	
	}
	
	return f_dy;
	
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
			B_squared                [i][j] = pow( B[0][i][j], 2 ) + pow( B[1][i][j], 2 ) + pow( B[2][i][j], 2.0 );
			E_squared                [i][j] = pow( E[0][i][j], 2 ) + pow( E[1][i][j], 2 ) + pow( E[2][i][j], 2.0 );
			B_squared_minus_E_squared[i][j] = B_squared[i][j] - E_squared[i][j];
			E_dot_B                  [i][j] = E[0][i][j]*B[0][i][j] + E[1][i][j]*B[1][i][j] + E[2][i][j]*B[2][i][j];
		}
	}
	
}




void ramp_up_electric_fields_for_rotation(){
	if( ( T_index >= T_index_rotation_ramp_start ) and ( T_index <= T_index_rotation_ramp_stop ) ){
	
		int n_steps_rampup = T_index_rotation_ramp_stop - T_index_rotation_ramp_start + 1;
		
		for( int i=0; i<n_r; i++ ){
			for( int j=1; j<n_t-1; j++ ){
				
				double factor = Omega.back() * r[i] * st[j] / (double) n_steps_rampup;
				
				E[0][i][j] += factor * B_t_function( r[i], t[j] );
				E[1][i][j] -= factor * B_r_function( r[i], t[j] );
				
			}
		}
		
		apply_inner_boundary_conditions();
		apply_outer_boundary_conditions();
	}
}




void calculate_n_gridpoints_within_FF_region(){
	
	//--- Reset values ---
	is_gridpoint_within_FF_region = std::vector< std::vector<int> > ( n_r, std::vector<int> ( n_t ) );
	n_gridpoints_within_FF_region = 0;
	
	//--- Perform calculation ---
	if( Omega[T_index] > 0 ){
		for( int i=0; i<n_r; i++ ){
			for( int j=1; j<n_t-1; j++ ){
				if( r[i] * st[j] < R_LC[T_index] ){
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
	
	//--- Reset values ---
	num_E_points_changed_for_second_FF_condition = 0;
	were_FF_conditions_applied = std::vector< std::vector<int> > ( n_r, std::vector<int> ( n_t ) );
	
	//--- Perform calculation ---
	calculate_n_gridpoints_within_FF_region();
	
	if( Omega[T_index] > 0 ){
		for( int i=0; i<n_r; i++ ){
			for( int j=1; j<n_t-1; j++ ){
				
				if( is_gridpoint_within_FF_region[i][j] ){
					
					double E_dot_B_tol = 1e-3;	// Arbitrarily chosen 20230921.
					
					if( ( E_dot_B[i][j] > E_dot_B_tol ) or ( B_squared_minus_E_squared[i][j] < 0 ) ){
						
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



std::vector< std::vector<double> > divergence_from_FD_derivatives( std::vector< std::vector< std::vector<double> > > v, std::vector< std::vector< std::vector<double> > > v_dr, std::vector< std::vector< std::vector<double> > > v_dt ){
	std::vector< std::vector<double> > v_div ( n_r, std::vector<double> ( n_t ) );
	
	for( int i=0; i<n_r; i++ ){
		for( int j=0; j<n_t; j++ ){
			v_div[i][j] = v_dr[0][i][j] + one_over_r[i] * ( 2.0*v[0][i][j] + v_dt[1][i][j] + ct_over_st[j] * v[1][i][j] );
		}
	}
	
	return v_div;
}


std::vector< std::vector< std::vector<double> > > curl_from_FD_derivatives( std::vector< std::vector< std::vector<double> > > v, std::vector< std::vector< std::vector<double> > > v_dr, std::vector< std::vector< std::vector<double> > > v_dt ){
	
	std::vector< std::vector< std::vector<double> > > v_curl ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );
	
	for( int i=0; i<n_r; i++ ){
		for( int j=0; j<n_t; j++ ){
			v_curl[0][i][j] = one_over_r[i] * ( v_dt[2][i][j] + ct_over_st[j] * v[2][i][j] );
			v_curl[1][i][j] = - v_dr[2][i][j] - one_over_r[i] * v[2][i][j];
			v_curl[2][i][j] = v_dr[1][i][j] + one_over_r[i] * ( v[1][i][j] - v_dt[0][i][j] );
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
			
			//--- J ---
			double B_dot_curl_B = B[0][i][j]*B_curl[0][i][j] + B[1][i][j]*B_curl[1][i][j] + B[2][i][j]*B_curl[2][i][j];
			double E_dot_curl_E = E[0][i][j]*E_curl[0][i][j] + E[1][i][j]*E_curl[1][i][j] + E[2][i][j]*E_curl[2][i][j];
			
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
		for( int j=1; j<n_t-1; j++ ){
			
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
	B_t_dt_inner.back() = ( B[1][0].back() - B[1][0][n_t-2] ) / delta_t;
        for( int j=1; j<n_t-1; j++ ){
	  B_t_dt_inner[j] = ( B[1][0][j+1] - B[1][0][j-1] ) / ( 2.0*delta_t );
	}

	for( int j=0; j<n_t; j++ ){
	  B_dT[0][0][j] = 0;
	  B_dT[1][0][j] = 0;
	  B_dT[2][0][j] = Omega[T_index] * ( 2.0*st[j]*B[0][0][j] + ct[j]*B[1][0][j] + st[j]*B_dr[0][0][j] + st[j]*B_t_dt_inner[j] );
	}
	
	//--- Set the time-derivatives of the components with boundary conditions to zero, to avoid unwanted evolution ---
	for( int j=0; j<n_t; j++ ){
	  //B_dT[0][0]    [j] = 0;
	  //B_dT[2][0]    [j] = 0;
		E_dT[1][0]    [j] = 0;
		E_dT[2][0]    [j] = 0;
		//B_dT[1].back()[j] = 0;
		//B_dT[2].back()[j] = 0;
		E_dT[1].back()[j] = 0;
		E_dT[2].back()[j] = 0;
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
		for( int j=1; j<n_t-1; j++ ){
			
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
		for( int j=1; j<n_t-1; j++ ){
			
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
		for( int j=1; j<n_t-1; j++ ){
			
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
				B[f][i][j] = B_minus1[f][i][j] + delta_T * ( (23.0/12.0) * B_dT_minus1[f][i][j] - (4.0/3.0) * B_dT_minus2[f][i][j] + (5.0/12.0) * B_dT_minus3[f][i][j] );
				E[f][i][j] = E_minus1[f][i][j] + delta_T * ( (23.0/12.0) * E_dT_minus1[f][i][j] - (4.0/3.0) * E_dT_minus2[f][i][j] + (5.0/12.0) * E_dT_minus3[f][i][j] );
			}
			
		}
	}
}




void integrate_vectors_wrt_time_Adams_Bashforth_Order_4(){
	// For the Adams-Bashforth method. When we increase the timestep, all of the n=0 values should become n=1 values, and all of the n=1 values should become n=2 values.
	
	//--- Reset values ---
	B_minus4    = std::vector< std::vector< std::vector<double> > > ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );
	E_minus4    = std::vector< std::vector< std::vector<double> > > ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );
	B_dT_minus4 = std::vector< std::vector< std::vector<double> > > ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );
	E_dT_minus4 = std::vector< std::vector< std::vector<double> > > ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );
	
	//--- Perform calculation ---
	for( int i=0; i<n_r; i++ ){
		for( int j=1; j<n_t-1; j++ ){
			
			//--- Setup the values four steps ago ---
			for( int f=0; f<3; f++ ){
				B_minus4   [f][i][j] = B_minus3   [f][i][j];
				E_minus4   [f][i][j] = E_minus3   [f][i][j];
				B_dT_minus4[f][i][j] = B_dT_minus3[f][i][j];
				E_dT_minus4[f][i][j] = E_dT_minus3[f][i][j];
			}
			
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
				B[f][i][j] = B_minus1[f][i][j] + delta_T * ( (55.0/24.0) * B_dT_minus1[f][i][j] - (59.0/24.0) * B_dT_minus2[f][i][j] + (37.0/24.0) * B_dT_minus3[f][i][j] - 0.375*B_dT_minus4[f][i][j] );
				E[f][i][j] = E_minus1[f][i][j] + delta_T * ( (55.0/24.0) * E_dT_minus1[f][i][j] - (59.0/24.0) * E_dT_minus2[f][i][j] + (37.0/24.0) * E_dT_minus3[f][i][j] - 0.375*E_dT_minus4[f][i][j] );
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
		bool conditioT_index_max_1 = ( T_i >= csv_profiles_write_T_min ) and ( T_i <= csv_profiles_write_T_max );
		bool conditioT_index_max_2 = ( T_i % csv_profiles_write_freq_T == 0 ) or ( T_i == T_index_max );
		if( conditioT_index_max_1 and conditioT_index_max_2 ){
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
	output_file_profiles << ",B,E,J,E_dot_B,B_squared_minus_E_squared,B_div,E_div,normalised_B_div,is_gridpoint_within_FF_region,were_FF_conditions_applied";
	output_file_profiles << ",B_x,E_x,B_y,E_y,B_z,E_z";
	output_file_profiles << ",J_r,J_x,J_t,J_y,J_p,J_z";
	output_file_profiles << ",B_r_dr,E_r_dr,B_t_dr,E_t_dr,B_p_dr,E_p_dr";
	output_file_profiles << ",curl_B_r,curl_E_r,curl_B_t,curl_E_t,curl_B_p,curl_E_p";
	output_file_profiles << ",B_r_dT,E_r_dT,B_t_dT,E_t_dT,B_p_dT,E_p_dT";
	output_file_profiles << ",B_r_dT_over_B_r,E_r_dT_over_E_r,B_t_dT_over_B_t,E_t_dT_over_E_t,B_p_dT_over_B_p,E_p_dT_over_E_p";
	output_file_profiles << ",B_t_dt_inner";
	
	output_file_profiles << "\n";
}




void output_to_csv_profiles( bool condition_override = false ){
	
	bool conditioT_index_max_1 = ( T_index >= csv_profiles_write_T_min ) and ( T_index <= csv_profiles_write_T_max );
	bool conditioT_index_max_2 = ( T_index % csv_profiles_write_freq_T == 0 ) or ( T_index == T_index_max );
	
	if( ( conditioT_index_max_1 and conditioT_index_max_2 ) or condition_override ){
	
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
																		   <<","<< B_div[i][j] <<","<< E_div[i][j] <<","<< normalised_B_div[i][j]
																		   <<","<< is_gridpoint_within_FF_region[i][j] <<","<< were_FF_conditions_applied[i][j];
					
					for( int f=0; f<3; f++ ){ output_file_profiles <<","<< B_Cartesian  [f][i][j] <<","<< E_Cartesian  [f][i][j]; }
					for( int f=0; f<3; f++ ){ output_file_profiles <<","<< J            [f][i][j] <<","<< J_Cartesian  [f][i][j]; }
					for( int f=0; f<3; f++ ){ output_file_profiles <<","<< B_dr         [f][i][j] <<","<< E_dr         [f][i][j]; }
					for( int f=0; f<3; f++ ){ output_file_profiles <<","<< B_curl       [f][i][j] <<","<< E_curl       [f][i][j]; }
					for( int f=0; f<3; f++ ){ output_file_profiles <<","<< B_dT         [f][i][j] <<","<< E_dT         [f][i][j]; }
					for( int f=0; f<3; f++ ){ output_file_profiles <<","<< B_dT_over_B  [f]       <<","<< E_dT_over_E  [f]      ; }

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
						<< ",E_pts_chngd,pct_E_pts_chngd"
						<< ",nans_B_r,nans_E_r,nans_B_t,nans_E_t,nans_B_p,nans_E_p"
						<< ",nans_tot,total_electric_energy,total_magnetic_energy,normalised_B_div_volume_integral"
						<< "\n";
	
}




void output_to_csv_history(){
	
	if( ( T_index % csv_history_write_freq_T == 0 ) or ( T_index == T_index_max ) ){
		
		output_file_history <<std::setprecision(precision_csv)<< T_index
		                                                      <<","<< B[0][cout_i][cout_j] <<","<< B[1][cout_i][cout_j] <<","<< B[2][cout_i][cout_j] <<","<< E[0][cout_i][cout_j] <<","<< E[1][cout_i][cout_j] <<","<< E[2][cout_i][cout_j]
															  <<","<< num_E_points_changed_for_second_FF_condition <<","<< percentage_E_points_points_changed_for_second_FF_condition;
		
		for( int f=0; f<3; f++ ){ output_file_history <<","<< B_nans[f]            <<","<< E_nans[f]; }
		
		output_file_history <<","<< nans_tot <<","<< total_electric_energy <<","<< total_magnetic_energy <<","<< normalised_B_div_volume_integral << "\n";
							
	}
}




void output_headers_to_csv_BCs(){
	// Split over separate lines for readability and ease of checking consistency with output_to_csv_profiles().
	
	// Remember to update int n_columns in void estimate_csv_filesize().
	
	output_file_BCs << "T_index,j";
	output_file_BCs << ",B_r_inner_before_BCs,E_r_inner_before_BCs,B_t_inner_before_BCs,E_t_inner_before_BCs,B_p_inner_before_BCs,E_p_inner_before_BCs";
	output_file_BCs << ",B_r_inner_after_BCs,E_r_inner_after_BCs,B_t_inner_after_BCs,E_t_inner_after_BCs,B_p_inner_after_BCs,E_p_inner_after_BCs";
	output_file_BCs << ",rel_diff_B_r_inner,rel_diff_E_r_inner,rel_diff_B_t_inner,rel_diff_E_t_inner,rel_diff_B_p_inner,rel_diff_E_p_inner";
	output_file_BCs << ",B_r_outer_before_BCs,E_r_outer_before_BCs,B_t_outer_before_BCs,E_t_outer_before_BCs,B_p_outer_before_BCs,E_p_outer_before_BCs";
	output_file_BCs << ",B_r_outer_after_BCs,E_r_outer_after_BCs,B_t_outer_after_BCs,E_t_outer_after_BCs,B_p_outer_after_BCs,E_p_outer_after_BCs";
	output_file_BCs << ",rel_diff_B_r_outer,rel_diff_E_r_outer,rel_diff_B_t_outer,rel_diff_E_t_outer,rel_diff_B_p_outer,rel_diff_E_p_outer";
	output_file_BCs << "\n";
	
}




void output_to_csv_BCs(){
	
	if( ( T_index % csv_BCs_write_freq_T == 0 ) or ( T_index == T_index_max ) ){
		
		for( int j=0; j<n_t; j++ ){
			
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
		
			output_file_BCs <<std::setprecision(precision_csv)<< T_index <<","<< j <<",";
			
			for( int f=0; f<3; f++ ){ output_file_BCs <<","<< B_inner_before_BCs[f][j] <<","<< E_inner_before_BCs[f][j]; }
			for( int f=0; f<3; f++ ){ output_file_BCs <<","<< B_inner_after_BCs [f][j] <<","<< E_inner_after_BCs [f][j]; }
			for( int f=0; f<3; f++ ){ output_file_BCs <<","<< B_rel_diff_inner  [f]    <<","<< E_rel_diff_inner  [f]   ; }
			for( int f=0; f<3; f++ ){ output_file_BCs <<","<< B_outer_before_BCs[f][j] <<","<< E_outer_before_BCs[f][j]; }
			for( int f=0; f<3; f++ ){ output_file_BCs <<","<< B_outer_after_BCs [f][j] <<","<< E_outer_after_BCs [f][j]; }
			for( int f=0; f<3; f++ ){ output_file_BCs <<","<< B_rel_diff_outer  [f]    <<","<< E_rel_diff_outer  [f]   ; }
			
			output_file_BCs << "\n";
		}
	}
	
}




void output_headers_to_csv_fields(){
	output_file_fields << "T_index,i,j,B_r,B_t,B_p,E_r,E_t,E_p"
	                   << ",B_r_minus1,B_t_minus1,B_p_minus1,E_r_minus1,E_t_minus1,E_p_minus1"
					   << ",B_r_minus2,B_t_minus2,B_p_minus2,E_r_minus2,E_t_minus2,E_p_minus2"
					   << ",B_r_minus3,B_t_minus3,B_p_minus3,E_r_minus3,E_t_minus3,E_p_minus3"
					   << ",B_r_dT,B_t_dT,B_p_dT,E_r_dT,E_t_dT,E_p_dT"
					   << ",B_r_dT_minus1,B_t_dT_minus1,B_p_dT_minus1,E_r_dT_minus1,E_t_dT_minus1,E_p_dT_minus1"
					   << ",B_r_dT_minus2,B_t_dT_minus2,B_p_dT_minus2,E_r_dT_minus2,E_t_dT_minus2,E_p_dT_minus2"
					   << ",B_r_dT_minus3,B_t_dT_minus3,B_p_dT_minus3,E_r_dT_minus3,E_t_dT_minus3,E_p_dT_minus3\n";
}




void output_to_csv_fields( bool condition_override = false ){
	
	bool condition = ( T_index >= csv_fields_write_T_min ) and ( ( T_index % csv_fields_write_freq_T == 0 ) or ( T_index == T_index_max ) );
	
	if( ( condition ) or ( condition_override ) ){
		for( int i=0; i<n_r; i++ ){
			for( int j=0; j<n_t; j++ ){
				output_file_fields <<std::setprecision(precision_csv)<< T_index <<","<< i <<","<< j;
				
				for( int f=0; f<3; f++ ){ output_file_fields <<","<< B          [f][i][j]; }
				for( int f=0; f<3; f++ ){ output_file_fields <<","<< E          [f][i][j]; }
				for( int f=0; f<3; f++ ){ output_file_fields <<","<< B_minus1   [f][i][j]; }
				for( int f=0; f<3; f++ ){ output_file_fields <<","<< E_minus1   [f][i][j]; }
				for( int f=0; f<3; f++ ){ output_file_fields <<","<< B_minus2   [f][i][j]; }
				for( int f=0; f<3; f++ ){ output_file_fields <<","<< E_minus2   [f][i][j]; }
				for( int f=0; f<3; f++ ){ output_file_fields <<","<< B_minus3   [f][i][j]; }
				for( int f=0; f<3; f++ ){ output_file_fields <<","<< E_minus3   [f][i][j]; }
				for( int f=0; f<3; f++ ){ output_file_fields <<","<< B_dT       [f][i][j]; }
				for( int f=0; f<3; f++ ){ output_file_fields <<","<< E_dT       [f][i][j]; }
				for( int f=0; f<3; f++ ){ output_file_fields <<","<< B_dT_minus1[f][i][j]; }
				for( int f=0; f<3; f++ ){ output_file_fields <<","<< E_dT_minus1[f][i][j]; }
				for( int f=0; f<3; f++ ){ output_file_fields <<","<< B_dT_minus2[f][i][j]; }
				for( int f=0; f<3; f++ ){ output_file_fields <<","<< E_dT_minus2[f][i][j]; }
				for( int f=0; f<3; f++ ){ output_file_fields <<","<< B_dT_minus3[f][i][j]; }
				for( int f=0; f<3; f++ ){ output_file_fields <<","<< E_dT_minus3[f][i][j]; }
				
				output_file_fields << "\n";
			}
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
