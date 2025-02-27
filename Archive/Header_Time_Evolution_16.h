/*
Header_Time_Evolution_16.h

Code snippets for the Time_Evolution_xx.cpp codes. Keeps them shorter and allows for two separate for-loops with Euler and Adams-Bashforth integration
without copying code.

V16: Combine V14 (up-to-date header file without ghost points) and V15 (initial version with ghost points).
     ABANDON GHOST POINTS. REVERT TO V14.
*/

//----- Global variables (don't touch) -----
double delta_r = ( r_max - r_min ) / ( (double) n_points_r - 1.0 );
double delta_t = pi / ( (double) n_points_t - 1.0 );

double ell_max_total = std::max( ell_max_rotation_off, ell_max_rotation_on );

std::vector<double> r     ( n_points_r );
std::vector<double> t     ( n_points_t );	// theta.
std::vector<double> st    ( n_points_t );	// sin(theta).
std::vector<double> ct    ( n_points_t );	// cos(theta).
std::vector<double> st_dt ( n_points_t );	// sin(theta) dtheta.

std::vector< std::vector<double> > x ( n_points_r, std::vector<double> ( n_points_t ) );
std::vector< std::vector<double> > z ( n_points_r, std::vector<double> ( n_points_t ) );

std::vector< std::vector<double> > P0 ( n_points_t, std::vector<double> ( ell_max_total+1 ) );
std::vector< std::vector<double> > P1 ( n_points_t, std::vector<double> ( ell_max_total+1 ) );

std::vector<double> sqrt_2Lplus1_over_4pi ( ell_max_total+1 );

std::vector< std::vector<double> > B_r_L    ( n_points_r, std::vector<double> ( ell_max_total+1 ) );
std::vector< std::vector<double> > B_1_L    ( n_points_r, std::vector<double> ( ell_max_total+1 ) );
std::vector< std::vector<double> > B_2_L    ( n_points_r, std::vector<double> ( ell_max_total+1 ) );
std::vector< std::vector<double> > E_r_L    ( n_points_r, std::vector<double> ( ell_max_total+1 ) );
std::vector< std::vector<double> > E_1_L    ( n_points_r, std::vector<double> ( ell_max_total+1 ) );
std::vector< std::vector<double> > E_2_L    ( n_points_r, std::vector<double> ( ell_max_total+1 ) );
std::vector< std::vector<double> > B_r_L_dr ( n_points_r, std::vector<double> ( ell_max_total+1 ) );
std::vector< std::vector<double> > B_1_L_dr ( n_points_r, std::vector<double> ( ell_max_total+1 ) );
std::vector< std::vector<double> > B_2_L_dr ( n_points_r, std::vector<double> ( ell_max_total+1 ) );
std::vector< std::vector<double> > E_r_L_dr ( n_points_r, std::vector<double> ( ell_max_total+1 ) );
std::vector< std::vector<double> > E_1_L_dr ( n_points_r, std::vector<double> ( ell_max_total+1 ) );
std::vector< std::vector<double> > E_2_L_dr ( n_points_r, std::vector<double> ( ell_max_total+1 ) );
std::vector< std::vector<double> > B_r_L_dr2( n_points_r, std::vector<double> ( ell_max_total+1 ) );

std::vector< std::vector< std::vector<double> > > B_r ( 3, std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) ) );	// Store the values at the previous few timesteps for use with Adams-Bashforth integration. Don't use new variables for each timestep because otherwise we have to pass the right arrays as arguments to all the functions and this becomes extremely slow. This way, we only need to pass the index n.
std::vector< std::vector< std::vector<double> > > B_t ( 3, std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) ) );
std::vector< std::vector< std::vector<double> > > B_p ( 3, std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) ) );
std::vector< std::vector< std::vector<double> > > E_r ( 3, std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) ) );
std::vector< std::vector< std::vector<double> > > E_t ( 3, std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) ) );
std::vector< std::vector< std::vector<double> > > E_p ( 3, std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) ) );

std::vector< std::vector< std::vector<double> > > B_r_dT ( 3, std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) ) );
std::vector< std::vector< std::vector<double> > > B_t_dT ( 3, std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) ) );
std::vector< std::vector< std::vector<double> > > B_p_dT ( 3, std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) ) );
std::vector< std::vector< std::vector<double> > > E_r_dT ( 3, std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) ) );
std::vector< std::vector< std::vector<double> > > E_t_dT ( 3, std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) ) );
std::vector< std::vector< std::vector<double> > > E_p_dT ( 3, std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) ) );

std::vector< std::vector<double> > curl_B_r     ( n_points_r, std::vector<double> ( n_points_t ) );		// Can assign to single variable within void calculate_time_derivatives() instead of global list to save memory, but want to output to CSV to check values 20230914.
std::vector< std::vector<double> > curl_B_t     ( n_points_r, std::vector<double> ( n_points_t ) );
std::vector< std::vector<double> > curl_B_p     ( n_points_r, std::vector<double> ( n_points_t ) );
std::vector< std::vector<double> > curl_E_r     ( n_points_r, std::vector<double> ( n_points_t ) );
std::vector< std::vector<double> > curl_E_t     ( n_points_r, std::vector<double> ( n_points_t ) );
std::vector< std::vector<double> > curl_E_p     ( n_points_r, std::vector<double> ( n_points_t ) );
std::vector< std::vector<double> > div_E        ( n_points_r, std::vector<double> ( n_points_t ) );
std::vector< std::vector<double> > B_x          ( n_points_r, std::vector<double> ( n_points_t ) );
std::vector< std::vector<double> > B_y          ( n_points_r, std::vector<double> ( n_points_t ) );
std::vector< std::vector<double> > B_z          ( n_points_r, std::vector<double> ( n_points_t ) );
std::vector< std::vector<double> > E_x          ( n_points_r, std::vector<double> ( n_points_t ) );
std::vector< std::vector<double> > E_y          ( n_points_r, std::vector<double> ( n_points_t ) );
std::vector< std::vector<double> > E_z          ( n_points_r, std::vector<double> ( n_points_t ) );
std::vector< std::vector<double> > B_r_dr       ( n_points_r, std::vector<double> ( n_points_t ) );		// Not used by the code, but useful within evaluate_radial_derivatives() for CSV output to check for strange behaviour in dB/dr and dE/dr.
std::vector< std::vector<double> > B_t_dr       ( n_points_r, std::vector<double> ( n_points_t ) );
std::vector< std::vector<double> > B_p_dr       ( n_points_r, std::vector<double> ( n_points_t ) );
std::vector< std::vector<double> > E_r_dr       ( n_points_r, std::vector<double> ( n_points_t ) );
std::vector< std::vector<double> > E_t_dr       ( n_points_r, std::vector<double> ( n_points_t ) );
std::vector< std::vector<double> > E_p_dr       ( n_points_r, std::vector<double> ( n_points_t ) );
std::vector< std::vector<double> > B_r_series   ( n_points_r, std::vector<double> ( n_points_t ) );		// To check where the VSH decomposition is in/accurate.
std::vector< std::vector<double> > B_t_series   ( n_points_r, std::vector<double> ( n_points_t ) );
std::vector< std::vector<double> > B_p_series   ( n_points_r, std::vector<double> ( n_points_t ) );
std::vector< std::vector<double> > E_r_series   ( n_points_r, std::vector<double> ( n_points_t ) );
std::vector< std::vector<double> > E_t_series   ( n_points_r, std::vector<double> ( n_points_t ) );
std::vector< std::vector<double> > E_p_series   ( n_points_r, std::vector<double> ( n_points_t ) );
std::vector< std::vector<double> > abs_dev_B_r  ( n_points_r, std::vector<double> ( n_points_t ) );		// To check where the VSH decomposition is in/accurate.
std::vector< std::vector<double> > abs_dev_B_t  ( n_points_r, std::vector<double> ( n_points_t ) );
std::vector< std::vector<double> > abs_dev_B_p  ( n_points_r, std::vector<double> ( n_points_t ) );
std::vector< std::vector<double> > abs_dev_E_r  ( n_points_r, std::vector<double> ( n_points_t ) );
std::vector< std::vector<double> > abs_dev_E_t  ( n_points_r, std::vector<double> ( n_points_t ) );
std::vector< std::vector<double> > abs_dev_E_p  ( n_points_r, std::vector<double> ( n_points_t ) );
std::vector< std::vector<double> > E_dot_B      ( n_points_r, std::vector<double> ( n_points_t ) );		// To check that force-free conditions are fulfilled.
std::vector< std::vector<double> > Bsq_minus_Esq( n_points_r, std::vector<double> ( n_points_t ) );		// To check that force-free conditions are fulfilled.
std::vector< std::vector<double> > alpha        ( n_points_r, std::vector<double> ( n_points_t ) );		// Don't have to be 2D array but might want to track in CSV file.
std::vector< std::vector<double> > beta         ( n_points_r, std::vector<double> ( n_points_t ) );		// Don't have to be 2D array but might want to track in CSV file.
std::vector< std::vector<double> > magnetic_energy_density ( n_points_r, std::vector<double> ( n_points_t ) );		// Energy density
std::vector< std::vector<double> > friction     ( n_points_r, std::vector<double> ( n_points_t ) );		// Sponge layer at outer boundary to absorb outgoing waves (see Parfrey, 2012, PhD thesis, S3.9).

double r_ghost_in  = r_min - delta_r;
double r_ghost_out = r_max + delta_r;
std::vector<double> B_r_L_ghost_in  ( ell_max_total+1 );
std::vector<double> B_1_L_ghost_in  ( ell_max_total+1 );
std::vector<double> B_2_L_ghost_in  ( ell_max_total+1 );
std::vector<double> E_r_L_ghost_in  ( ell_max_total+1 );
std::vector<double> E_1_L_ghost_in  ( ell_max_total+1 );
std::vector<double> E_2_L_ghost_in  ( ell_max_total+1 );
std::vector<double> B_r_L_ghost_out ( ell_max_total+1 );
std::vector<double> B_1_L_ghost_out ( ell_max_total+1 );
std::vector<double> B_2_L_ghost_out ( ell_max_total+1 );
std::vector<double> E_r_L_ghost_out ( ell_max_total+1 );
std::vector<double> E_1_L_ghost_out ( ell_max_total+1 );
std::vector<double> E_2_L_ghost_out ( ell_max_total+1 );

std::vector<double> integral_r_t_dr ( n_points_r );	// Radial  part of finite element for integrals over r^2 sin(theta) dr dtheta.
std::vector<double> integral_r_t_dt ( n_points_t );	// Angular part of finite element for integrals over r^2 sin(theta) dr dtheta.

double stdev_B_r       = 0;
double stdev_B_t       = 0;
double stdev_B_p       = 0;
double stdev_E_r       = 0;
double stdev_E_t       = 0;
double stdev_E_p       = 0;
double max_abs_dev_B_r = 0;
double max_abs_dev_B_t = 0;
double max_abs_dev_B_p = 0;
double max_abs_dev_E_r = 0;
double max_abs_dev_E_t = 0;
double max_abs_dev_E_p = 0;

double total_magnetic_energy = 0;

int nans_B_r = 0;
int nans_B_t = 0;
int nans_B_p = 0;
int nans_E_r = 0;
int nans_E_t = 0;
int nans_E_p = 0;
int nans_tot = 0;

std::vector<double> VSH_integrand_r ( n_points_t );
std::vector<double> VSH_integrand_1 ( n_points_t );
std::vector<double> VSH_integrand_2 ( n_points_t );

std::ofstream output_file_profiles;
std::ofstream output_file_history;
std::ofstream output_file_VSH_coeffs;
std::ofstream output_file_log;

double time_start_seconds = std::chrono::high_resolution_clock::now().time_since_epoch().count() * 1e-9;
time_t time_start_integer = time( 0 );
char*  time_start_string  = ctime( &time_start_integer );

double T = 0;

int    n_gridpoints_within_force_free_region = 0;
int    num_E_points_changed_for_second_force_free_condition               = 0;	// Output to CSV to keep track of how active the FF condition application is.
double percentage_E_points_points_changed_for_second_force_free_condition = 0;

double Omega                     = 0;	    // Will be updated by the code at each timestep.
double star_rotation_angle       = 0;	    // To keep track of how far the star has rotated.
bool   is_rotation_on            = false;	// Used by functions with different output depending whether rotation is on or off.
int    ell_max_for_this_timestep = 0;		// Set at each timestep to either ell_max_rotation_on or ell_max_rotation_off.




//----- Functions -----


void calculate_gridpoints(){
	
	for( int i=0; i<n_points_r; i++ ){
		r[i] = r_min + i * delta_r;
	}
	for( int j=0; j<n_points_t; j++ ){
		t [j] = t_min + j * delta_t;
		st[j] = sin( t[j] );
		ct[j] = cos( t[j] );
	}
	for( int j=0; j<n_points_t-1; j++ ){
		st_dt[j] = ct[j+1] - ct[j];
	}
	
	for( int i=0; i<n_points_r; i++ ){
		for( int j=0; j<n_points_t; j++ ){
			x[i][j] = r[i] * st[j];
			z[i][j] = r[i] * ct[j];
		}
	}
	
	// Infinitesimals dr and dtheta for integral over r,theta.
	for( int i=0; i<n_points_r-1; i++ ){
		integral_r_t_dr[i] = ( pow(r[i+1],3) - pow(r[i],3) ) / 3.0;
	}
	for( int j=0; j<n_points_t-1; j++ ){
		integral_r_t_dt[j] = ct[j] - ct[j+1];
	}
	
}




void calculate_associated_legendre_functions(){
	
	for( int j=0; j<n_points_t; j++ ){
		
		double x = ct[j];	// For readability.
		
		//--- ell=0 ---
		P0[j][0] = 1.0;
		P0[j][1] = x;
		
		for( int ell=2; ell<=ell_max_total; ell++ ){
			P0[j][ell] = ( P0[j][ell-1] * x * ( 2.0*ell-1.0 ) - P0[j][ell-2] * ( ell-1.0 ) ) / ( (double)ell );
		}
		
		//--- ell=1 ---
		P1[j][1] = sqrt( 1.0 - x*x ) * -1.0;
		P1[j][2] = sqrt( 1.0 - x*x ) * -3.0 * x;
		
		for( int ell=3; ell<=ell_max_total; ell++ ){
			P1[j][ell] = ( P1[j][ell-1] * x * ( 2.0*ell-1.0 ) - P1[j][ell-2] * (double)ell ) / ( ell-1.0 );
		}
		
	}
	
}




void calculate_sqrt_2Lplus1_over_4pi(){
	for( int ell=0; ell<=ell_max_total; ell++ ){
		sqrt_2Lplus1_over_4pi[ell] = sqrt( (2.0*ell+1.0) / ( 4.0*pi ) );
	}
}




double integral_trapezium_sintheta_dtheta( std::vector<double> f ){
	// Trapezium rule for function only of theta, where the infinitesimal is sin(theta) dtheta instead of dtheta.
	
	double ret = 0;
	
	for( int k=1; k<f.size(); k++ ){
		ret += ( f[k-1] + f[k] ) * ( ct[k-1] - ct[k] );
	}
	
	return 0.5 * ret;
}




double integral_2d( std::vector< std::vector<double> > f, std::vector<double> du1, std::vector<double> du2 ){
	// Integral of a 2D function f(u1,u2) given as a list of pre-evaluated numbers, with the area elements du1 and du2 already calculated.
	int    N1  = du1.size();
	int    N2  = du2.size();
	double ret = 0;
	
	for( int i=1; i<N1-1; i++ ){
		ret += ( f[i][0] * du2[0] + f[i].back() * du2[N2-2] ) * ( du1[i-1] + du1[i] );
		for( int j=1; j<N2-1; j++ ){
			ret += f[i][j] * ( ( du1[i-1] + du1[i] ) * ( du2[j-1] + du2[j] ) );
		}
	}
	for( int j=1; j<N2-1; j++ ){
		ret += ( f[0][j] * du1[0] + f.back()[j] * du1[N1-2] ) * ( du2[j-1] + du2[j] );
	}
	ret += ( f[0][0] * du2[0] + f[0].back() * du2[N2-2] ) * du1[0] + ( f.back()[0] * du2[0] + f.back().back() * du2[N2-2] ) * du1[N1-2];
	
	return ret * 0.25;
}




void determine_whether_rotation_on_and_set_ell_max_for_this_timestep( int T_index ){
	
	if( T_index < T_index_Omega_ramp_start ){
		is_rotation_on = false;
		ell_max_for_this_timestep = ell_max_rotation_off;
	}
	else {
		is_rotation_on = true;
		ell_max_for_this_timestep = ell_max_rotation_on;
	}
	
}




void VSH_decomposition( int n,
                        std::vector< std::vector< std::vector<double> > >  A_r, std::vector< std::vector< std::vector<double> > > A_t, std::vector< std::vector< std::vector<double> > > A_p,
                        std::vector< std::vector<double> > &A_r_L, std::vector< std::vector<double> > &A_1_L, std::vector< std::vector<double> > &A_2_L ){
	// First, built lists of integrand values. Then, numerically integrate and put the values into an array.
	// Two separate loops due to two separate ranges on ell.
	
	for( int i=0; i<n_points_r; i++ ){
		
		//--- r-coefficients go from ell=0 ---
		for( int ell=0; ell<=ell_max_for_this_timestep; ell++ ){
			
			for( int j=0; j<n_points_t; j++ ){
				VSH_integrand_r[j] = A_r[n][i][j] * P0[j][ell];
			}
			
			A_r_L[i][ell] = integral_trapezium_sintheta_dtheta( VSH_integrand_r ) * sqrt( (2.0*ell+1.0)*pi );
			
		}
		
		
		//--- (1)- and (2)-coefficients go from ell=1 ---
		for( int ell=1; ell<=ell_max_for_this_timestep; ell++ ){
		
			for( int j=0; j<n_points_t; j++ ){
				VSH_integrand_1[j] = A_t[n][i][j] * P1[j][ell];
				VSH_integrand_2[j] = A_p[n][i][j] * P1[j][ell];
			}
			
			A_1_L[i][ell] = integral_trapezium_sintheta_dtheta( VSH_integrand_1 ) * sqrt( (2.0*ell+1.0)*pi ) / ( (double)ell*(ell+1.0) );
			A_2_L[i][ell] = integral_trapezium_sintheta_dtheta( VSH_integrand_2 ) * sqrt( (2.0*ell+1.0)*pi ) / ( (double)ell*(ell+1.0) );
			
		}
	}
}




void evaluate_radial_derivatives(){
	// Not used by the code, but useful for CSV output to check for strange behaviour in dB/dr and dE/dr.
	// To do <20230824: Second derivatives too.
	
	//--- Reset values to zero ---
	B_r_dr = std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) );
	B_t_dr = std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) );
	B_p_dr = std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) );
	E_r_dr = std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) );
	E_t_dr = std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) );
	E_p_dr = std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) );
	
	//--- Sum over ell for each gridpoint ---
	for( int i=0; i<n_points_r; i++ ){
		for( int j=0; j<n_points_t; j++ ){
			for( int ell=0; ell<=ell_max_for_this_timestep; ell++ ){
		
				double sqrt_factor = sqrt( (2.0*ell+1.0) / ( 4.0 * pi ) );
				
				B_r_dr[i][j] += sqrt_factor * P0[j][ell] * B_r_L_dr[i][ell];
				B_t_dr[i][j] += sqrt_factor * P1[j][ell] * B_1_L_dr[i][ell];
				B_p_dr[i][j] += sqrt_factor * P1[j][ell] * B_2_L_dr[i][ell];
				
				E_r_dr[i][j] += sqrt_factor * P0[j][ell] * E_r_L_dr[i][ell];
				E_t_dr[i][j] += sqrt_factor * P1[j][ell] * E_1_L_dr[i][ell];
				E_p_dr[i][j] += sqrt_factor * P1[j][ell] * E_2_L_dr[i][ell];
				
			}
		}
	}
}




void calculate_Cartesian_vector_components( int n ){
	
	for( int i=0; i<n_points_r; i++ ){
		for( int j=0; j<n_points_t; j++ ){
			
			B_x[i][j] = B_r[n][i][j] * st[j] + B_t[n][i][j] * ct[j];
			B_y[i][j] = B_p[n][i][j];									// Extra variable assignment not necessary, but use for now to avoid confusion.
			B_z[i][j] = B_r[n][i][j] * ct[j] - B_t[n][i][j] * st[j];
			
			E_x[i][j] = E_r[n][i][j] * st[j] + E_t[n][i][j] * ct[j];
			E_y[i][j] = E_p[n][i][j];									// Extra variable assignment not necessary, but use for now to avoid confusion.
			E_z[i][j] = E_r[n][i][j] * ct[j] - E_t[n][i][j] * st[j];
			
		}
	}
}




void apply_initial_field_value_and_VSH_decomposition_for_B(){
	
	for( int i=0; i<n_points_r; i++ ){
		
		for( int j=0; j<n_points_t; j++ ){
			B_r[0][i][j] = B_r_function( r[i], t[j] );
			B_t[0][i][j] = B_t_function( r[i], t[j] );
			B_p[0][i][j] = B_p_function( r[i], t[j] );
		}
		
		for( int L=0; L<=ell_max_for_this_timestep; L++ ){
			
			B_r_L[i][L] = B_r_L_function( r[i], L );
			B_1_L[i][L] = B_1_L_function( r[i], L );
			B_2_L[i][L] = B_2_L_function( r[i], L );
			
			B_r_L_ghost_in[L] = B_r_L_function( r_ghost_in, L );
			B_1_L_ghost_in[L] = B_1_L_function( r_ghost_in, L );
			B_2_L_ghost_in[L] = B_2_L_function( r_ghost_in, L );
			
			B_r_L_ghost_out[L] = B_r_L_function( r_ghost_out, L );
			B_1_L_ghost_out[L] = B_1_L_function( r_ghost_out, L );
			B_2_L_ghost_out[L] = B_2_L_function( r_ghost_out, L );
		}
		
	}
}




void set_field_values_at_ghost_points_for_E(){
	// THIS IS ONLY ACCURATE FOR A NON-ROTATING OBJECT! A NEW CONSIDERATION IS NEEDED WHEN ROTATION STARTS!
	for( int L=0; L<=ell_max_for_this_timestep; L++ ){
			
			E_r_L_ghost_in[L] = 0;
			E_1_L_ghost_in[L] = 0;
			E_2_L_ghost_in[L] = 0;
			
			E_r_L_ghost_out[L] = 0;
			E_1_L_ghost_out[L] = 0;
			E_2_L_ghost_out[L] = 0;
		}
}




void calculate_n_gridpoints_within_force_free_region(){
	
	n_gridpoints_within_force_free_region = 0;
	
	for( int i=0; i<n_points_r; i++ ){
		for( int j=1; j<n_points_t; j++ ){
			if( r[i] * pow( st[j], -2 ) < 1.0 / Omega ){
				n_gridpoints_within_force_free_region++;
			}
		}
	}
}



void apply_force_free_conditions_within_force_free_region( int n ){
	// Update the electric field at all points within the last closed field line of a rotating dipole, except the poles,
	// everywhere the first or second force-free condition is violated, with an arbitrary tolerance for the value of E dot B.
	// Petri 2012, Section 3.2.
	
	calculate_n_gridpoints_within_force_free_region();
	
	num_E_points_changed_for_second_force_free_condition = 0;
	
	for( int i=0; i<n_points_r; i++ ){
		for( int j=1; j<n_points_t-1; j++ ){
			if( r[i] * pow( st[j], -2 ) < 1.0 / Omega ){
				
				double E_dot_B = E_r[n][i][j]*B_r[n][i][j] + E_t[n][i][j]*B_t[n][i][j] + E_p[n][i][j]*B_p[n][i][j];
				double Bsq     = B_r[n][i][j]*B_r[n][i][j] + B_t[n][i][j]*B_t[n][i][j] + B_p[n][i][j]*B_p[n][i][j];
				double Esq     = E_r[n][i][j]*E_r[n][i][j] + E_t[n][i][j]*E_t[n][i][j] + E_p[n][i][j]*E_p[n][i][j];
				
				double E_dot_B_tol = 1e-3;	// Arbitrarily chosen 20230921.
				
				if( ( E_dot_B > E_dot_B_tol ) or ( Bsq - Esq < 0 ) ){
					
					num_E_points_changed_for_second_force_free_condition += 1;
					
					double E_dot_B_over_Bsq = E_dot_B / Bsq;
					double E_prime_r = E_r[n][i][j] - E_dot_B_over_Bsq * B_r[n][i][j];
					double E_prime_t = E_t[n][i][j] - E_dot_B_over_Bsq * B_t[n][i][j];
					double E_prime_p = E_p[n][i][j] - E_dot_B_over_Bsq * B_p[n][i][j];
					double sqrt_Bsq_over_Eprimesq = sqrt( Bsq / ( E_prime_r*E_prime_r + E_prime_t*E_prime_t + E_prime_p*E_prime_p ) );
					
					E_r[n][i][j] = E_prime_r * sqrt_Bsq_over_Eprimesq;
					E_t[n][i][j] = E_prime_t * sqrt_Bsq_over_Eprimesq;
					E_p[n][i][j] = E_prime_p * sqrt_Bsq_over_Eprimesq;
					
				}
			}
		}
	}
	
	// Calculate percentage of gridpoints at which E changed. Below format makes it exactly one decimal place.
	percentage_E_points_points_changed_for_second_force_free_condition = round( (double)1000 * num_E_points_changed_for_second_force_free_condition / n_gridpoints_within_force_free_region ) / 10;
}




void calculate_time_derivatives( int T_index, int n ){
									 
	for( int i=0; i<n_points_r; i++ ){
		for( int j=0; j<n_points_t; j++ ){
			
			//--- Spatial derivatives ---
			div_E   [i][j] = 0;
			curl_E_r[i][j] = 0;
			curl_E_t[i][j] = 0;
			curl_E_p[i][j] = 0;
			curl_B_r[i][j] = 0;
			curl_B_t[i][j] = 0;
			curl_B_p[i][j] = 0;
			
			//- L=0 -
			div_E[i][j] += sqrt_2Lplus1_over_4pi[0] * ( E_r_L_dr[i][0] + 2.0*pow(r[i],-1) * E_r_L[i][0] );
			
			//- L>0 -
			for( int L=1; L<=ell_max_for_this_timestep; L++ ){
				
				div_E[i][j] += sqrt_2Lplus1_over_4pi[L] * ( E_r_L_dr[i][L] + 2.0*pow(r[i],-1) * E_r_L[i][L] - L*(L+1.0)*pow(r[i],-1) * E_1_L[i][L] );
				
				curl_E_r[i][j] += sqrt_2Lplus1_over_4pi[L] * P0[j][L] * -L*(L+1) * E_2_L[i][L] / r[i];
				curl_E_t[i][j] += sqrt_2Lplus1_over_4pi[L] * P1[j][L] * ( -E_2_L_dr[i][L] - E_2_L[i][L] / r[i] );
				curl_E_p[i][j] += sqrt_2Lplus1_over_4pi[L] * P1[j][L] * ( E_1_L_dr[i][L] + ( E_1_L[i][L] - E_r_L[i][L] ) / r[i] );
				
				curl_B_r[i][j] += sqrt_2Lplus1_over_4pi[L] * P0[j][L] * -L*(L+1) * B_2_L[i][L] / r[i];
				curl_B_t[i][j] += sqrt_2Lplus1_over_4pi[L] * P1[j][L] * ( -B_2_L_dr[i][L] - B_2_L[i][L] / r[i] );
				curl_B_p[i][j] += sqrt_2Lplus1_over_4pi[L] * P1[j][L] * ( r[i] * B_r_L_dr2[i][L] + 4.0 * B_r_L_dr[i][L] - (L+2.0)*(L-1.0) * B_r_L[i][L] / r[i] ) / ( (double) L*(L+1.0) );
				
			}
			
			//--- alpha and beta ---
			double B_dot_curl_B = B_r[n][i][j]*curl_B_r   [i][j] + B_t[n][i][j]*curl_B_t   [i][j] + B_p[n][i][j]*curl_B_p   [i][j];
			double E_dot_curl_E = E_r[n][i][j]*curl_E_r   [i][j] + E_t[n][i][j]*curl_E_t   [i][j] + E_p[n][i][j]*curl_E_p   [i][j];
			double B_squared    = B_r[n][i][j]*B_r     [n][i][j] + B_t[n][i][j]*     B_t[n][i][j] + B_p[n][i][j]*     B_p[n][i][j];
			double E_squared    = E_r[n][i][j]*E_r     [n][i][j] + E_t[n][i][j]*     E_t[n][i][j] + E_p[n][i][j]*     E_p[n][i][j];
			
			alpha[i][j] = ( E_dot_curl_E - B_dot_curl_B ) / B_squared;
			beta [i][j] = - div_E[i][j] / B_squared;
			
			//--- E cross B ---
			double E_cross_B_r = E_t[n][i][j] * B_p[n][i][j] - E_p[n][i][j] * B_t[n][i][j];
			double E_cross_B_t = E_p[n][i][j] * B_r[n][i][j] - E_r[n][i][j] * B_p[n][i][j];
			double E_cross_B_p = E_r[n][i][j] * B_t[n][i][j] - E_t[n][i][j] * B_r[n][i][j];
			
			//--- Frictional term, applied to all components except B_r_dT ---
			friction[i][j] = 0;
			if( use_outer_sponge_layer and ( r[i] >= r_sponge ) ){
				friction[i][j] = pow( ( r[i] - r_sponge ) / ( r_max - r_sponge ), friction_beta );
				friction[i][j] = friction_sigma_0 * ( 1.0 - exp( - friction_gamma * friction[i][j] ) );
			}
			
			//--- Time derivatives ---
			int n_for_dT = std::min( abs( n-1 ), 1 );
			
			B_r_dT[n_for_dT][i][j] = - curl_E_r[i][j];
			B_t_dT[n_for_dT][i][j] = - curl_E_t[i][j] - friction[i][j] * B_t[n][i][j];
			B_p_dT[n_for_dT][i][j] = - curl_E_p[i][j] - friction[i][j] * B_p[n][i][j];
			
			E_r_dT[n_for_dT][i][j] = curl_B_r[i][j] + alpha[i][j] * B_r[n][i][j] + beta[i][j] * E_cross_B_r - friction[i][j] * E_r[n][i][j];
			E_t_dT[n_for_dT][i][j] = curl_B_t[i][j] + alpha[i][j] * B_t[n][i][j] + beta[i][j] * E_cross_B_t - friction[i][j] * E_t[n][i][j];
			E_p_dT[n_for_dT][i][j] = curl_B_p[i][j] + alpha[i][j] * B_p[n][i][j] + beta[i][j] * E_cross_B_p - friction[i][j] * E_p[n][i][j];
			
			if( ( T_index >= T_index_Omega_ramp_start ) and ( T_index <= T_index_Omega_ramp_stop ) ){
				E_r_dT[n_for_dT][i][j] += E_r_dT_function( r[i], t[j], dOmega_by_dT );
				E_t_dT[n_for_dT][i][j] += E_t_dT_function( r[i], t[j], dOmega_by_dT );
				E_p_dT[n_for_dT][i][j] += E_p_dT_function( r[i], t[j], dOmega_by_dT );
			}
			
		}
	}
}




void calculate_radial_derivatives_of_VSH_coeffs(){
	
	for( int L=0; L<=ell_max_for_this_timestep; L++ ){
		
		// Inner boundary (using inner ghost points).
		B_r_L_dr[0][L] = ( B_r_L[1][L] - B_r_L_ghost_in[L] ) / ( 2.0*delta_r );
		B_2_L_dr[0][L] = ( B_2_L[1][L] - B_2_L_ghost_in[L] ) / ( 2.0*delta_r );
		E_r_L_dr[0][L] = ( E_r_L[1][L] - E_r_L_ghost_in[L] ) / ( 2.0*delta_r );
		E_1_L_dr[0][L] = ( E_1_L[1][L] - E_1_L_ghost_in[L] ) / ( 2.0*delta_r );
		E_2_L_dr[0][L] = ( E_2_L[1][L] - E_2_L_ghost_in[L] ) / ( 2.0*delta_r );
		
		B_r_L_dr2[0][L] = ( B_r_L[2][L] - 2.0*B_r_L[0][L] + B_r_L_ghost_in[L] ) / ( delta_r*delta_r );
		
		// Symmetric derivatives for all intermediate points.
		for( int i=1; i<n_points_r-1; i++ ){
			B_r_L_dr [i][L] = ( B_r_L[i+1][L] - B_r_L[i-1][L] ) / ( 2.0*delta_r );
			B_2_L_dr [i][L] = ( B_2_L[i+1][L] - B_2_L[i-1][L] ) / ( 2.0*delta_r );
			E_r_L_dr [i][L] = ( E_r_L[i+1][L] - E_r_L[i-1][L] ) / ( 2.0*delta_r );
			E_1_L_dr [i][L] = ( E_1_L[i+1][L] - E_1_L[i-1][L] ) / ( 2.0*delta_r );
			E_2_L_dr [i][L] = ( E_2_L[i+1][L] - E_2_L[i-1][L] ) / ( 2.0*delta_r );
			
			B_r_L_dr2[i][L] = ( B_r_L[i+1][L] - 2.0*B_r_L[i][L] + B_r_L[i-1][L] ) / ( delta_r*delta_r );
		}
		
		// Outewr boundary (using outer ghost points).
		B_r_L_dr.back()[L] = ( B_r_L_ghost_out[L] - B_r_L[n_points_r-2][L] ) / ( 2.0*delta_r );
		B_2_L_dr.back()[L] = ( B_2_L_ghost_out[L] - B_2_L[n_points_r-2][L] ) / ( 2.0*delta_r );
		E_r_L_dr.back()[L] = ( E_r_L_ghost_out[L] - E_r_L[n_points_r-2][L] ) / ( 2.0*delta_r );
		E_1_L_dr.back()[L] = ( E_1_L_ghost_out[L] - E_1_L[n_points_r-2][L] ) / ( 2.0*delta_r );
		E_2_L_dr.back()[L] = ( E_2_L_ghost_out[L] - E_2_L[n_points_r-2][L] ) / ( 2.0*delta_r );
		
		B_r_L_dr2.back()[L] = ( B_r_L_ghost_out[L]  - 2.0*B_r_L.back()[L] + B_r_L[n_points_r-2][L] ) / ( delta_r*delta_r );
		
	}

}




void calculate_B_1_L_and_B_1_L_dr_initial(){
	// This version of the function is for T=0, at which time we DON'T want to update B^{(1),L} because it has been read-in from a header file.
	
	for( int L=1; L<=ell_max_for_this_timestep; L++ ){
		for( int i=0; i<n_points_r; i++ ){
			B_1_L_dr[i][L] = ( r[i] * B_r_L_dr2[i][L] + 3.0 * B_r_L_dr[i][L] ) / ( L*(L+1.0) );
		}
	}
	
}
void calculate_B_1_L_and_B_1_L_dr(){
	// This version of the function is for T>0, at which time we want to update B^{(1),L}.
	
	for( int L=1; L<=ell_max_for_this_timestep; L++ ){
		for( int i=0; i<n_points_r; i++ ){
			B_1_L   [i][L] = ( r[i] * B_r_L_dr [i][L] + 2.0 * B_r_L   [i][L] ) / ( L*(L+1.0) );
			B_1_L_dr[i][L] = ( r[i] * B_r_L_dr2[i][L] + 3.0 * B_r_L_dr[i][L] ) / ( L*(L+1.0) );
		}
	}
	
}
	
	


void calculate_force_free_conditions( int n ){
	// For output to CSV file only.
	for( int i=0; i<n_points_r; i++ ){
		for( int j=0; j<n_points_t; j++ ){
			E_dot_B      [i][j] = E_r[n][i][j]*B_r[n][i][j] + E_t[n][i][j]*B_t[n][i][j] + E_p[n][i][j]*B_p[n][i][j];
			Bsq_minus_Esq[i][j] = B_r[n][i][j]*B_r[n][i][j] + B_t[n][i][j]*B_t[n][i][j] + B_p[n][i][j]*B_p[n][i][j] - E_r[n][i][j]*E_r[n][i][j] - E_t[n][i][j]*E_t[n][i][j] - E_p[n][i][j]*E_p[n][i][j];
		}
	}
}




void integrate_vectors_wrt_time_Adams_Bashforth_Order_1(){
	// For use at the first timestep.
	for( int i=0; i<n_points_r; i++ ){
		for( int j=0; j<n_points_t; j++ ){
			
			B_r[1][i][j] = B_r[0][i][j] + delta_T * B_r_dT[0][i][j];
			B_t[1][i][j] = B_t[0][i][j] + delta_T * B_t_dT[0][i][j];
			B_p[1][i][j] = B_p[0][i][j] + delta_T * B_p_dT[0][i][j];
			E_r[1][i][j] = E_r[0][i][j] + delta_T * E_r_dT[0][i][j];
			E_t[1][i][j] = E_t[0][i][j] + delta_T * E_t_dT[0][i][j];
			E_p[1][i][j] = E_p[0][i][j] + delta_T * E_p_dT[0][i][j];
			
			B_r[0][i][j] = B_r[1][i][j];
			B_t[0][i][j] = B_t[1][i][j];
			B_p[0][i][j] = B_p[1][i][j];
			E_r[0][i][j] = E_r[1][i][j];
			E_t[0][i][j] = E_t[1][i][j];
			E_p[0][i][j] = E_p[1][i][j];
			
		}
	}
}




void integrate_vectors_wrt_time_Adams_Bashforth_Order_2(){
	// For use at the second and all subsequent timesteps.
	for( int i=0; i<n_points_r; i++ ){
		for( int j=0; j<n_points_t; j++ ){
			
			B_r[2][i][j] = B_r[1][i][j] + delta_T * ( 1.5 * B_r_dT[1][i][j] - 0.5 * B_r_dT[0][i][j] );
			B_t[2][i][j] = B_t[1][i][j] + delta_T * ( 1.5 * B_t_dT[1][i][j] - 0.5 * B_t_dT[0][i][j] );
			B_p[2][i][j] = B_p[1][i][j] + delta_T * ( 1.5 * B_p_dT[1][i][j] - 0.5 * B_p_dT[0][i][j] );
			E_r[2][i][j] = E_r[1][i][j] + delta_T * ( 1.5 * E_r_dT[1][i][j] - 0.5 * E_r_dT[0][i][j] );
			E_t[2][i][j] = E_t[1][i][j] + delta_T * ( 1.5 * E_t_dT[1][i][j] - 0.5 * E_t_dT[0][i][j] );
			E_p[2][i][j] = E_p[1][i][j] + delta_T * ( 1.5 * E_p_dT[1][i][j] - 0.5 * E_p_dT[0][i][j] );
			
			B_r[0][i][j] = B_r[1][i][j];
			B_t[0][i][j] = B_t[1][i][j];
			B_p[0][i][j] = B_p[1][i][j];
			E_r[0][i][j] = E_r[1][i][j];
			E_t[0][i][j] = E_t[1][i][j];
			E_p[0][i][j] = E_p[1][i][j];
			
			B_r[1][i][j] = B_r[2][i][j];
			B_t[1][i][j] = B_t[2][i][j];
			B_p[1][i][j] = B_p[2][i][j];
			E_r[1][i][j] = E_r[2][i][j];
			E_t[1][i][j] = E_t[2][i][j];
			E_p[1][i][j] = E_p[2][i][j];
			
			B_r_dT[0][i][j] = B_r_dT[1][i][j];
			B_t_dT[0][i][j] = B_t_dT[1][i][j];
			B_p_dT[0][i][j] = B_p_dT[1][i][j];
			E_r_dT[0][i][j] = E_r_dT[1][i][j];
			E_t_dT[0][i][j] = E_t_dT[1][i][j];
			E_p_dT[0][i][j] = E_p_dT[1][i][j];
			
			
		}
	}
}




void apply_inner_boundary_conditions( int n ){
	// Petri 2012, \S3.3.
	
	if( is_rotation_on ){
		for( int j=0; j<n_points_t; j++ ){
			B_r[n][0][j] = B_r_function( r_min, t[j] );
			E_t[n][0][j] = - Omega * r[0] * st[j] * B_r_function( r_min, t[j] );
			E_p[n][0][j] = 0;
		}
	}
	
	else {
		for( int j=0; j<n_points_t; j++ ){
			B_r[n][0][j] = B_r_function( r_min, t[j] );
			B_t[n][0][j] = B_t_function( r_min, t[j] );
			B_p[n][0][j] = B_p_function( r_min, t[j] );
			E_r[n][0][j] = 0;
			E_p[n][0][j] = 0;
			E_t[n][0][j] = 0;
		}
	}

}




void apply_outer_boundary_conditions( int n ){
	// Petri 2012, \S3.3.
	
	if( is_rotation_on ){
		for( int j=0; j<n_points_t; j++ ){
			B_t[n].back()[j] = 0.5 * ( B_t[n].back()[j] - E_p[n].back()[j] );
			B_p[n].back()[j] = 0.5 * ( E_t[n].back()[j] + B_p[n].back()[j] );
			E_t[n].back()[j] =   B_p[n].back()[j];
			E_p[n].back()[j] = - B_t[n].back()[j];
		}
	}
	
	else{
		for( int j=0; j<n_points_t; j++ ){
			B_r[n].back()[j] = B_r_function( r_max, t[j] );
			B_t[n].back()[j] = B_t_function( r_max, t[j] );
			B_p[n].back()[j] = B_p_function( r_max, t[j] );
			E_r[n].back()[j] = 0;
			E_p[n].back()[j] = 0;
			E_t[n].back()[j] = 0;
		}
	}
	
}




void calculate_stdev( int n ){
	
	//--- Reset values to zero ---
	B_r_series = std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) );
	B_t_series = std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) );
	B_p_series = std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) );
	E_r_series = std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) );
	E_t_series = std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) );
	E_p_series = std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) );
	stdev_B_r = 0;
	stdev_B_t = 0;
	stdev_B_p = 0;
	stdev_E_r = 0;
	stdev_E_t = 0;
	stdev_E_p = 0;
	
	//--- Evaluate VSH series by summing over ell for each gridpoint ---
	for( int i=0; i<n_points_r; i++ ){
		for( int j=0; j<n_points_t; j++ ){
			for( int ell=0; ell<=ell_max_for_this_timestep; ell++ ){
				
				B_r_series[i][j] += sqrt_2Lplus1_over_4pi[ell] * P0[j][ell] * B_r_L[i][ell];
				B_t_series[i][j] += sqrt_2Lplus1_over_4pi[ell] * P1[j][ell] * B_1_L[i][ell];
				B_p_series[i][j] += sqrt_2Lplus1_over_4pi[ell] * P1[j][ell] * B_2_L[i][ell];
				
				E_r_series[i][j] += sqrt_2Lplus1_over_4pi[ell] * P0[j][ell] * E_r_L[i][ell];
				E_t_series[i][j] += sqrt_2Lplus1_over_4pi[ell] * P1[j][ell] * E_1_L[i][ell];
				E_p_series[i][j] += sqrt_2Lplus1_over_4pi[ell] * P1[j][ell] * E_2_L[i][ell];
				
			}
		}
	}
	
	//--- Calculate stdev and build lists of square deviations along the way ---
	for( int i=0; i<n_points_r; i++ ){
		for( int j=0; j<n_points_t; j++ ){
			
			abs_dev_B_r[i][j] = abs( B_r[n][i][j] - B_r_series[i][j] );
			abs_dev_B_t[i][j] = abs( B_t[n][i][j] - B_t_series[i][j] );
			abs_dev_B_p[i][j] = abs( B_p[n][i][j] - B_p_series[i][j] );
			abs_dev_E_r[i][j] = abs( E_r[n][i][j] - E_r_series[i][j] );
			abs_dev_E_t[i][j] = abs( E_t[n][i][j] - E_t_series[i][j] );
			abs_dev_E_p[i][j] = abs( E_p[n][i][j] - E_p_series[i][j] );
			
			stdev_B_r += pow( abs_dev_B_r[i][j], 2 );
			stdev_B_t += pow( abs_dev_B_t[i][j], 2 );
			stdev_B_p += pow( abs_dev_B_p[i][j], 2 );
			stdev_E_r += pow( abs_dev_E_r[i][j], 2 );
			stdev_E_t += pow( abs_dev_E_t[i][j], 2 );
			stdev_E_p += pow( abs_dev_E_p[i][j], 2 );
			
			if( abs_dev_B_r[i][j] > max_abs_dev_B_r ){ max_abs_dev_B_r = abs_dev_B_r[i][j]; }
			if( abs_dev_B_t[i][j] > max_abs_dev_B_t ){ max_abs_dev_B_t = abs_dev_B_t[i][j]; }
			if( abs_dev_B_p[i][j] > max_abs_dev_B_p ){ max_abs_dev_B_p = abs_dev_B_p[i][j]; }
			if( abs_dev_E_r[i][j] > max_abs_dev_E_r ){ max_abs_dev_E_r = abs_dev_E_r[i][j]; }
			if( abs_dev_E_t[i][j] > max_abs_dev_E_t ){ max_abs_dev_E_t = abs_dev_E_t[i][j]; }
			if( abs_dev_E_p[i][j] > max_abs_dev_E_p ){ max_abs_dev_E_p = abs_dev_E_p[i][j]; }
			
		}
	}
	
	stdev_B_r = sqrt( stdev_B_r / ( (double) n_points_r * n_points_t ) );
	stdev_B_t = sqrt( stdev_B_t / ( (double) n_points_r * n_points_t ) );
	stdev_B_p = sqrt( stdev_B_p / ( (double) n_points_r * n_points_t ) );
	stdev_E_r = sqrt( stdev_E_r / ( (double) n_points_r * n_points_t ) );
	stdev_E_t = sqrt( stdev_E_t / ( (double) n_points_r * n_points_t ) );
	stdev_E_p = sqrt( stdev_E_p / ( (double) n_points_r * n_points_t ) );
	
}




void calculate_magnetic_energy( int n ){
	
	for( int i=0; i<n_points_r; i++ ){
		for( int j=0; j<n_points_t; j++ ){
			magnetic_energy_density[i][j] = pow(B_r[n][i][j],2) + pow(B_t[n][i][j],2) + pow(B_p[n][i][j],2);
		}
	}
	
	total_magnetic_energy = integral_2d( magnetic_energy_density, integral_r_t_dr, integral_r_t_dt ) * 2.0 * pi;
	
}




double estimate_csv_profiles_filesize(){
	
	int n_columns = 63;		// You are responsible for keeping this value updated.
	int n_rows = 0;
	double fudge_factor = 11.0;
	
	for( int T_index=0; T_index<=n_timesteps; T_index++ ){
		bool condition_T_1 = ( T_index >= csv_profiles_write_T_min ) and ( T_index <= csv_profiles_write_T_max );
		bool condition_T_2 = ( T_index % csv_profiles_write_freq_T == 0 ) or ( T_index == n_timesteps );
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
	
	output_file_profiles << "T_index,T,T_sec,i,j,r,t,x,z";
	output_file_profiles << ",B_x,B_y,B_z,B,E_x,E_y,E_z,E";
	output_file_profiles << ",alpha,beta";
	output_file_profiles << ",E_dot_B,Bsq_minus_Esq,sgn(Bsq_minus_Esq)";
	output_file_profiles << ",B_r_dr,B_t_dr,B_p_dr,E_r_dr,E_t_dr,E_p_dr";
	output_file_profiles << ",B_r,B_t,B_p,E_r,E_t,E_p";
	output_file_profiles << ",B_r_series,B_t_series,B_p_series,E_r_series,E_t_series,E_p_series";
	output_file_profiles << ",abs_dev_B_r,abs_dev_B_t,abs_dev_B_p,abs_dev_E_r,abs_dev_E_t,abs_dev_E_p";
	output_file_profiles << ",magnetic_energy_density";
	output_file_profiles << ",curl_B_r,curl_B_t,curl_B_p";
	output_file_profiles << ",B_r_dT,B_t_dT,B_p_dT,E_r_dT,E_t_dT,E_p_dT";
	output_file_profiles << ",B_r_dT_over_B_r,B_t_dT_over_B_t,B_p_dT_over_B_p,E_r_dT_over_E_r,E_t_dT_over_E_t,E_p_dT_over_E_p";
	output_file_profiles << ",friction";
	output_file_profiles << "\n";
}


int sgn( double x ){
	//signum function with sgn(x)=0 for x=0.
	return ( x > 0 ) - ( x < 0 );
}

void output_to_csv_profiles( int T_index, int n ){
	
	bool condition_T_1 = ( T_index >= csv_profiles_write_T_min ) and ( T_index <= csv_profiles_write_T_max );
	bool condition_T_2 = ( T_index % csv_profiles_write_freq_T == 0 ) or ( T_index == n_timesteps );
	
	if( condition_T_1 and condition_T_2 ){
	
		for( int i=csv_profiles_write_i_min; i<n_points_r; i++ ){
			for( int j=0; j<n_points_t; j++ ){
				
				bool condition_r = ( i <= csv_profiles_write_i_max ) and ( ( i % csv_profiles_write_freq_r == 0 ) or ( i == n_points_r - 1 ) );
				bool condition_t = ( j % csv_profiles_write_freq_t == 0 ) or ( j == n_points_t - 1 );
		
				if( condition_r and condition_t ){
					
					double B = sqrt( pow( B_r[n][i][j], 2 ) + pow( B_t[n][i][j], 2 ) + pow( B_p[n][i][j], 2 ) );
					double E = sqrt( pow( E_r[n][i][j], 2 ) + pow( E_t[n][i][j], 2 ) + pow( E_p[n][i][j], 2 ) );
					
					double tolerance_for_zero = 1e-13;
					double B_r_dT_over_B_r = 0;
					double B_t_dT_over_B_t = 0;
					double B_p_dT_over_B_p = 0;
					double E_r_dT_over_E_r = 0;
					double E_t_dT_over_E_t = 0;
					double E_p_dT_over_E_p = 0;
					if( abs( B_r[n][i][j] ) > tolerance_for_zero ){ B_r_dT_over_B_r = B_r_dT[n][i][j] / B_r[n][i][j]; }
					if( abs( B_t[n][i][j] ) > tolerance_for_zero ){ B_t_dT_over_B_t = B_t_dT[n][i][j] / B_t[n][i][j]; }
					if( abs( B_p[n][i][j] ) > tolerance_for_zero ){ B_p_dT_over_B_p = B_p_dT[n][i][j] / B_p[n][i][j]; }
					if( abs( E_r[n][i][j] ) > tolerance_for_zero ){ E_r_dT_over_E_r = E_r_dT[n][i][j] / E_r[n][i][j]; }
					if( abs( E_t[n][i][j] ) > tolerance_for_zero ){ E_t_dT_over_E_t = E_t_dT[n][i][j] / E_t[n][i][j]; }
					if( abs( E_p[n][i][j] ) > tolerance_for_zero ){ E_p_dT_over_E_p = E_p_dT[n][i][j] / E_p[n][i][j]; }
					
					
					output_file_profiles <<std::setprecision(csv_precision)<< T_index <<","<< T <<","<< T*T_factor <<","<< i <<","<< j <<","<< r[i] <<","<< t[j] <<","<< x[i][j] <<","<< z[i][j]
					                                                       <<","<< B_x[i][j] <<","<< B_y[i][j] <<","<< B_z[i][j] <<","<< B
																		   <<","<< E_x[i][j] <<","<< E_y[i][j] <<","<< E_z[i][j] <<","<< E
																		   <<","<< alpha[i][j] <<","<< beta[i][j]
																		   <<","<< E_dot_B[i][j] <<","<< Bsq_minus_Esq[i][j] <<","<< sgn( Bsq_minus_Esq[i][j] )
																		   <<","<< B_r_dr[i][j] <<","<< B_t_dr[i][j] <<","<<  B_p_dr[i][j] <<","<< E_r_dr[i][j] <<","<< E_t_dr[i][j] <<","<<  E_p_dr[i][j]
																		   <<","<< B_r    [n][i][j] <<","<< B_t    [n][i][j] <<","<< B_p    [n][i][j] <<","<< E_r    [n][i][j] <<","<< E_t    [n][i][j] <<","<< E_p    [n][i][j]
																		   <<","<< B_r_series[i][j] <<","<< B_t_series[i][j] <<","<< B_p_series[i][j] <<","<< E_r_series[i][j] <<","<< E_t_series[i][j] <<","<< E_p_series[i][j]
																		   <<","<< abs_dev_B_r [i][j] <<","<< abs_dev_B_t [i][j] <<","<< abs_dev_B_p [i][j] <<","<< abs_dev_E_r [i][j] <<","<< abs_dev_E_t [i][j] <<","<< abs_dev_E_p [i][j]
																		   <<","<< magnetic_energy_density[i][j]
																		   <<","<< curl_B_r[i][j] <<","<< curl_B_t[i][j] <<","<< curl_B_p[i][j]
																		   <<","<< B_r_dT[n][i][j] <<","<< B_t_dT[n][i][j] <<","<< B_p_dT[n][i][j] <<","<< E_r_dT[n][i][j] <<","<< E_t_dT[n][i][j] <<","<< E_p_dT[n][i][j]
																		   <<","<< B_r_dT_over_B_r <<","<< B_t_dT_over_B_t <<","<< B_p_dT_over_B_p <<","<< E_r_dT_over_E_r <<","<< E_t_dT_over_E_t <<","<< E_p_dT_over_E_p
																		   <<","<< friction[i][j]
																		   << "\n";
					
				}
			}
		}
	}
}




void count_nans( int n ){
	nans_B_r = 0;
	nans_B_t = 0;
	nans_B_p = 0;
	nans_E_r = 0;
	nans_E_t = 0;
	nans_E_p = 0;
	
	for( int i=0; i<n_points_r; i++ ){
		for( int j=0; j<n_points_t; j++ ){
			if( ( isnan( B_r[n][i][j] ) ) or ( isinf( B_r[n][i][j] ) ) ){ nans_B_r += 1; }
			if( ( isnan( B_t[n][i][j] ) ) or ( isinf( B_t[n][i][j] ) ) ){ nans_B_t += 1; }
			if( ( isnan( B_p[n][i][j] ) ) or ( isinf( B_p[n][i][j] ) ) ){ nans_B_p += 1; }
			if( ( isnan( E_r[n][i][j] ) ) or ( isinf( E_r[n][i][j] ) ) ){ nans_E_r += 1; }
			if( ( isnan( E_t[n][i][j] ) ) or ( isinf( E_t[n][i][j] ) ) ){ nans_E_t += 1; }
			if( ( isnan( E_p[n][i][j] ) ) or ( isinf( E_p[n][i][j] ) ) ){ nans_E_p += 1; }
		}
	}
	
	nans_tot = nans_B_r + nans_B_t + nans_B_p + nans_E_r + nans_E_t + nans_E_p;
}




void output_headers_to_csv_history(){
	// Split over separate lines for readability and ease of checking consistency with output_to_csv_history().
	
	output_file_history << "T_index,T,T_sec";
	output_file_history << ",Omega,R_LC,star_rotation_angle";
	output_file_history << ",ell_max";
	output_file_history << ",E_pts_chngd,pct_E_pts_chngd";
	output_file_history << ",stdev_B_r,stdev_B_t,stdev_B_p,stdev_E_r,stdev_E_t,stdev_E_p";
	output_file_history << ",max_abs_dev_B_r,max_abs_dev_B_t,max_abs_dev_B_p,max_abs_dev_E_r,max_abs_dev_E_t,max_abs_dev_E_p";
	output_file_history << ",nans_B_r,nans_B_t,nans_B_p,nans_E_r,nans_E_t,nans_E_p,nans_tot";
	output_file_history << ",total_magnetic_energy";
	output_file_history << "\n";
}




void output_to_csv_history( int T_index ){
	
	bool condition_T = ( T_index % csv_history_write_freq_T == 0 ) or ( T_index == n_timesteps );
	
	if( condition_T ){
		
		double R_LC = 0;
		if( Omega > 0 ){
			R_LC = 1.0 / Omega;
		}
		
		output_file_history <<std::setprecision(csv_precision)<< T_index <<","<< T <<","<< T*T_factor
		                                                      <<","<< Omega <<","<< R_LC <<","<< star_rotation_angle
		                                                      <<","<< ell_max_for_this_timestep
															  <<","<< num_E_points_changed_for_second_force_free_condition <<","<< percentage_E_points_points_changed_for_second_force_free_condition
															  <<","<< stdev_B_r     <<","<< stdev_B_t     <<","<< stdev_B_p     <<","<< stdev_E_r     <<","<< stdev_E_t     <<","<< stdev_E_p
															  <<","<< max_abs_dev_B_r <<","<< max_abs_dev_B_t <<","<< max_abs_dev_B_p <<","<< max_abs_dev_E_r <<","<< max_abs_dev_E_t <<","<< max_abs_dev_E_p
															  <<","<< nans_B_r <<","<< nans_B_t <<","<< nans_B_p <<","<< nans_E_r <<","<< nans_E_t <<","<< nans_E_p <<","<< nans_tot
															  <<","<< total_magnetic_energy
															  <<"\n";
							
	}
}




void output_headers_to_csv_VSH_coeffs(){
	// Split over separate lines for readability and ease of checking consistency with output_to_csv_history().
	
	output_file_VSH_coeffs << "T_index,i,r,ell,B_r_L,B_1_L,B_2_L,E_r_L,E_1_L,E_2_L\n";
	
}




void output_to_csv_VSH_coeffs( int T_index ){
	
	bool condition_T = ( T_index % csv_VSH_coeffs_write_freq_T == 0 ) or ( T_index == n_timesteps );
	
	if( condition_T ){
		for( int i=0; i<n_points_r; i++ ){
			for( int ell=0; ell<=ell_max_for_this_timestep; ell++ ){
				output_file_VSH_coeffs <<std::setprecision(csv_precision)<< T_index <<","<< i <<","<< r[i] <<","<< ell
				                                                         <<","<< B_r_L[i][ell] <<","<< B_1_L[i][ell] <<","<< B_2_L[i][ell]
																		 <<","<< E_r_L[i][ell] <<","<< E_1_L[i][ell] <<","<< E_2_L[i][ell]
																		 <<"\n";
			}
		}
	}
	
}




void output_headers_to_screen(){
	std::cout <<std::left<<std::setw(w)<< "T_index" <<std::left<<std::setw(w)<< "T" <<"|\t"
	          <<std::left<<std::setw(w)<< "B_r" <<std::left<<std::setw(w)<< "B_t" <<std::left<<std::setw(w)<< "B_p" <<"|\t"
			  <<std::left<<std::setw(w)<< "E_r" <<std::left<<std::setw(w)<< "E_t" <<std::left<<std::setw(w)<< "E_p" <<"|\t"
			  <<std::left<<std::setw(w)<< "E dot B" <<std::left<<std::setw(w)<< "B^2 - E^2"
			  <<std::left<<std::setw(w)<< "% E pts chngd" <<"|\t"
			  <<std::left<<std::setw(w)<< "nans_tot" <<std::left<<std::setw(w)<< "% done" <<std::left<<std::setw(w)<< "Est. time left"
			  << std::endl;
}




void output_to_screen( int T_index, int n ){
	
	if( ( T_index % cout_freq_T == 0 ) or ( T_index == n_timesteps ) ){
	
		double percent_done = round( (double)1000*T_index/n_timesteps )/10; // This format makes it exactly one decimal place.
		
		std::cout <<std::left<<std::setw(w)<< T_index <<std::left<<std::setw(w)<< T <<"|\t"
				  <<std::left<<std::setw(w)<< B_r[n][cout_i][cout_j] <<std::left<<std::setw(w)<< B_t[n][cout_i][cout_j] <<std::left<<std::setw(w)<< B_p[n][cout_i][cout_j] <<"|\t"
				  <<std::left<<std::setw(w)<< E_r[n][cout_i][cout_j] <<std::left<<std::setw(w)<< E_t[n][cout_i][cout_j] <<std::left<<std::setw(w)<< E_p[n][cout_i][cout_j] <<"|\t"
				  <<std::left<<std::setw(w)<< E_dot_B[cout_i][cout_j] <<std::left<<std::setw(w)<< Bsq_minus_Esq[cout_i][cout_j]
				  <<std::left<<std::setw(w)<< percentage_E_points_points_changed_for_second_force_free_condition <<"|\t"
				  <<std::left<<std::setw(w)<< nans_tot <<std::left<<std::setw(w)<< percent_done;
		
		if( T_index == 0 ){
			std::cout << std::endl;
		} else {
			double time_now_seconds = std::chrono::high_resolution_clock::now().time_since_epoch().count() * 1e-9;
			double runtime = time_now_seconds - time_start_seconds;
			double time_left = ( (double) n_timesteps / T_index - 1.0 ) * runtime;
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




void output_parameters_to_log_file(){
	
	output_file_log << log_file_comment << "\n\n";
	
	output_file_log << "Time start:\t" << time_start_string << "\n";
	
	output_file_log << "Parameter values (given in code units, then in SI units)\n";
	output_file_log << "Omega_max              :\t" <<std::left<<std::setw(16)<< Omega_max << "\n";
	output_file_log << "R_LC_max               :\t" <<std::left<<std::setw(16)<< 1.0 / Omega_max << "\n";
	output_file_log << "r_min, r_max           :\t" <<std::left<<std::setw(16)<< r_min   <<std::left<<std::setw(16)<< r_max << "\n";
	output_file_log << "n_points_r, delta_r    :\t" <<std::left<<std::setw(16)<< n_points_r <<std::left<<std::setw(16)<< delta_r << std::endl;
	output_file_log << "n_points_t, delta_t    :\t" <<std::left<<std::setw(16)<< n_points_t <<std::left<<std::setw(16)<< delta_t << std::endl;
	output_file_log << "Timestep               :\t" <<std::left<<std::setw(16)<< delta_T <<std::left<<std::setw(16)<< delta_T * T_factor << "\n";
	output_file_log << "CFL max timestep       :\t" <<std::left<<std::setw(16)<< delta_r <<std::left<<std::setw(16)<< delta_r * T_factor << "\n";
	output_file_log << "Length of simulation   :\t" <<std::left<<std::setw(16)<< T_max << T_max * T_factor <<"\tsteps: " << n_timesteps << "\n";
	output_file_log << "Est. CSV size (MB)     :\t" << estimate_csv_profiles_filesize() << "\n";
	output_file_log << "ell_max_rotation_off/on:\t" << ell_max_rotation_off <<"\t"<< ell_max_rotation_on << "\n";
	output_file_log << "Printing values at ( r[" << cout_i <<"], theta[" << cout_j <<"] ) = ( " << r[cout_i] <<", " << t[cout_j] << " ) = ( " << r[cout_i] <<", " << t[cout_j]/pi << " pi )." << "\n";
	output_file_log << "dOmega_by_dT_index     :\t" << dOmega_by_dT_index <<"\tdOmega_by_dT:\t"<< dOmega_by_dT << "\n";
	output_file_log << "Ramp start             :\t" << T_index_Omega_ramp_start <<"\t"<< T_index_Omega_ramp_start * delta_T << "\n";
	output_file_log << "Ramp stop              :\t" << T_index_Omega_ramp_stop  <<"\t"<< T_index_Omega_ramp_stop  * delta_T << "\n";
	
	output_file_log << "\nuse_outer_sponge_layer:\t" << use_outer_sponge_layer << "\n";
	output_file_log << "sigma_0, gamma, beta    :\t" << friction_sigma_0 <<"\t"<< friction_gamma <<"\t"<< friction_beta << "\n";
	
	output_file_log << "\nCSV file saved:\tCSV/" << output_filename << "_1_profiles.csv"   << "\n";
	output_file_log <<   "CSV file saved:\tCSV/" << output_filename << "_2_history.csv"    << "\n";
	output_file_log <<   "CSV file saved:\tCSV/" << output_filename << "_3_VSH_coeffs.csv" << "\n";
}




void output_execution_time(){
	
	double time_stop_seconds = std::chrono::high_resolution_clock::now().time_since_epoch().count() * 1e-9;
	time_t time_stop_integer = time( 0 );
	char*  time_stop_string  = ctime( &time_stop_integer );

	double exec_time   = time_stop_seconds - time_start_seconds;
	int    exec_time_h = exec_time / 3600;
	int    exec_time_m = exec_time / 60 - exec_time_h * 60;
	int    exec_time_s = exec_time - exec_time_h * 3600 - exec_time_m * 60;
	
	std::cout       <<                              "\nTime stop:\t" << time_stop_string << "Execution time (h:mm:ss):\t" << exec_time_h;
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