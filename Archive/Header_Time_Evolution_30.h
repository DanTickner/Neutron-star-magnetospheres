/*
Header_Time_Evolution_30.h

Code snippets for the Time_Evolution_xx.cpp codes. Keeps them shorter and allows for two separate for-loops with Euler and Adams-Bashforth integration
without copying code.

V30: Add boundary condition on d/dT B at the inner surface, at the cost of the boundary condition on B_phi at the inner surface that I used in V29.
*/

//----- Global variables (don't touch) -----
double delta_r = ( r_max - r_min ) / ( (double) n_points_r - 1.0 );
double delta_t = pi / ( (double) n_points_t - 1.0 );

std::vector<double> r     ( n_points_r );
std::vector<double> t     ( n_points_t );	// theta.
std::vector<double> st    ( n_points_t );	// sin(theta).
std::vector<double> ct    ( n_points_t );	// cos(theta).
std::vector<double> st_dt ( n_points_t );	// sin(theta) dtheta.

std::vector< std::vector<double> > x ( n_points_r, std::vector<double> ( n_points_t ) );
std::vector< std::vector<double> > z ( n_points_r, std::vector<double> ( n_points_t ) );

std::vector< std::vector<double> > P0 ( n_points_t, std::vector<double> ( ell_max+1 ) );
std::vector< std::vector<double> > P1 ( n_points_t, std::vector<double> ( ell_max+1 ) );

std::vector<double> sqrt_2Lplus1_over_4pi        ( ell_max+1 );
std::vector<double> sqrt_2Lplus1_pi              ( ell_max+1 );
std::vector<double> sqrt_2Lplus1_pi_over_L_Lplus1( ell_max+1 );

std::vector< std::vector< std::vector<double> > > B_VSH_coeffs     ( 3, std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( ell_max+1 ) ) );
std::vector< std::vector< std::vector<double> > > B_VSH_coeffs_dr  ( 3, std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( ell_max+1 ) ) );
std::vector< std::vector< std::vector<double> > > B_VSH_coeffs_dr2 ( 3, std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( ell_max+1 ) ) );
std::vector< std::vector< std::vector<double> > > E_VSH_coeffs     ( 3, std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( ell_max+1 ) ) );
std::vector< std::vector< std::vector<double> > > E_VSH_coeffs_dr  ( 3, std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( ell_max+1 ) ) );
std::vector< std::vector< std::vector<double> > > E_VSH_coeffs_dr2 ( 3, std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( ell_max+1 ) ) );

std::vector< std::vector< std::vector<double> > > B             ( 3, std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) ) );
std::vector< std::vector< std::vector<double> > > B_minus1      ( 3, std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) ) );
std::vector< std::vector< std::vector<double> > > B_minus2      ( 3, std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) ) );
std::vector< std::vector< std::vector<double> > > B_minus3      ( 3, std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) ) );
std::vector< std::vector< std::vector<double> > > B_dT          ( 3, std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) ) );
std::vector< std::vector< std::vector<double> > > B_dT_minus1   ( 3, std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) ) );
std::vector< std::vector< std::vector<double> > > B_dT_minus2   ( 3, std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) ) );
std::vector< std::vector< std::vector<double> > > B_dT_minus3   ( 3, std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) ) );
std::vector< std::vector< std::vector<double> > > B_curl        ( 3, std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) ) );
std::vector< std::vector< std::vector<double> > > B_Cartesian   ( 3, std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) ) );
std::vector< std::vector< std::vector<double> > > B_dr          ( 3, std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) ) );
std::vector< std::vector< std::vector<double> > > B_dr2         ( 3, std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) ) );
std::vector< std::vector< std::vector<double> > > B_series      ( 3, std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) ) );	// Values of the VSH decomposition of the field components.
std::vector< std::vector< std::vector<double> > > B_abs_rel_dev ( 3, std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) ) );

std::vector< std::vector< std::vector<double> > > E             ( 3, std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) ) );
std::vector< std::vector< std::vector<double> > > E_minus1      ( 3, std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) ) );
std::vector< std::vector< std::vector<double> > > E_minus2      ( 3, std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) ) );
std::vector< std::vector< std::vector<double> > > E_minus3      ( 3, std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) ) );
std::vector< std::vector< std::vector<double> > > E_dT          ( 3, std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) ) );
std::vector< std::vector< std::vector<double> > > E_dT_minus1   ( 3, std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) ) );
std::vector< std::vector< std::vector<double> > > E_dT_minus2   ( 3, std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) ) );
std::vector< std::vector< std::vector<double> > > E_dT_minus3   ( 3, std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) ) );
std::vector< std::vector< std::vector<double> > > E_curl        ( 3, std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) ) );
std::vector< std::vector< std::vector<double> > > E_Cartesian   ( 3, std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) ) );
std::vector< std::vector< std::vector<double> > > E_dr          ( 3, std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) ) );
std::vector< std::vector< std::vector<double> > > E_dr2         ( 3, std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) ) );
std::vector< std::vector< std::vector<double> > > E_series      ( 3, std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) ) );	// Values of the VSH decomposition of the field components.
std::vector< std::vector< std::vector<double> > > E_abs_rel_dev ( 3, std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) ) );

std::vector< std::vector< std::vector<double> > > J             ( 3, std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) ) );	// Current density vector.
std::vector< std::vector< std::vector<double> > > J_Cartesian   ( 3, std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) ) );


std::vector< std::vector<double> > B_squared                     ( n_points_r, std::vector<double> ( n_points_t ) );
std::vector< std::vector<double> > E_squared                     ( n_points_r, std::vector<double> ( n_points_t ) );
std::vector< std::vector<double> > B_squared_minus_E_squared     ( n_points_r, std::vector<double> ( n_points_t ) );
std::vector< std::vector<double> > E_dot_B                       ( n_points_r, std::vector<double> ( n_points_t ) );
std::vector< std::vector<double> > E_div                         ( n_points_r, std::vector<double> ( n_points_t ) );
std::vector< std::vector<double> > magnetic_energy_density       ( n_points_r, std::vector<double> ( n_points_t ) );		// Energy density
std::vector< std::vector<double> > friction                      ( n_points_r, std::vector<double> ( n_points_t ) );		// Sponge layer at outer boundary to absorb outgoing waves (see Parfrey, 2012, PhD thesis, S3.9).
std::vector< std::vector<int   > > is_gridpoint_within_FF_region ( n_points_r, std::vector<int   > ( n_points_t ) );
std::vector< std::vector<int   > > were_FF_conditions_applied    ( n_points_r, std::vector<int   > ( n_points_t ) );		// +1: FF conditions applied here. -1: FF conditions not applied here. 0: Outside FF region.

std::vector<double> integral_r_t_dr ( n_points_r );	// Radial  part of finite element for integrals over r^2 sin(theta) dr dtheta.
std::vector<double> integral_r_t_dt ( n_points_t );	// Angular part of finite element for integrals over r^2 sin(theta) dr dtheta.

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

std::ofstream output_file_profiles;
std::ofstream output_file_history;
std::ofstream output_file_VSH_coeffs;
std::ofstream output_file_BCs;
std::ofstream output_file_log;

double time_start_seconds = std::chrono::high_resolution_clock::now().time_since_epoch().count() * 1e-9;
time_t time_start_integer = time( 0 );
char*  time_start_string  = ctime( &time_start_integer );

int    T_index = 0;
double T       = 0;

int    n_gridpoints_within_FF_region                              = 0;
int    num_E_points_changed_for_second_FF_condition               = 0;	// Output to CSV to keep track of how active the FF condition application is.
double percentage_E_points_points_changed_for_second_FF_condition = 0;

double Omega                     = 0;	    // Will be updated by the code at each timestep.
double Omega_prev                = 0;		// The value of Omega at the previous timestep.
double star_rotation_angle       = 0;	    // To keep track of how far the star has rotated.
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
		
		//--- m = 0 ---
		P0[j][0] = 1.0;
		P0[j][1] = x;
		
		for( int ell=1; ell<ell_max; ell++ ){
			P0[j][ell+1] = ( ( 2.0*ell+1.0 ) * x * P0[j][ell] - ell * P0[j][ell-1] ) / ( (double) ell+1.0 );
		}
		
		//--- m = 1 ---
		P1[j][1] = sqrt( 1.0 - x*x ) * -1.0;
		P1[j][2] = sqrt( 1.0 - x*x ) * -3.0 * x;
		
		for( int ell=2; ell<ell_max; ell++ ){
			P1[j][ell+1] = ( ( 2.0*ell+1.0 ) * x * P1[j][ell] - ( ell+1.0 ) * P1[j][ell-1] ) / ( (double) ell );
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




void apply_inner_boundary_conditions( bool record_before_and_after = false ){
	// Petri 2012, \S3.3.
	
	//--- Record values before inner BCs applied ---
	if( record_before_and_after ){
		for( int j=0; j<n_points_t; j++ ){
			for( int f=0; f<3; f++ ){
				B_inner_before_BCs[f][j] = B[f][0][j];
				E_inner_before_BCs[f][j] = E[f][0][j];
			}
		}
	}
	
	//--- Apply inner BCs ---
	for( int j=0; j<n_points_t; j++ ){
		B[0][0][j] = B_r_function( r_min, t[j] );
		//B[2][0][j] = 2.0 * delta_T * Omega * st[j] * cos(t[j]);					// SPECIFIC TO A DIPOLE AND NOT SURE IT'S THEORETICALLY GROUNDED; JUST SEEMS TO FIT THE DATA.
		E[1][0][j] = B_r_function( r_min, t[j] ) * - Omega * r_min * st[j];
		E[2][0][j] = 0;
	}
	
	//B[2][0][ (int) n_points_t / 2 - 1 ] = 0;	// Set B_phi ( R_NS, pi/2 ) = 0 to avoid evolution of noise. Suggestion from Sam 20231208.
	
	//--- Record values after inner BCs applied ---
	if( record_before_and_after ){
		for( int j=0; j<n_points_t; j++ ){
			for( int f=0; f<3; f++ ){
				B_inner_after_BCs[f][j] = B[f][0][j];
				E_inner_after_BCs[f][j] = E[f][0][j];
			}
		}
	}

}




void apply_outer_boundary_conditions( bool record_before_and_after = false ){
	// Petri 2012, \S3.3.
	
	//--- Record values before outer BCs applied ---
	if( record_before_and_after ){
		for( int j=0; j<n_points_t; j++ ){
			for( int f=0; f<3; f++ ){
				B_outer_before_BCs[f][j] = B[f].back()[j];
				E_outer_before_BCs[f][j] = E[f].back()[j];
			}
		}
	}
	
	//--- Apply outer BCs ---
	if( T_index < T_index_Omega_ramp_start ){
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
			B[2].back()[j] = 0.5 * ( E[1].back()[j] + B[2].back()[j] );
			E[1].back()[j] =   B[2].back()[j];
			E[2].back()[j] = - B[1].back()[j];
		}
	}
	
	//--- Record values after outer BCs applied ---
	if( record_before_and_after ){
		for( int j=0; j<n_points_t; j++ ){
			for( int f=0; f<3; f++ ){
				B_outer_after_BCs[f][j] = B[f].back()[j];
				E_outer_after_BCs[f][j] = E[f].back()[j];
			}
		}
	}
	
}




void VSH_decomposition_for_E(){
	// Decomposition for vectors which in general are NOT divergence-free. In our code, this is only the E-field, so we may as well make it specific for E.
	// First, built lists of integrand values. Then, numerically integrate and put the values into an array.
	// Two separate loops due to two separate ranges on ell.
	
	for( int i=0; i<n_points_r; i++ ){
		
		//--- r-coefficients go from ell=0 ---
		for( int ell=0; ell<=ell_max; ell++ ){
			
			for( int j=0; j<n_points_t; j++ ){
				VSH_integrand_r[j] = E[0][i][j] * P0[j][ell];
			}
			
			E_VSH_coeffs[0][i][ell] = integral_trapezium_sintheta_dtheta( VSH_integrand_r ) * sqrt_2Lplus1_pi[ell];
			
		}
		
		
		//--- (1)- and (2)-coefficients go from ell=1 ---
		for( int ell=1; ell<=ell_max; ell++ ){
		
			for( int j=0; j<n_points_t; j++ ){
				VSH_integrand_1[j] = E[1][i][j] * P1[j][ell];
				VSH_integrand_2[j] = E[2][i][j] * P1[j][ell];
			}
			
			E_VSH_coeffs[1][i][ell] = integral_trapezium_sintheta_dtheta( VSH_integrand_1 ) * sqrt_2Lplus1_pi_over_L_Lplus1[ell];
			E_VSH_coeffs[2][i][ell] = integral_trapezium_sintheta_dtheta( VSH_integrand_2 ) * sqrt_2Lplus1_pi_over_L_Lplus1[ell];
			
		}
	}
}

/* ADDED 20231206 TO TEST WHY VSH DECOMPOSITION WAS INACCURATE. FOUND ALTERNATIVE CAUSE, BUT DON'T DELETE THIS UNTIL SURE.
double trapezium( std::vector<double> f, std::vector<double> x ) {
	// Trapezium rule for function only of x.
	double ret = f.back()*x.back() - f[0]*x[0];
	
	for( int i=0; i<x.size()-1; i++ ){
		ret += f[i]*x[i+1] - f[i+1]*x[i];
	}
	
	return ret * 0.5;
}
*/




void VSH_decomposition_for_B(){
	// Specific decomposition for divergence-free vectors, where the (1)-coeffs need not be calculated by integration. Instead, they are given in terms of the r-coeffs, but this is only calculate once the radial derivatives are calculated. In our code, this is only the B-field, so we may as well make it specific for B.
	// THIS DOES NOT CALCULATE THE (1)-COEFFICIENTS!
	// First, built lists of integrand values. Then, numerically integrate and put the values into an array.
	// Two separate loops due to two separate ranges on ell.
	
	for( int i=0; i<n_points_r; i++ ){
		
		//--- r-coefficients go from ell=0 ---
		for( int ell=0; ell<=ell_max; ell++ ){
			
			for( int j=0; j<n_points_t; j++ ){
				VSH_integrand_r[j] = B[0][i][j] * P0[j][ell];
				//VSH_integrand_r[j] = B[0][i][j] * P0[j][ell] * st[j];
			}
			
			B_VSH_coeffs[0][i][ell] = integral_trapezium_sintheta_dtheta( VSH_integrand_r ) * sqrt_2Lplus1_pi[ell];
			//B_VSH_coeffs[0][i][ell] = trapezium( VSH_integrand_r, t ) * sqrt_2Lplus1_pi[ell];
			
		}
		
		
		//--- (2)-coefficients go from ell=1 ---
		for( int ell=1; ell<=ell_max; ell++ ){
		
			for( int j=0; j<n_points_t; j++ ){
				VSH_integrand_2[j] = B[2][i][j] * P1[j][ell];
			}
			
			B_VSH_coeffs[2][i][ell] = integral_trapezium_sintheta_dtheta( VSH_integrand_2 ) * sqrt_2Lplus1_pi_over_L_Lplus1[ell];
			
		}
		
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
				B_dr [0][i][j] += sqrt_2Lplus1_over_4pi[ell] * P0[j][ell] * B_VSH_coeffs_dr [0][i][ell];				
				E_dr [0][i][j] += sqrt_2Lplus1_over_4pi[ell] * P0[j][ell] * E_VSH_coeffs_dr [0][i][ell];
				B_dr2[0][i][j] += sqrt_2Lplus1_over_4pi[ell] * P0[j][ell] * B_VSH_coeffs_dr2[0][i][ell];
			}
			
			for( int ell=1; ell<=ell_max; ell++ ){
				B_dr [1][i][j] += sqrt_2Lplus1_over_4pi[ell] * P1[j][ell] * B_VSH_coeffs_dr [1][i][ell];
				B_dr [2][i][j] += sqrt_2Lplus1_over_4pi[ell] * P1[j][ell] * B_VSH_coeffs_dr [2][i][ell];
				E_dr [1][i][j] += sqrt_2Lplus1_over_4pi[ell] * P1[j][ell] * E_VSH_coeffs_dr [1][i][ell];
				E_dr [2][i][j] += sqrt_2Lplus1_over_4pi[ell] * P1[j][ell] * E_VSH_coeffs_dr [2][i][ell];
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




void apply_initial_field_values_for_B(){
	
	for( int i=0; i<n_points_r; i++ ){
		for( int j=0; j<n_points_t; j++ ){
			B[0][i][j] = B_r_function( r[i], t[j] );
			B[1][i][j] = B_t_function( r[i], t[j] );
			B[2][i][j] = B_p_function( r[i], t[j] );
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
	if( ( T_index >= T_index_Omega_ramp_start ) and ( T_index <= T_index_Omega_ramp_stop ) ){
		
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
				E[0][i][j] += ( Omega - Omega_prev ) * pow( r[i], -2 ) * pow( st[j], 2 );
				E[1][i][j] += ( Omega - Omega_prev ) * pow( r[i], -2 ) * -2.0 * st[j] * cos(t[j]);
			}
		}
		
		
		apply_inner_boundary_conditions();
		apply_outer_boundary_conditions();
	}
}




void calculate_n_gridpoints_within_FF_region(){
	
	is_gridpoint_within_FF_region = std::vector< std::vector<int> > ( n_points_r, std::vector<int> ( n_points_t ) );
	n_gridpoints_within_FF_region = 0;
	
	if( Omega > 0 ){
		for( int i=0; i<n_points_r; i++ ){
			for( int j=1; j<n_points_t-1; j++ ){
				if( r[i] * pow( st[j], -2 ) < 1.0 / Omega ){
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
	
	if( Omega > 0 ){
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
			E_div[i][j] += sqrt_2Lplus1_over_4pi[0] * ( E_VSH_coeffs_dr[0][i][0] + 2.0*pow(r[i],-1) * E_VSH_coeffs[0][i][0] );
			
			//- L>0 -
			for( int L=1; L<=ell_max; L++ ){
				
				E_div[i][j] += sqrt_2Lplus1_over_4pi[L] * P0[j][L] * ( E_VSH_coeffs_dr[0][i][L] + 2.0*pow(r[i],-1) * E_VSH_coeffs[0][i][L] - L*(L+1.0)*pow(r[i],-1) * E_VSH_coeffs[1][i][L] );
				
				B_curl[0][i][j] += sqrt_2Lplus1_over_4pi[L] * P0[j][L] * -L*(L+1.0) * B_VSH_coeffs[2][i][L] / r[i];
				B_curl[1][i][j] += sqrt_2Lplus1_over_4pi[L] * P1[j][L] * ( -B_VSH_coeffs_dr[2][i][L] - B_VSH_coeffs[2][i][L] / r[i] );
				B_curl[2][i][j] += sqrt_2Lplus1_over_4pi[L] * P1[j][L] * ( r[i] * B_VSH_coeffs_dr2[0][i][L] + 4.0 * B_VSH_coeffs_dr[0][i][L] - (L+2.0)*(L-1.0) * B_VSH_coeffs[0][i][L] / r[i] ) / ( (double) L*(L+1.0) );
				
				E_curl[0][i][j] += sqrt_2Lplus1_over_4pi[L] * P0[j][L] * -L*(L+1.0) * E_VSH_coeffs[2][i][L] / r[i];
				E_curl[1][i][j] += sqrt_2Lplus1_over_4pi[L] * P1[j][L] * ( -E_VSH_coeffs_dr[2][i][L] - E_VSH_coeffs[2][i][L] / r[i] );
				E_curl[2][i][j] += sqrt_2Lplus1_over_4pi[L] * P1[j][L] * ( E_VSH_coeffs_dr[1][i][L] + ( - E_VSH_coeffs[0][i][L] + E_VSH_coeffs[1][i][L] ) / r[i] );
				
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
	  B_dT[2][0][j] = Omega * ( 2.0*st[j]*B[0][0][j] + ct[j]*B[1][0][j] + st[j]*B_dr[0][0][j] + st[j]*B_t_dt_inner[j] );
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




void calculate_radial_derivatives_of_VSH_coeffs(){
	
	for( int L=0; L<=ell_max; L++ ){
		
		// Right derivatives at inner boundary.
		B_VSH_coeffs_dr [0][0][L] = ( B_VSH_coeffs[0][1][L] - B_VSH_coeffs[0][0][L] ) / delta_r;
		B_VSH_coeffs_dr [2][0][L] = ( B_VSH_coeffs[2][1][L] - B_VSH_coeffs[2][0][L] ) / delta_r;
		E_VSH_coeffs_dr [0][0][L] = ( E_VSH_coeffs[0][1][L] - E_VSH_coeffs[0][0][L] ) / delta_r;
		E_VSH_coeffs_dr [1][0][L] = ( E_VSH_coeffs[1][1][L] - E_VSH_coeffs[1][0][L] ) / delta_r;
		E_VSH_coeffs_dr [2][0][L] = ( E_VSH_coeffs[2][1][L] - E_VSH_coeffs[2][0][L] ) / delta_r;
		
		B_VSH_coeffs_dr2[0][0][L] = ( B_VSH_coeffs[0][2][L] - 2.0*B_VSH_coeffs[0][1][L] + B_VSH_coeffs[0][0][L] ) / ( delta_r*delta_r );
		
		// Symmetric derivatives for all intermediate points.
		for( int i=1; i<n_points_r-1; i++ ){
			
			B_VSH_coeffs_dr [0][i][L] = ( B_VSH_coeffs[0][i+1][L] - B_VSH_coeffs[0][i-1][L] ) / ( 2.0*delta_r );
			B_VSH_coeffs_dr [2][i][L] = ( B_VSH_coeffs[2][i+1][L] - B_VSH_coeffs[2][i-1][L] ) / ( 2.0*delta_r );
			E_VSH_coeffs_dr [0][i][L] = ( E_VSH_coeffs[0][i+1][L] - E_VSH_coeffs[0][i-1][L] ) / ( 2.0*delta_r );
			E_VSH_coeffs_dr [1][i][L] = ( E_VSH_coeffs[1][i+1][L] - E_VSH_coeffs[1][i-1][L] ) / ( 2.0*delta_r );
			E_VSH_coeffs_dr [2][i][L] = ( E_VSH_coeffs[2][i+1][L] - E_VSH_coeffs[2][i-1][L] ) / ( 2.0*delta_r );
					
			B_VSH_coeffs_dr2[0][i][L] = ( B_VSH_coeffs[0][i+1][L] - 2.0*B_VSH_coeffs[0][i][L] + B_VSH_coeffs[0][i-1][L] ) / ( delta_r*delta_r );
		}
		
		// Left derivatives at outer boundary.
		B_VSH_coeffs_dr[0] .back()[L] = ( B_VSH_coeffs[0].back()[L] - B_VSH_coeffs[0][n_points_r-2][L] ) / delta_r;
		B_VSH_coeffs_dr[2] .back()[L] = ( B_VSH_coeffs[2].back()[L] - B_VSH_coeffs[2][n_points_r-2][L] ) / delta_r;
		E_VSH_coeffs_dr[0] .back()[L] = ( E_VSH_coeffs[0].back()[L] - E_VSH_coeffs[0][n_points_r-2][L] ) / delta_r;
		E_VSH_coeffs_dr[1] .back()[L] = ( E_VSH_coeffs[1].back()[L] - E_VSH_coeffs[1][n_points_r-2][L] ) / delta_r;
		E_VSH_coeffs_dr[2] .back()[L] = ( E_VSH_coeffs[2].back()[L] - E_VSH_coeffs[2][n_points_r-2][L] ) / delta_r;
		
		B_VSH_coeffs_dr2[0].back()[L] = ( B_VSH_coeffs[0].back()[L]  - 2.0*B_VSH_coeffs[0][n_points_r-2][L] + B_VSH_coeffs[0][n_points_r-3][L] ) / ( delta_r*delta_r );
		
	}

}





void calculate_B_1_L_and_B_1_L_dr(){
	// This version of the function is for T>0, at which time we want to update B^{(1),L}.
	
	for( int L=1; L<=ell_max; L++ ){
		for( int i=0; i<n_points_r; i++ ){
			B_VSH_coeffs   [1][i][L] = ( r[i] * B_VSH_coeffs_dr [0][i][L] + 2.0 * B_VSH_coeffs   [0][i][L] ) / ( (double) L*(L+1.0) );
			B_VSH_coeffs_dr[1][i][L] = ( r[i] * B_VSH_coeffs_dr2[0][i][L] + 3.0 * B_VSH_coeffs_dr[0][i][L] ) / ( (double) L*(L+1.0) );
		}
	}
	
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




void calculate_stdev(){
	
	//--- Reset values to zero ---
	B_series          = std::vector< std::vector< std::vector<double> > > ( 3, std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) ) );
	E_series          = std::vector< std::vector< std::vector<double> > > ( 3, std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) ) );
	B_abs_rel_dev     = std::vector< std::vector< std::vector<double> > > ( 3, std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) ) );
	E_abs_rel_dev     = std::vector< std::vector< std::vector<double> > > ( 3, std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) ) );
	B_stdev           = std::vector<double> ( 3 );
	E_stdev           = std::vector<double> ( 3 );
	B_max_abs_rel_dev = std::vector<double> ( 3 );
	E_max_abs_rel_dev = std::vector<double> ( 3 );
	
	//--- Calculate stdev ---
	for( int i=0; i<n_points_r; i++ ){
		for( int j=0; j<n_points_t; j++ ){
			
			//--- Evaluate VSH series by summing over ell for each gridpoint ---
			for( int ell=0; ell<=ell_max; ell++ ){
				B_series[0][i][j] += sqrt_2Lplus1_over_4pi[ell] * P0[j][ell] * B_VSH_coeffs[0][i][ell];
				E_series[0][i][j] += sqrt_2Lplus1_over_4pi[ell] * P0[j][ell] * E_VSH_coeffs[0][i][ell];
			}
			
			for( int ell=1; ell<=ell_max; ell++ ){
				B_series[1][i][j] += sqrt_2Lplus1_over_4pi[ell] * P1[j][ell] * B_VSH_coeffs[1][i][ell];
				B_series[2][i][j] += sqrt_2Lplus1_over_4pi[ell] * P1[j][ell] * B_VSH_coeffs[2][i][ell];
				E_series[1][i][j] += sqrt_2Lplus1_over_4pi[ell] * P1[j][ell] * E_VSH_coeffs[1][i][ell];
				E_series[2][i][j] += sqrt_2Lplus1_over_4pi[ell] * P1[j][ell] * E_VSH_coeffs[2][i][ell];
			}
			
			//--- Calculate stdev and build lists of square deviations along the way ---
			for( int f=0; f<3; f++ ){
				
				B_stdev[f] += pow( B[f][i][j] - B_series[f][i][j], 2 );
				E_stdev[f] += pow( E[f][i][j] - E_series[f][i][j], 2 );
				
				if( B[f][i][j] == 0 ){
					B_abs_rel_dev[f][i][j] = 0;
				}
				else {
					B_abs_rel_dev[f][i][j] = abs( 1.0 - B_series[f][i][j] / B[f][i][j] );
				}
				
				if( E[f][i][j] == 0 ){
					E_abs_rel_dev[f][i][j] = 0;
				}
				else {
					E_abs_rel_dev[f][i][j] = abs( 1.0 - E_series[f][i][j] / E[f][i][j] );
				}
				
				if( B_abs_rel_dev[f][i][j] > B_max_abs_rel_dev[f] ){
					B_max_abs_rel_dev[f] = B_abs_rel_dev[f][i][j];
				}
				
				if( E_abs_rel_dev[f][i][j] > E_max_abs_rel_dev[f] ){
					E_max_abs_rel_dev[f] = E_abs_rel_dev[f][i][j];
				}
				
			}
			
		}
	}
	
	for( int f=0; f<3; f++ ){
		B_stdev[f] = sqrt( B_stdev[f] / ( (double) n_points_total ) );
		E_stdev[f] = sqrt( E_stdev[f] / ( (double) n_points_total ) );
	}
	
}




void calculate_magnetic_energy(){
	
	for( int i=0; i<n_points_r; i++ ){
		for( int j=0; j<n_points_t; j++ ){
			magnetic_energy_density[i][j] = pow( B[0][i][j], 2 ) + pow( B[1][i][j], 2 ) + pow( B[2][i][j], 2 );
		}
	}
	
	total_magnetic_energy = integral_2d( magnetic_energy_density, integral_r_t_dr, integral_r_t_dt ) * 2.0 * pi;
	
}




double estimate_csv_profiles_filesize(){
	
	double fudge_factor = 17.0;		// Free parameter to tune final result.
	int n_columns = 19 + 8 * 6;		// You are responsible for keeping this value updated based on how many values are specified in void output_headers_to_csv_profiles().
	int n_rows = 0;					// Calculated by the for-loop below.
	
	for( int T_index_here=0; T_index_here<=n_timesteps; T_index_here++ ){
		bool condition_T_1 = ( T_index_here >= csv_profiles_write_T_min ) and ( T_index_here <= csv_profiles_write_T_max );
		bool condition_T_2 = ( T_index_here % csv_profiles_write_freq_T == 0 ) or ( T_index_here == n_timesteps );
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
	bool condition_T_2 = ( T_index % csv_profiles_write_freq_T == 0 ) or ( T_index == n_timesteps );
	
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
					
					
					output_file_profiles <<std::setprecision(csv_precision)<< T_index <<","<< T <<","<< T*T_factor <<","<< i <<","<< j <<","<< r[i] <<","<< t[j] <<","<< x[i][j] <<","<< z[i][j]
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
	
	output_file_history << "T_index,T,T_sec";
	output_file_history << ",Omega,R_LC,star_rotation_angle";
	output_file_history << ",E_pts_chngd,pct_E_pts_chngd";
	output_file_history << ",stdev_B_r,stdev_E_r,stdev_B_t,stdev_E_t,stdev_B_p,stdev_E_p";
	output_file_history << ",max_abs_rel_dev_B_r,max_abs_rel_dev_E_r,max_abs_rel_dev_B_t,max_abs_rel_dev_E_t,max_abs_rel_dev_B_p,max_abs_rel_dev_E_p";
	output_file_history << ",nans_B_r,nans_E_r,nans_B_t,nans_E_t,nans_B_p,nans_E_p";
	output_file_history << ",total_magnetic_energy";
	output_file_history << "\n";
}




void output_to_csv_history(){
	
	bool condition_T = ( T_index % csv_history_write_freq_T == 0 ) or ( T_index == n_timesteps );
	
	if( condition_T ){
		
		double R_LC = 0;
		if( Omega > 0 ){
			R_LC = 1.0 / Omega;
		}
		
		output_file_history <<std::setprecision(csv_precision)<< T_index <<","<< T <<","<< T*T_factor
		                                                      <<","<< Omega <<","<< R_LC <<","<< star_rotation_angle
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
	
	bool condition_T = ( T_index % csv_VSH_coeffs_write_freq_T == 0 ) or ( T_index == n_timesteps );
	
	if( condition_T ){
		for( int i=0; i<n_points_r; i++ ){
			for( int ell=0; ell<=ell_max; ell++ ){
				
				output_file_VSH_coeffs <<std::setprecision(csv_precision)<< T_index <<","<< i <<","<< r[i] <<","<< ell;
				
				for( int f=0; f<3; f++ ){ output_file_VSH_coeffs <<std::setprecision(csv_precision)<<","<< B_VSH_coeffs   [f][i][ell] <<std::setprecision(csv_precision)<<","<< E_VSH_coeffs   [f][i][ell]; }
				for( int f=0; f<3; f++ ){ output_file_VSH_coeffs <<std::setprecision(csv_precision)<<","<< B_VSH_coeffs_dr[f][i][ell] <<std::setprecision(csv_precision)<<","<< E_VSH_coeffs_dr[f][i][ell]; }
				
				output_file_VSH_coeffs <<std::setprecision(csv_precision)<< ","<< B_VSH_coeffs_dr2[0][i][ell] << "\n";
				
			}
		}
	}
	
}




void output_headers_to_csv_BCs(){
	// Split over separate lines for readability and ease of checking consistency with output_to_csv_profiles().
	
	// Remember to update int n_columns in void estimate_csv_filesize().
	
	output_file_BCs << "T_index,T,T_sec,j,t";
	output_file_BCs << ",B_r_inner_before_BCs,E_r_inner_before_BCs,B_t_inner_before_BCs,E_t_inner_before_BCs,B_p_inner_before_BCs,E_p_inner_before_BCs";
	output_file_BCs << ",B_r_inner_after_BCs,E_r_inner_after_BCs,B_t_inner_after_BCs,E_t_inner_after_BCs,B_p_inner_after_BCs,E_p_inner_after_BCs";
	output_file_BCs << ",rel_diff_B_r_inner,rel_diff_E_r_inner,rel_diff_B_t_inner,rel_diff_E_t_inner,rel_diff_B_p_inner,rel_diff_E_p_inner";
	output_file_BCs << ",B_r_outer_before_BCs,E_r_outer_before_BCs,B_t_outer_before_BCs,E_t_outer_before_BCs,B_p_outer_before_BCs,E_p_outer_before_BCs";
	output_file_BCs << ",B_r_outer_after_BCs,E_r_outer_after_BCs,B_t_outer_after_BCs,E_t_outer_after_BCs,B_p_outer_after_BCs,E_p_outer_after_BCs";
	output_file_BCs << ",rel_diff_B_r_outer,rel_diff_E_r_outer,rel_diff_B_t_outer,rel_diff_E_t_outer,rel_diff_B_p_outer,rel_diff_E_p_outer";
	output_file_BCs << "\n";
	
}




void output_to_csv_BCs(){
	
	bool condition_T = ( T_index % csv_BCs_write_freq_T == 0 ) or ( T_index == n_timesteps );
	
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
		
			output_file_BCs <<std::setprecision(csv_precision)<< T_index <<","<< T <<","<< T*T_factor <<","<< j <<","<< t[j];
			
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
	
	if( ( T_index % cout_freq_T == 0 ) or ( T_index == n_timesteps ) ){
	
		double percent_done = round( (double)1000*T_index/n_timesteps )/10; // This format makes it exactly one decimal place.
		
		std::cout <<std::left<<std::setw(w)<< T_index <<std::left<<std::setw(w)<< T <<"|\t"
				  <<std::left<<std::setw(w)<< B[0][cout_i][cout_j] <<std::left<<std::setw(w)<< B[1][cout_i][cout_j] <<std::left<<std::setw(w)<< B[2][cout_i][cout_j] <<"|\t"
				  <<std::left<<std::setw(w)<< E[0][cout_i][cout_j] <<std::left<<std::setw(w)<< E[1][cout_i][cout_j] <<std::left<<std::setw(w)<< E[2][cout_i][cout_j] <<"|\t"
				  <<std::left<<std::setw(w)<< percentage_E_points_points_changed_for_second_FF_condition <<"|\t"
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




void output_parameters_to_screen(){
	// No need to estimate the size of the history CSV since it's far smaller than the profiles CSV.
	// For the future, perhaps split into multiple profiles CSVs. Then, add them all together to get an estimated total space required. A nice feature would be to fail if that much space isn't available, or if it exceeds some specified threshold.
	// NOTE: Number of rotation periods needs updating to include the ramp. Calculate on paper by integrating dOmega_by_dT.
	// Or add a variable rotation_angle which just gets Omega * delta_T added to it at each timestep, and can be recorded in the _history.csv file.
	
	std::cout << log_file_comment << "\n" << std::endl;
	
	std::cout << "Time start:\t" << time_start_string << std::endl;
	
	std::cout << "Parameter values (given in code units, then in SI units)" << std::endl;
	std::cout << "Omega_max              :\t" <<std::left<<std::setw(16)<< Omega_max<< std::endl;
	std::cout << "R_LC_max               :\t" <<std::left<<std::setw(16)<< 1.0 / Omega_max    << std::endl;
	std::cout << "r_min, r_max           :\t" <<std::left<<std::setw(16)<< r_min   <<std::left<<std::setw(16)<< r_max              << std::endl;
	std::cout << "n_points_r, delta_r    :\t" <<std::left<<std::setw(16)<< n_points_r <<std::left<<std::setw(16)<< delta_r << std::endl;
	std::cout << "n_points_t, delta_t    :\t" <<std::left<<std::setw(16)<< n_points_t <<std::left<<std::setw(16)<< delta_t << std::endl;
	std::cout << "Timestep               :\t" <<std::left<<std::setw(16)<< delta_T <<std::left<<std::setw(16)<< delta_T * T_factor << std::endl;
	std::cout << "CFL max timestep       :\t" <<std::left<<std::setw(16)<< delta_r <<std::left<<std::setw(16)<< delta_r * T_factor << std::endl;
	std::cout << "Length of simulation   :\t" <<std::left<<std::setw(16)<< T_max << T_max * T_factor <<"\tsteps: " << n_timesteps << std::endl;
	std::cout << "Est. CSV size (MB)     :\t" << estimate_csv_profiles_filesize() << std::endl;
	std::cout << "ell_max                :\t" << ell_max << std::endl;
	std::cout << "Printing values at ( r[" << cout_i <<"], theta[" << cout_j <<"] ) = ( " << r[cout_i] <<", " << t[cout_j] << " ) = ( " << r[cout_i] <<", " << t[cout_j]/pi << " pi )." << std::endl;
	std::cout << "dOmega_by_dT_index     :\t" << dOmega_by_dT_index <<"\tdOmega_by_dT:\t"<< dOmega_by_dT << std::endl;
	std::cout << "Ramp start             :\t" << T_index_Omega_ramp_start <<"\t"<< T_index_Omega_ramp_start * delta_T << std::endl;
	std::cout << "Ramp stop              :\t" << T_index_Omega_ramp_stop  <<"\t"<< T_index_Omega_ramp_stop  * delta_T << std::endl;
	
	std::cout << "\nuse_outer_sponge_layer:\t" << use_outer_sponge_layer << std::endl;
	std::cout << "sigma_0, gamma, beta    :\t" << friction_sigma_0 <<"\t"<< friction_gamma <<"\t"<< friction_beta << std::endl;
	
	std::cout << "\nCSV file saved:\tCSV/"  << output_filename << "_1_profiles.csv"   << std::endl;
	std::cout <<   "CSV file saved:\tCSV/"  << output_filename << "_2_history.csv"    << std::endl;
	std::cout <<   "CSV file saved:\tCSV/"  << output_filename << "_3_VSH_coeffs.csv" << std::endl;
	std::cout <<   "CSV file saved:\tCSV/"  << output_filename << "_4_BCs.csv"        << std::endl;
	std::cout <<   "Log file saved:\tLogs/" << output_filename << ".txt"              << std::endl;
	
	std::cout << "\nTo cut evolution short, press CTRL+C at any time. This is perfectly safe and the CSV file will still be saved up to that point.\n" << std::endl;
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
	output_file_log << "ell_max                :\t" << ell_max << "\n";
	output_file_log << "Printing values at ( r[" << cout_i <<"], theta[" << cout_j <<"] ) = ( " << r[cout_i] <<", " << t[cout_j] << " ) = ( " << r[cout_i] <<", " << t[cout_j]/pi << " pi )." << "\n";
	output_file_log << "dOmega_by_dT_index     :\t" << dOmega_by_dT_index <<"\tdOmega_by_dT:\t"<< dOmega_by_dT << "\n";
	output_file_log << "Ramp start             :\t" << T_index_Omega_ramp_start <<"\t"<< T_index_Omega_ramp_start * delta_T << "\n";
	output_file_log << "Ramp stop              :\t" << T_index_Omega_ramp_stop  <<"\t"<< T_index_Omega_ramp_stop  * delta_T << "\n";
	
	output_file_log << "\nuse_outer_sponge_layer:\t" << use_outer_sponge_layer << "\n";
	output_file_log << "sigma_0, gamma, beta    :\t" << friction_sigma_0 <<"\t"<< friction_gamma <<"\t"<< friction_beta << "\n";
	
	output_file_log << "\nCSV file saved:\tCSV/" << output_filename << "_1_profiles.csv"   << "\n";
	output_file_log <<   "CSV file saved:\tCSV/" << output_filename << "_2_history.csv"    << "\n";
	output_file_log <<   "CSV file saved:\tCSV/" << output_filename << "_3_VSH_coeffs.csv" << "\n";
	output_file_log <<   "CSV file saved:\tCSV/" << output_filename << "_4_BCs.csv"        << "\n";
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
