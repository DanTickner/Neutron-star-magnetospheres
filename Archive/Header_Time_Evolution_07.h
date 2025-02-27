/*
Header_Time_Evolution_07.h

Code snippets for the Time_Evolution_xx.cpp codes. Keeps them shorter and allows for two separate for-loops with Euler and Adams-Bashforth integration
without copying code.

V07: Force-free condition is not valid everywhere within the light cylinder; it must also be within the last closed field line.
Add this extra condition and change variable names to reflect this.
n_gridpoints_inside_light_cylinder -> n_gridpoints_within_force_free_region
void calculate_n_gridpoints_inside_light_cylinder() -> void calculate_n_gridpoints_within_force_free_region()
void apply_second_force_free_condition_within_light_cylinder() -> void apply_second_force_free_condition_within_force_free_region()
As a result of this other condition, we can't have theta=0, so adjust the for-loop to account for this (start at j=1).
Of course, this was true for the original condition, but didn't represent a possible division by zero error so had been left in.

NOTICED ERROR IN EXPRESSION FOR E_cross_B, which has been corrected in this and all future header files.
*/

//----- Global variables (don't touch) -----
double delta_r = ( r_max - r_min ) / ( (double) n_points_r - 1.0 );
double delta_t = pi / ( (double) n_points_t - 1.0 );

std::vector<double> r  ( n_points_r );
std::vector<double> t  ( n_points_t );
std::vector<double> st ( n_points_t );
std::vector<double> ct ( n_points_t );

std::vector< std::vector<double> > x ( n_points_r, std::vector<double> ( n_points_t ) );
std::vector< std::vector<double> > z ( n_points_r, std::vector<double> ( n_points_t ) );

std::vector< std::vector<double> > P0 ( n_points_t, std::vector<double> ( ell_max+1 ) );
std::vector< std::vector<double> > P1 ( n_points_t, std::vector<double> ( ell_max+1 ) );

std::vector<double> sqrt_2Lplus1_over_4pi ( ell_max+1 );

std::vector< std::vector<double> > B_r_L    ( n_points_r, std::vector<double> ( ell_max+1 ) );
std::vector< std::vector<double> > B_1_L    ( n_points_r, std::vector<double> ( ell_max+1 ) );
std::vector< std::vector<double> > B_2_L    ( n_points_r, std::vector<double> ( ell_max+1 ) );
std::vector< std::vector<double> > E_r_L    ( n_points_r, std::vector<double> ( ell_max+1 ) );
std::vector< std::vector<double> > E_1_L    ( n_points_r, std::vector<double> ( ell_max+1 ) );
std::vector< std::vector<double> > E_2_L    ( n_points_r, std::vector<double> ( ell_max+1 ) );
std::vector< std::vector<double> > B_r_L_dr ( n_points_r, std::vector<double> ( ell_max+1 ) );
std::vector< std::vector<double> > B_1_L_dr ( n_points_r, std::vector<double> ( ell_max+1 ) );
std::vector< std::vector<double> > B_2_L_dr ( n_points_r, std::vector<double> ( ell_max+1 ) );
std::vector< std::vector<double> > E_r_L_dr ( n_points_r, std::vector<double> ( ell_max+1 ) );
std::vector< std::vector<double> > E_1_L_dr ( n_points_r, std::vector<double> ( ell_max+1 ) );
std::vector< std::vector<double> > E_2_L_dr ( n_points_r, std::vector<double> ( ell_max+1 ) );
std::vector< std::vector<double> > B_r_L_dr2( n_points_r, std::vector<double> ( ell_max+1 ) );

std::vector< std::vector<double> > B_r          ( n_points_r, std::vector<double> ( n_points_t ) );
std::vector< std::vector<double> > B_t          ( n_points_r, std::vector<double> ( n_points_t ) );
std::vector< std::vector<double> > B_p          ( n_points_r, std::vector<double> ( n_points_t ) );
std::vector< std::vector<double> > E_r          ( n_points_r, std::vector<double> ( n_points_t ) );
std::vector< std::vector<double> > E_t          ( n_points_r, std::vector<double> ( n_points_t ) );
std::vector< std::vector<double> > E_p          ( n_points_r, std::vector<double> ( n_points_t ) );
std::vector< std::vector<double> > curl_B_r     ( n_points_r, std::vector<double> ( n_points_t ) );		// Can assign to single variable within void calculate_time_derivatives() instead of global list to save memory, but want to output to CSV to check values 20230914.
std::vector< std::vector<double> > curl_B_t     ( n_points_r, std::vector<double> ( n_points_t ) );
std::vector< std::vector<double> > curl_B_p     ( n_points_r, std::vector<double> ( n_points_t ) );
std::vector< std::vector<double> > curl_E_r     ( n_points_r, std::vector<double> ( n_points_t ) );
std::vector< std::vector<double> > curl_E_t     ( n_points_r, std::vector<double> ( n_points_t ) );
std::vector< std::vector<double> > curl_E_p     ( n_points_r, std::vector<double> ( n_points_t ) );
std::vector< std::vector<double> > div_E        ( n_points_r, std::vector<double> ( n_points_t ) );

std::vector< std::vector<double> > B_r_dT       ( n_points_r, std::vector<double> ( n_points_t ) );
std::vector< std::vector<double> > B_t_dT       ( n_points_r, std::vector<double> ( n_points_t ) );
std::vector< std::vector<double> > B_p_dT       ( n_points_r, std::vector<double> ( n_points_t ) );
std::vector< std::vector<double> > E_r_dT       ( n_points_r, std::vector<double> ( n_points_t ) );
std::vector< std::vector<double> > E_t_dT       ( n_points_r, std::vector<double> ( n_points_t ) );
std::vector< std::vector<double> > E_p_dT       ( n_points_r, std::vector<double> ( n_points_t ) );
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
std::vector< std::vector<double> > sqdev_B_r    ( n_points_r, std::vector<double> ( n_points_t ) );		// To check where the VSH decomposition is in/accurate.
std::vector< std::vector<double> > sqdev_B_t    ( n_points_r, std::vector<double> ( n_points_t ) );
std::vector< std::vector<double> > sqdev_B_p    ( n_points_r, std::vector<double> ( n_points_t ) );
std::vector< std::vector<double> > sqdev_E_r    ( n_points_r, std::vector<double> ( n_points_t ) );
std::vector< std::vector<double> > sqdev_E_t    ( n_points_r, std::vector<double> ( n_points_t ) );
std::vector< std::vector<double> > sqdev_E_p    ( n_points_r, std::vector<double> ( n_points_t ) );
std::vector< std::vector<double> > E_dot_B      ( n_points_r, std::vector<double> ( n_points_t ) );		// To check that force-free conditions are fulfilled.
std::vector< std::vector<double> > Bsq_minus_Esq( n_points_r, std::vector<double> ( n_points_t ) );		// To check that force-free conditions are fulfilled.
std::vector< std::vector<double> > alpha        ( n_points_r, std::vector<double> ( n_points_t ) );		// Don't have to be 2D array but might want to track in CSV file.
std::vector< std::vector<double> > beta         ( n_points_r, std::vector<double> ( n_points_t ) );		// Don't have to be 2D array but might want to track in CSV file.
std::vector< std::vector<double> > magnetic_energy_density ( n_points_r, std::vector<double> ( n_points_t ) );		// Energy density


std::vector<double> integral_r_t_dr ( n_points_r );	// Radial  part of finite element for integrals over r^2 sin(theta) dr dtheta.
std::vector<double> integral_r_t_dt ( n_points_t );	// Angular part of finite element for integrals over r^2 sin(theta) dr dtheta.

double stdev_B_r = 0;
double stdev_B_t = 0;
double stdev_B_p = 0;
double stdev_E_r = 0;
double stdev_E_t = 0;
double stdev_E_p = 0;
double max_sqdev_B_r = 0;
double max_sqdev_B_t = 0;
double max_sqdev_B_p = 0;
double max_sqdev_E_r = 0;
double max_sqdev_E_t = 0;
double max_sqdev_E_p = 0;
double total_magnetic_energy = 0;

std::vector<double> VSH_integrand_r ( n_points_t );
std::vector<double> VSH_integrand_1 ( n_points_t );
std::vector<double> VSH_integrand_2 ( n_points_t );

std::ofstream output_file_profiles;
std::ofstream output_file_history;

double time_start_seconds = 0;

double T = 0;

int n_gridpoints_within_force_free_region = 0;
double percentage_E_points_points_changed_for_second_force_free_condition = 0;	// Output to CSV to keep track of how active the FF condition application is.




//----- Functions -----


void calculate_gridpoints(){
	
	for( int i=0; i<n_points_r; i++ ){
		r[i] = r_min + i * delta_r;
	}
	for( int j=0; j<n_points_t; j++ ){
		t [j] = j * delta_t;
		st[j] = sin( t[j] );
		ct[j] = cos( t[j] );
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
		
		P0[j][0] = 1.0;
		P0[j][1] = x;
		P0[j][2] = 0.5 * ( 3.0*x*x - 1.0 );		// Hardcode so that we can use the same for-loop for P0 and P1, starting at ell=3.
		
		P1[j][1] = sqrt( 1.0 - x*x ) * -1.0;
		P1[j][2] = sqrt( 1.0 - x*x ) * -3.0 * x;
		
		for( int ell=3; ell<=ell_max; ell++ ){
			P0[j][ell] = ( P0[j][ell-1] * x * ( 2.0*ell-1.0 ) - P0[j][ell-2] * ( ell-1.0 ) ) / ( (double)ell );
			P1[j][ell] = ( P1[j][ell-1] * x * ( 2.0*ell-1.0 ) - P1[j][ell-2] * (double)ell ) / ( ell-1.0 );
		}
		
	}
	
}




void calculate_sqrt_2Lplus1_over_4pi(){
	for( int ell=0; ell<=ell_max; ell++ ){
		sqrt_2Lplus1_over_4pi[ell] = sqrt( (2.0*ell+1.0) / ( 4.0*pi ) );
	}
}




double integral_trapezium( std::vector<double> f, std::vector<double> x ) {
	// Trapezium rule for function only of x.
	double ret = f.back()*x.back() - f[0]*x[0];
	
	for( int i=0; i<x.size()-1; i++ ){
		ret += f[i]*x[i+1] - f[i+1]*x[i];
	}
	
	return ret * 0.5;
}




double integral_2d( std::vector< std::vector<double> > f, std::vector<double> du1, std::vector<double> du2 ){
	// Integral of a 2D function given as a list, with the area elements already calculated.
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




void VSH_decomposition_for_E(){
	for( int i=0; i<n_points_r; i++ ){
		for( int ell=0; ell<=ell_max; ell++ ){
		
			//--- Build lists of integrand values ---
			for( int j=0; j<n_points_t; j++ ){
				VSH_integrand_r[j] = E_r[i][j] * P0[j][ell] * st[j];
				VSH_integrand_1[j] = E_t[i][j] * P1[j][ell] * st[j];
				VSH_integrand_2[j] = E_p[i][j] * P1[j][ell] * st[j];
			}
			
			//--- Integrate and put values into array ---
			E_r_L[i][ell] = integral_trapezium( VSH_integrand_r, t ) * sqrt( (2.0*ell+1.0)*pi );
			E_1_L[i][ell] = integral_trapezium( VSH_integrand_1, t ) * sqrt( (2.0*ell+1.0)*pi ) / ( (double)ell*(ell+1.0) );
			E_2_L[i][ell] = integral_trapezium( VSH_integrand_2, t ) * sqrt( (2.0*ell+1.0)*pi ) / ( (double)ell*(ell+1.0) );
			
			
			if( ell == 0 ){
				E_1_L[i][ell] = 0;
				E_2_L[i][ell] = 0;
			}
			
		}
	}
}

void VSH_decomposition_for_B(){
	for( int i=0; i<n_points_r; i++ ){
		for( int ell=0; ell<=ell_max; ell++ ){
		
			//--- Build lists of integrand values ---
			for( int j=0; j<n_points_t; j++ ){
				VSH_integrand_r[j] = B_r[i][j] * P0[j][ell] * st[j];
				VSH_integrand_1[j] = B_t[i][j] * P1[j][ell] * st[j];
				VSH_integrand_2[j] = B_p[i][j] * P1[j][ell] * st[j];
			}
			
			//--- Integrate and put values into array ---
			B_r_L[i][ell] = integral_trapezium( VSH_integrand_r, t ) * sqrt( (2.0*ell+1.0)*pi );
			B_1_L[i][ell] = integral_trapezium( VSH_integrand_1, t ) * sqrt( (2.0*ell+1.0)*pi ) / ( (double)ell*(ell+1.0) );
			B_2_L[i][ell] = integral_trapezium( VSH_integrand_2, t ) * sqrt( (2.0*ell+1.0)*pi ) / ( (double)ell*(ell+1.0) );
			
			
			if( ell == 0 ){
				B_1_L[i][ell] = 0;
				B_2_L[i][ell] = 0;
			}
			
		}
	}
}




void evaluate_VSH_series(){
	
	//--- Reset values to zero ---
	B_r = std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) );
	B_t = std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) );
	B_p = std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) );
	E_r = std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) );
	E_t = std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) );
	E_p = std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) );
	
	//--- Sum over ell for each gridpoint ---
	for( int i=0; i<n_points_r; i++ ){
		for( int j=0; j<n_points_t; j++ ){
			for( int ell=0; ell<=ell_max; ell++ ){
				
				B_r[i][j] += sqrt_2Lplus1_over_4pi[i] * P0[j][ell] * B_r_L[i][ell];
				B_t[i][j] += sqrt_2Lplus1_over_4pi[i] * P1[j][ell] * B_1_L[i][ell];
				B_t[i][j] += sqrt_2Lplus1_over_4pi[i] * P1[j][ell] * B_2_L[i][ell];
				
				E_r[i][j] += sqrt_2Lplus1_over_4pi[i] * P0[j][ell] * E_r_L[i][ell];
				E_t[i][j] += sqrt_2Lplus1_over_4pi[i] * P1[j][ell] * E_1_L[i][ell];
				E_t[i][j] += sqrt_2Lplus1_over_4pi[i] * P1[j][ell] * E_2_L[i][ell];
				
			}
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
			for( int ell=0; ell<=ell_max; ell++ ){
		
				double sqrt_factor = sqrt( (2.0*ell+1.0) / ( 4.0 * pi ) );
				
				B_r_dr[i][j] += sqrt_factor * P0[j][ell] * B_r_L_dr[i][ell];
				B_t_dr[i][j] += sqrt_factor * P1[j][ell] * B_1_L_dr[i][ell];
				B_t_dr[i][j] += sqrt_factor * P1[j][ell] * B_2_L_dr[i][ell];
				
				E_r_dr[i][j] += sqrt_factor * P0[j][ell] * E_r_L_dr[i][ell];
				E_t_dr[i][j] += sqrt_factor * P1[j][ell] * E_1_L_dr[i][ell];
				E_t_dr[i][j] += sqrt_factor * P1[j][ell] * E_2_L_dr[i][ell];
				
			}
		}
	}
}




void calculate_Cartesian_vector_components(){
	
	for( int i=0; i<n_points_r; i++ ){
		for( int j=0; j<n_points_t; j++ ){
			
			B_x[i][j] = B_r[i][j] * st[j] + B_t[i][j] * ct[j];
			B_y[i][j] = B_p[i][j];									// Extra variable assignment not necessary, but use for now to avoid confusion.
			B_z[i][j] = B_r[i][j] * ct[j] - B_t[i][j] * st[j];
			
			E_x[i][j] = E_r[i][j] * st[j] + E_t[i][j] * ct[j];
			E_y[i][j] = E_p[i][j];									// Extra variable assignment not necessary, but use for now to avoid confusion.
			E_z[i][j] = E_r[i][j] * ct[j] - E_t[i][j] * st[j];
			
		}
	}
}




void apply_initial_field_value_and_VSH_decomposition_for_B(){
	
	for( int i=0; i<n_points_r; i++ ){
		
		for( int j=0; j<n_points_t; j++ ){
			B_r[i][j] = B_r_function( r[i], t[j] );
			B_t[i][j] = B_t_function( r[i], t[j] );
			B_p[i][j] = B_p_function( r[i], t[j] );
		}
		
		for( int L=0; L<=ell_max; L++ ){
			B_r_L[i][L] = B_r_L_function( r[i], L );
			B_1_L[i][L] = B_1_L_function( r[i], L );
			B_2_L[i][L] = B_2_L_function( r[i], L );
		}
		
	}
}
void apply_initial_field_value_and_VSH_decomposition_for_E(){
	
	for( int i=0; i<n_points_r; i++ ){
		
		for( int j=0; j<n_points_t; j++ ){
			E_r[i][j] = E_r_function( r[i], t[j] );
			E_t[i][j] = E_t_function( r[i], t[j] );
			E_p[i][j] = E_p_function( r[i], t[j] );
		}
		
		for( int L=0; L<=ell_max; L++ ){
			E_r_L[i][L] = E_r_L_function( r[i], L );
			E_1_L[i][L] = E_1_L_function( r[i], L );
			E_2_L[i][L] = E_2_L_function( r[i], L );
		}
		
	}
}




void calculate_n_gridpoints_within_force_free_region(){
	for( int i=0; i<n_points_r; i++ ){
		for( int j=1; j<n_points_t; j++ ){
			if( r[i] * pow( st[j], -2 ) < R_LC ){
				n_gridpoints_within_force_free_region++;
			}
		}
	}
	std::cout << "Total number of gridpoints         :\t" << n_points_r * n_points_t << std::endl;
	std::cout << "Gridpoints within force-free region:\t" << n_gridpoints_within_force_free_region << "\t(" << 100 * n_gridpoints_within_force_free_region / ( (double) n_points_r * n_points_t ) << " %)" << std::endl;
}



void apply_second_force_free_condition_within_force_free_region(){
	// Petri 2012, Section 3.2.
	
	//--- Update the electric field at all points within the light cylinder if the second FF condition B^2-E^2<0 is violated ---
	int num_E_points_changed_for_second_force_free_condition = 0;
	
	for( int i=0; i<n_points_r; i++ ){
		for( int j=1; j<n_points_t; j++ ){
			if( r[i] * pow( st[j], -2 ) < R_LC ){
				
				double Bsq = B_r[i][j]*B_r[i][j] + B_t[i][j]*B_t[i][j] + B_p[i][j]*B_p[i][j];
				double Esq = E_r[i][j]*E_r[i][j] + E_t[i][j]*E_t[i][j] + E_p[i][j]*E_p[i][j];
				if( Bsq - Esq < 0 ){
					
					num_E_points_changed_for_second_force_free_condition += 1;
					
					double E_dot_B_over_Bsq = ( E_r[i][j]*B_r[i][j]+E_t[i][j]*B_t[i][j]+E_p[i][j]*B_p[i][j] ) / Bsq;
					double E_prime_r = E_r[i][j] - E_dot_B_over_Bsq * B_r[i][j];
					double E_prime_t = E_t[i][j] - E_dot_B_over_Bsq * B_t[i][j];
					double E_prime_p = E_p[i][j] - E_dot_B_over_Bsq * B_p[i][j];
					double sqrt_Bsq_over_Eprimesq = sqrt( Bsq / ( E_prime_r*E_prime_r + E_prime_t*E_prime_t + E_prime_p*E_prime_p ) );
					
					E_r[i][j] = E_prime_r * sqrt_Bsq_over_Eprimesq;
					E_t[i][j] = E_prime_t * sqrt_Bsq_over_Eprimesq;
					E_p[i][j] = E_prime_p * sqrt_Bsq_over_Eprimesq;
					
				}
			}
		}
	}
	
	// Calculate percentage of gridpoints at which E changed. Below format makes it exactly one decimal place.
	percentage_E_points_points_changed_for_second_force_free_condition = round( (double)1000 * num_E_points_changed_for_second_force_free_condition / n_gridpoints_within_force_free_region ) / 10;
}




void calculate_time_derivatives(){
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
			for( int L=1; L<=ell_max; L++ ){
				
				div_E[i][j] += sqrt_2Lplus1_over_4pi[L] * ( E_r_L_dr[i][L] + 2.0*pow(r[i],-1) * E_r_L[i][L] - L*(L+1.0)*pow(r[i],-1) * E_1_L[i][L] );
				
				curl_E_r[i][j] += sqrt_2Lplus1_over_4pi[L] * P0[j][L] * -L*(L+1) * E_2_L[i][L] / r[i];
				curl_E_t[i][j] += sqrt_2Lplus1_over_4pi[L] * P1[j][L] * ( -E_2_L_dr[i][L] - E_2_L[i][L] / r[i] );
				curl_E_p[i][j] += sqrt_2Lplus1_over_4pi[L] * P1[j][L] * ( E_1_L_dr[i][L] + ( E_1_L[i][L] - E_r_L[i][L] ) / r[i] );
				
				curl_B_r[i][j] += sqrt_2Lplus1_over_4pi[L] * P0[j][L] * -L*(L+1) * B_2_L[i][L] / r[i];
				curl_B_t[i][j] += sqrt_2Lplus1_over_4pi[L] * P1[j][L] * ( -B_2_L_dr[i][L] - B_2_L[i][L] / r[i] );
				curl_B_p[i][j] += sqrt_2Lplus1_over_4pi[L] * P1[j][L] * ( r[i] * B_r_L_dr2[i][L] + 4.0 * B_r_L_dr[i][L] - (L+2.0)*(L-1.0) * B_r_L[i][L] / r[i] ) / ( L*(L+1.0) );
				
			}
			
			//--- alpha and beta ---
			double B_dot_curl_B = B_r[i][j]*curl_B_r[i][j] + B_t[i][j]*curl_B_t[i][j] + B_p[i][j]*curl_B_p[i][j];
			double E_dot_curl_E = E_r[i][j]*curl_E_r[i][j] + E_t[i][j]*curl_E_t[i][j] + E_p[i][j]*curl_E_p[i][j];
			double B_squared    = B_r[i][j]*B_r     [i][j] + B_t[i][j]*     B_t[i][j] + B_p[i][j]*     B_p[i][j];
			double E_squared    = E_r[i][j]*E_r     [i][j] + E_t[i][j]*     E_t[i][j] + E_p[i][j]*     E_p[i][j];
			
			alpha[i][j] = ( E_dot_curl_E - B_dot_curl_B ) / B_squared;
			beta [i][j]  = - div_E[i][j] / B_squared;
			
			//--- E cross B ---
			double E_cross_B_r = E_t[i][j] * B_p[i][j] - E_p[i][j] * B_t[i][j];
			double E_cross_B_t = E_p[i][j] * B_r[i][j] - E_r[i][j] * B_p[i][j];
			double E_cross_B_p = E_r[i][j] * B_t[i][j] - E_t[i][j] * B_r[i][j];
			
			//--- Time derivatives ---
			B_r_dT[i][j] = - curl_E_r[i][j];
			B_t_dT[i][j] = - curl_E_t[i][j];
			B_p_dT[i][j] = - curl_E_p[i][j];
			
			E_r_dT[i][j] = curl_B_r[i][j] + alpha[i][j] * B_r[i][j] + beta[i][j] * E_cross_B_r;
			E_t_dT[i][j] = curl_B_t[i][j] + alpha[i][j] * B_t[i][j] + beta[i][j] * E_cross_B_t;
			E_p_dT[i][j] = curl_B_p[i][j] + alpha[i][j] * B_p[i][j] + beta[i][j] * E_cross_B_p;
			
		}
	}
}




void calculate_radial_derivatives_of_VSH_coeffs(){
	
	for( int L=0; L<=ell_max; L++ ){
		
		// Right derivatives at inner boundary.
		E_r_L_dr [0][L] = ( E_r_L[1][L] - E_r_L[0][L] ) / delta_r;
		E_1_L_dr [0][L] = ( E_1_L[1][L] - E_1_L[0][L] ) / delta_r;
		E_2_L_dr [0][L] = ( E_2_L[1][L] - E_2_L[0][L] ) / delta_r;
		B_r_L_dr [0][L] = ( B_r_L[1][L] - B_r_L[0][L] ) / delta_r;
		B_2_L_dr [0][L] = ( B_2_L[1][L] - B_2_L[0][L] ) / delta_r;
		
		B_r_L_dr2[0][L] = ( B_r_L[2][L] - 2.0*B_r_L[1][L] + B_r_L[0][L] ) / ( delta_r*delta_r );
		
		// Symmetric derivatives for all intermediate points.
		for( int i=1; i<n_points_r-1; i++ ){
			E_r_L_dr [i][L] = ( E_r_L[i+1][L] - E_r_L[i-1][L] ) / ( 2.0*delta_r );
			E_1_L_dr [i][L] = ( E_1_L[i+1][L] - E_1_L[i-1][L] ) / ( 2.0*delta_r );
			E_2_L_dr [i][L] = ( E_2_L[i+1][L] - E_2_L[i-1][L] ) / ( 2.0*delta_r );
			B_r_L_dr [i][L] = ( B_r_L[i+1][L] - B_r_L[i-1][L] ) / ( 2.0*delta_r );
			B_2_L_dr [i][L] = ( B_2_L[i+1][L] - B_2_L[i-1][L] ) / ( 2.0*delta_r );
			
			B_r_L_dr2[i][L] = ( B_r_L[i+1][L] - 2.0*B_r_L[i][L] + B_r_L[i-1][L] ) / ( delta_r*delta_r );
		}
		
		// Left derivatives at outer boundary.
		E_r_L_dr .back()[L] = ( E_r_L.back()[L] - E_r_L[n_points_r-2][L] ) / delta_r;
		E_1_L_dr .back()[L] = ( E_1_L.back()[L] - E_1_L[n_points_r-2][L] ) / delta_r;
		E_2_L_dr .back()[L] = ( E_2_L.back()[L] - E_2_L[n_points_r-2][L] ) / delta_r;
		B_r_L_dr .back()[L] = ( B_r_L.back()[L] - B_r_L[n_points_r-2][L] ) / delta_r;
		B_2_L_dr .back()[L] = ( B_2_L.back()[L] - B_2_L[n_points_r-2][L] ) / delta_r;
		
		B_r_L_dr2.back()[L] = ( B_r_L.back()[L]  - 2.0*B_r_L[n_points_r-2][L] + B_r_L[n_points_r-3][L] ) / ( delta_r*delta_r );
		
	}

}




void calculate_B_1_L_and_B_1_L_dr_initial(){
	// This version of the function is for T=0, at which time we DON'T want to update B^{(1),L} because it has been read-in from a header file.
	
	for( int L=1; L<=ell_max; L++ ){
		for( int i=0; i<n_points_r; i++ ){
			B_1_L_dr[i][L] = ( r[i] * B_r_L_dr2[i][L] + 3.0 * B_r_L[i][L] ) / ( L*(L+1.0) );
		}
	}
	
}
void calculate_B_1_L_and_B_1_L_dr(){
	// This version of the function is for T>0, at which time we want to update B^{(1),L}.
	
	for( int L=1; L<=ell_max; L++ ){
		for( int i=0; i<n_points_r; i++ ){
			B_1_L[i][L] = ( r[i] * B_r_L_dr[i][L] + 2.0 * B_r_L[i][L] ) / ( L*(L+1.0) );
			B_1_L_dr[i][L] = ( r[i] * B_r_L_dr2[i][L] + 3.0 * B_r_L[i][L] ) / ( L*(L+1.0) );
		}
	}
	
}
	
	


void calculate_force_free_conditions(){
	// For output to CSV file only.
	for( int i=0; i<n_points_r; i++ ){
		for( int j=0; j<n_points_t; j++ ){
			E_dot_B      [i][j] = E_r[i][j]*B_r[i][j] + E_t[i][j]*B_t[i][j] + E_p[i][j]*B_p[i][j];
			Bsq_minus_Esq[i][j] = B_r[i][j]*B_r[i][j] + B_t[i][j]*B_t[i][j] + B_p[i][j]*B_p[i][j] - E_r[i][j]*E_r[i][j] - E_t[i][j]*E_t[i][j] - E_p[i][j]*E_p[i][j];
		}
	}
}




void integrate_vectors_wrt_time(){
	// For now, just use the Euler method. Update to Adams-Bashforth later.
	for( int i=0; i<n_points_r; i++ ){
		for( int j=0; j<n_points_t; j++ ){
			B_r[i][j] += B_r_dT[i][j] * delta_T;
			B_t[i][j] += B_t_dT[i][j] * delta_T;
			B_p[i][j] += B_p_dT[i][j] * delta_T;
			E_r[i][j] += E_r_dT[i][j] * delta_T;
			E_t[i][j] += E_t_dT[i][j] * delta_T;
			E_p[i][j] += E_p_dT[i][j] * delta_T;
		}
	}
}




void apply_inner_boundary_conditions(){
	// Petri 2012, \S3.3.
	for( int j=0; j<n_points_t; j++ ){
		B_r[0][j] = 2.0 * ct[j];
		E_t[0][j] = - Omega * 2.0 * st[j] * ct[j];
		E_p[0][j] = 0;
	}
}
void apply_outer_boundary_conditions(){
	// Petri 2012, \S3.3.
	for( int j=0; j<n_points_t; j++ ){
		B_t.back()[j] = 0.5 * ( B_t.back()[j] - E_p.back()[j] );
		B_p.back()[j] = 0.5 * ( E_t.back()[j] + B_p.back()[j] );
		E_t.back()[j] =   B_p.back()[j];
		E_p.back()[j] = - B_t.back()[j];
	}
}




void calculate_stdev(){
	
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
			for( int ell=0; ell<=ell_max; ell++ ){
				
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
			
			sqdev_B_r[i][j] = pow( B_r[i][j] - B_r_series[i][j], 2 );
			sqdev_B_t[i][j] = pow( B_t[i][j] - B_t_series[i][j], 2 );
			sqdev_B_p[i][j] = pow( B_p[i][j] - B_p_series[i][j], 2 );
			sqdev_E_r[i][j] = pow( E_r[i][j] - E_r_series[i][j], 2 );
			sqdev_E_t[i][j] = pow( E_t[i][j] - E_t_series[i][j], 2 );
			sqdev_E_p[i][j] = pow( E_p[i][j] - E_p_series[i][j], 2 );
			
			stdev_B_r += sqdev_B_r[i][j];
			stdev_B_t += sqdev_B_t[i][j];
			stdev_B_p += sqdev_B_p[i][j];
			stdev_E_r += sqdev_E_r[i][j];
			stdev_E_t += sqdev_E_t[i][j];
			stdev_E_p += sqdev_E_p[i][j];
			
			if( sqdev_B_r[i][j] > max_sqdev_B_r ){ max_sqdev_B_r = sqdev_B_r[i][j]; }
			if( sqdev_B_t[i][j] > max_sqdev_B_t ){ max_sqdev_B_t = sqdev_B_t[i][j]; }
			if( sqdev_B_p[i][j] > max_sqdev_B_p ){ max_sqdev_B_p = sqdev_B_p[i][j]; }
			if( sqdev_E_r[i][j] > max_sqdev_E_r ){ max_sqdev_E_r = sqdev_E_r[i][j]; }
			if( sqdev_E_t[i][j] > max_sqdev_E_t ){ max_sqdev_E_t = sqdev_E_t[i][j]; }
			if( sqdev_E_p[i][j] > max_sqdev_E_p ){ max_sqdev_E_p = sqdev_E_p[i][j]; }
			
		}
	}
	
	stdev_B_r = sqrt( stdev_B_r / ( (double) n_points_r * n_points_t ) );
	stdev_B_t = sqrt( stdev_B_t / ( (double) n_points_r * n_points_t ) );
	stdev_B_p = sqrt( stdev_B_p / ( (double) n_points_r * n_points_t ) );
	stdev_E_r = sqrt( stdev_E_r / ( (double) n_points_r * n_points_t ) );
	stdev_E_t = sqrt( stdev_E_t / ( (double) n_points_r * n_points_t ) );
	stdev_E_p = sqrt( stdev_E_p / ( (double) n_points_r * n_points_t ) );
	
}




void calculate_magnetic_energy(){
	
	for( int i=0; i<n_points_r; i++ ){
		for( int j=0; j<n_points_t; j++ ){
			magnetic_energy_density[i][j] = B_r[i][j]*B_r[i][j] + B_t[i][j]*B_t[i][j] + B_p[i][j]*B_p[i][j];
		}
	}
	
	total_magnetic_energy = integral_2d( magnetic_energy_density, integral_r_t_dr, integral_r_t_dt ) * 2.0 * pi;
	
}




double estimate_csv_profiles_filesize(){
	
	int n_columns = 50;		// You are responsible for keeping this value updated.
	int n_rows = 0;
	double fudge_factor = 10.0;
	
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
	output_file_profiles << ",sqdev_B_r,sqdev_B_t,sqdev_B_p,sqdev_E_r,sqdev_E_t,sqdev_E_p";
	output_file_profiles << ",magnetic_energy_density";
	output_file_profiles << ",curl_B_r,curl_B_t,curl_B_p";
	output_file_profiles << "\n";
}


int sgn( double x ){
	//signum function with sgn(x)=0 for x=0.
	return ( x > 0 ) - ( x < 0 );
}

void output_to_csv_profiles( int T_index ){
	
	bool condition_T_1 = ( T_index >= csv_profiles_write_T_min ) and ( T_index <= csv_profiles_write_T_max );
	bool condition_T_2 = ( T_index % csv_profiles_write_freq_T == 0 ) or ( T_index == n_timesteps );
	
	if( condition_T_1 and condition_T_2 ){
	
		for( int i=csv_profiles_write_i_min; i<n_points_r; i++ ){
			for( int j=0; j<n_points_t; j++ ){
				
				bool condition_r = ( i <= csv_profiles_write_i_max ) and ( ( i % csv_profiles_write_freq_r == 0 ) or ( i == n_points_r - 1 ) );
				bool condition_t = ( j % csv_profiles_write_freq_t == 0 ) or ( j == n_points_t - 1 );
		
				if( condition_r and condition_t ){
					
					double B = sqrt( pow( B_r[i][j], 2 ) + pow( B_t[i][j], 2 ) + pow( B_p[i][j], 2 ) );
					double E = sqrt( pow( E_r[i][j], 2 ) + pow( E_t[i][j], 2 ) + pow( E_p[i][j], 2 ) );
					
					output_file_profiles << T_index <<","<< T <<","<< T*T_factor <<","<< i <<","<< j <<","<< r[i] <<","<< t[j] <<","<< x[i][j] <<","<< z[i][j]
										 <<","<< B_x[i][j] <<","<< B_y[i][j] <<","<< B_z[i][j] <<","<< B
										 <<","<< E_x[i][j] <<","<< E_y[i][j] <<","<< E_z[i][j] <<","<< E
										 <<","<< alpha[i][j] <<","<< beta[i][j]
										 <<","<< E_dot_B[i][j] <<","<< Bsq_minus_Esq[i][j] <<","<< sgn( Bsq_minus_Esq[i][j] )
										 <<","<< B_r_dr[i][j] <<","<< B_t_dr[i][j] <<","<<  B_p_dr[i][j] <<","<< E_r_dr[i][j] <<","<< E_t_dr[i][j] <<","<<  E_p_dr[i][j]
										 <<","<< B_r       [i][j] <<","<< B_t       [i][j] <<","<< B_p       [i][j] <<","<< E_r       [i][j] <<","<< E_t       [i][j] <<","<< E_p       [i][j]
										 <<","<< B_r_series[i][j] <<","<< B_t_series[i][j] <<","<< B_p_series[i][j] <<","<< E_r_series[i][j] <<","<< E_t_series[i][j] <<","<< E_p_series[i][j]
										 <<","<< sqdev_B_r [i][j] <<","<< sqdev_B_t [i][j] <<","<< sqdev_B_p [i][j] <<","<< sqdev_E_r [i][j] <<","<< sqdev_E_t [i][j] <<","<< sqdev_E_p [i][j]
										 <<","<< magnetic_energy_density[i][j]
										 <<","<< curl_B_r[i][j] <<","<< curl_B_t[i][j] <<","<< curl_B_p[i][j]
										 << "\n";
					
				}
			}
		}
	}
}




void output_headers_to_csv_history(){
	// Split over separate lines for readability and ease of checking consistency with output_to_csv_history().
	
	output_file_history << "T_index,T,T_sec";
	output_file_history << ",pct_E_pts_chngd";
	output_file_history << ",stdev_B_r,stdev_B_t,stdev_B_p,stdev_E_r,stdev_E_t,stdev_E_p";
	output_file_history << ",max_sqdev_B_r,max_sqdev_B_t,max_sqdev_B_p,max_sqdev_E_r,max_sqdev_E_t,max_sqdev_E_p";
	output_file_history << ",total_magnetic_energy";
	output_file_history << "\n";
}




void output_to_csv_history( int T_index ){
	
	bool condition_T = ( T_index % csv_history_write_freq_T == 0 ) or ( T_index == n_timesteps );
	
	if( condition_T ){
		
		output_file_history << T_index <<","<< T <<","<< T*T_factor
							<<","<< percentage_E_points_points_changed_for_second_force_free_condition
							<<","<< stdev_B_r     <<","<< stdev_B_t     <<","<< stdev_B_p     <<","<< stdev_E_r     <<","<< stdev_E_t     <<","<< stdev_E_p
							<<","<< max_sqdev_B_r <<","<< max_sqdev_B_t <<","<< max_sqdev_B_p <<","<< max_sqdev_E_r <<","<< max_sqdev_E_t <<","<< max_sqdev_E_p
							<<","<< total_magnetic_energy
							<<"\n";
							
	}
}




void output_headers_to_screen(){
	std::cout <<std::left<<std::setw(w)<< "T_index" <<std::left<<std::setw(w)<< "T" <<"|\t"
	          <<std::left<<std::setw(w)<< "B_r" <<std::left<<std::setw(w)<< "B_t" <<std::left<<std::setw(w)<< "B_p" <<"|\t"
			  <<std::left<<std::setw(w)<< "E_r" <<std::left<<std::setw(w)<< "E_t" <<std::left<<std::setw(w)<< "E_p" <<"|\t"
			  <<std::left<<std::setw(w)<< "E dot B" <<std::left<<std::setw(w)<< "B^2 - E^2"
			  <<std::left<<std::setw(w)<< "% E pts chngd" <<"|\t"
			  <<std::left<<std::setw(w)<< "% done" <<std::left<<std::setw(w)<< "Est. time left"
			  << std::endl;
}




void output_to_screen( int T_index ){
	
	if( ( T_index % cout_freq_T == 0 ) or ( T_index == n_timesteps ) ){
	
		double percent_done = round( (double)1000*T_index/n_timesteps )/10; // This format makes it exactly one decimal place.
		
		std::cout <<std::left<<std::setw(w)<< T_index <<std::left<<std::setw(w)<< T <<"|\t"
				  <<std::left<<std::setw(w)<< B_r[cout_i][cout_j] <<std::left<<std::setw(w)<< B_t[cout_i][cout_j] <<std::left<<std::setw(w)<< B_p[cout_i][cout_j] <<"|\t"
				  <<std::left<<std::setw(w)<< E_r[cout_i][cout_j] <<std::left<<std::setw(w)<< E_t[cout_i][cout_j] <<std::left<<std::setw(w)<< E_p[cout_i][cout_j] <<"|\t"
				  <<std::left<<std::setw(w)<< E_dot_B[cout_i][cout_j] <<std::left<<std::setw(w)<< Bsq_minus_Esq[cout_i][cout_j]
				  <<std::left<<std::setw(w)<< percentage_E_points_points_changed_for_second_force_free_condition <<"|\t"
				  <<std::left<<std::setw(w)<< percent_done;
		
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




void output_execution_time(){
	
	double   time_stop_seconds   = std::chrono::high_resolution_clock::now().time_since_epoch().count() * 1e-9;
	double exec_time   = time_stop_seconds - time_start_seconds;
	int    exec_time_h = exec_time / 3600;
	int    exec_time_m = exec_time / 60 - exec_time_h * 60;
	int    exec_time_s = exec_time - exec_time_h * 3600 - exec_time_m * 60;
	std::cout << "\nExecution time (h:mm:ss):\t" << exec_time_h;
	if( exec_time_m < 10 ){
		std::cout << ":0" << exec_time_m;
	} else {
		std::cout << ":"  << exec_time_m;
	}
	if( exec_time_s < 10 ){
		std::cout << ":0" << exec_time_s << std::endl;
	} else {
		std::cout << ":"  << exec_time_s << std::endl;
	}
	
}