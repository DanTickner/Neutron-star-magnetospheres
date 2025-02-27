/*
Header_Time_Evolution_03.h

Code snippets for the Time_Evolution_xx.cpp codes. Keeps them shorter and allows for two separate for-loops with Euler and Adams-Bashforth integration
without copying code.

V03: Apply second force-free condition before VSH decomposition. Then, we don't need an extra step in that algorithm to redo it.


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

//std::vector< std::vector< std::vector<double> > > Y_L1_cross_Psi_L2_coeff_2_L   ( ell_max+1, std::vector< std::vector<double> > ( ell_max+1, std::vector<double> ( ell_max+1 ) ) );
//std::vector< std::vector< std::vector<double> > > Y_L1_cross_Phi_L2_coeff_1_L   ( ell_max+1, std::vector< std::vector<double> > ( ell_max+1, std::vector<double> ( ell_max+1 ) ) );
//std::vector< std::vector< std::vector<double> > > Psi_L1_cross_Phi_L2_coeff_r_L ( ell_max+1, std::vector< std::vector<double> > ( ell_max+1, std::vector<double> ( ell_max+1 ) ) );

std::vector< std::vector<double> > B_r_L   ( n_points_r, std::vector<double> ( ell_max+1 ) );
std::vector< std::vector<double> > B_1_L   ( n_points_r, std::vector<double> ( ell_max+1 ) );
std::vector< std::vector<double> > B_2_L   ( n_points_r, std::vector<double> ( ell_max+1 ) );
std::vector< std::vector<double> > E_r_L   ( n_points_r, std::vector<double> ( ell_max+1 ) );
std::vector< std::vector<double> > E_1_L   ( n_points_r, std::vector<double> ( ell_max+1 ) );
std::vector< std::vector<double> > E_2_L   ( n_points_r, std::vector<double> ( ell_max+1 ) );
//std::vector< std::vector<double> > ExB_r_L ( n_points_r, std::vector<double> ( ell_max+1 ) );
//std::vector< std::vector<double> > ExB_1_L ( n_points_r, std::vector<double> ( ell_max+1 ) );
//std::vector< std::vector<double> > ExB_2_L ( n_points_r, std::vector<double> ( ell_max+1 ) );

std::vector< std::vector<double> > B_r_L_dr  ( n_points_r, std::vector<double> ( ell_max+1 ) );
std::vector< std::vector<double> > B_1_L_dr  ( n_points_r, std::vector<double> ( ell_max+1 ) );
std::vector< std::vector<double> > B_2_L_dr  ( n_points_r, std::vector<double> ( ell_max+1 ) );
std::vector< std::vector<double> > E_r_L_dr  ( n_points_r, std::vector<double> ( ell_max+1 ) );
std::vector< std::vector<double> > E_1_L_dr  ( n_points_r, std::vector<double> ( ell_max+1 ) );
std::vector< std::vector<double> > E_2_L_dr  ( n_points_r, std::vector<double> ( ell_max+1 ) );
std::vector< std::vector<double> > B_r_L_dr2 ( n_points_r, std::vector<double> ( ell_max+1 ) );

std::vector< std::vector<double> > B_r    ( n_points_r, std::vector<double> ( n_points_t ) );
std::vector< std::vector<double> > B_t    ( n_points_r, std::vector<double> ( n_points_t ) );
std::vector< std::vector<double> > B_p    ( n_points_r, std::vector<double> ( n_points_t ) );
std::vector< std::vector<double> > E_r    ( n_points_r, std::vector<double> ( n_points_t ) );
std::vector< std::vector<double> > E_t    ( n_points_r, std::vector<double> ( n_points_t ) );
std::vector< std::vector<double> > E_p    ( n_points_r, std::vector<double> ( n_points_t ) );
std::vector< std::vector<double> > B_r_dT ( n_points_r, std::vector<double> ( n_points_t ) );
std::vector< std::vector<double> > B_t_dT ( n_points_r, std::vector<double> ( n_points_t ) );
std::vector< std::vector<double> > B_p_dT ( n_points_r, std::vector<double> ( n_points_t ) );
std::vector< std::vector<double> > E_r_dT ( n_points_r, std::vector<double> ( n_points_t ) );
std::vector< std::vector<double> > E_t_dT ( n_points_r, std::vector<double> ( n_points_t ) );
std::vector< std::vector<double> > E_p_dT ( n_points_r, std::vector<double> ( n_points_t ) );
std::vector< std::vector<double> > B_x    ( n_points_r, std::vector<double> ( n_points_t ) );
std::vector< std::vector<double> > B_y    ( n_points_r, std::vector<double> ( n_points_t ) );
std::vector< std::vector<double> > B_z    ( n_points_r, std::vector<double> ( n_points_t ) );
std::vector< std::vector<double> > E_x    ( n_points_r, std::vector<double> ( n_points_t ) );
std::vector< std::vector<double> > E_y    ( n_points_r, std::vector<double> ( n_points_t ) );
std::vector< std::vector<double> > E_z    ( n_points_r, std::vector<double> ( n_points_t ) );
std::vector< std::vector<double> > B_r_dr ( n_points_r, std::vector<double> ( n_points_t ) );	// Not used by the code, but useful within evaluate_radial_derivatives() for CSV output to check for strange behaviour in dB/dr and dE/dr.
std::vector< std::vector<double> > B_t_dr ( n_points_r, std::vector<double> ( n_points_t ) );
std::vector< std::vector<double> > B_p_dr ( n_points_r, std::vector<double> ( n_points_t ) );
std::vector< std::vector<double> > E_r_dr ( n_points_r, std::vector<double> ( n_points_t ) );
std::vector< std::vector<double> > E_t_dr ( n_points_r, std::vector<double> ( n_points_t ) );
std::vector< std::vector<double> > E_p_dr ( n_points_r, std::vector<double> ( n_points_t ) );

std::vector< std::vector<double> > E_dot_B      ( n_points_r, std::vector<double> ( n_points_t ) );		// To check that force-free conditions are fulfilled.
std::vector< std::vector<double> > Bsq_minus_Esq( n_points_r, std::vector<double> ( n_points_t ) );		// To check that force-free conditions are fulfilled.

std::vector< std::vector<double> > alpha             ( n_points_r, std::vector<double> ( n_points_t ) );	// Don't have to be 2D array but might want to track in CSV file.
std::vector< std::vector<double> > beta              ( n_points_r, std::vector<double> ( n_points_t ) );	// Don't have to be 2D array but might want to track in CSV file.
//std::vector< std::vector<double> > E_dot_B_over_Bsq  ( n_points_r, std::vector<double> ( n_points_t ) );	// Don't have to be 2D array but might want to track in CSV file.

//std::vector<double> alpha_integrand_0 ( n_points_t );
//std::vector<double> alpha_integrand_1 ( n_points_t );
//std::vector<double> beta_integrand_0  ( n_points_t );
//std::vector<double> beta_integrand_1  ( n_points_t );

std::vector<double> VSH_integrand_r ( n_points_t );
std::vector<double> VSH_integrand_1 ( n_points_t );
std::vector<double> VSH_integrand_2 ( n_points_t );

//std::vector<double> E_VSH_integrand_r ( n_points_t );
//std::vector<double> E_VSH_integrand_1 ( n_points_t );
//std::vector<double> E_VSH_integrand_2 ( n_points_t );
//std::vector< std::vector<double> > alpha_integrals_0( ell_max+1, std::vector<double> ( ell_max+1 ) );
//std::vector< std::vector<double> > alpha_integrals_1( ell_max+1, std::vector<double> ( ell_max+1 ) );
//std::vector< std::vector<double> > beta_integrals_0 ( ell_max+1, std::vector<double> ( ell_max+1 ) );
//std::vector< std::vector<double> > beta_integrals_1 ( ell_max+1, std::vector<double> ( ell_max+1 ) );

//std::vector<double> B_r_L_dT_i ( ell_max+1 );
//std::vector<double> B_2_L_dT_i ( ell_max+1 );
//std::vector<double> E_r_L_dT_i ( ell_max+1 );
//std::vector<double> E_1_L_dT_i ( ell_max+1 );
//std::vector<double> E_2_L_dT_i ( ell_max+1 );

//std::vector<double> B_tilde_integrand_0     ( n_points_t );
//std::vector<double> B_tilde_integrand_1     ( n_points_t );
//std::vector<double> E_corrected_integrand_0 ( n_points_t );
//std::vector<double> E_corrected_integrand_1 ( n_points_t );
//std::vector< std::vector<double> > B_tilde_integrals_0    ( ell_max+1, std::vector<double> ( ell_max+1 ) );
//std::vector< std::vector<double> > B_tilde_integrals_1    ( ell_max+1, std::vector<double> ( ell_max+1 ) );
//std::vector< std::vector<double> > E_corrected_integrals_0( ell_max+1, std::vector<double> ( ell_max+1 ) );
//std::vector< std::vector<double> > E_corrected_integrals_1( ell_max+1, std::vector<double> ( ell_max+1 ) );

//std::vector< std::vector<double> > E_tilde_r_L     ( n_points_r, std::vector<double> ( ell_max+1 ) );
//std::vector< std::vector<double> > E_tilde_1_L     ( n_points_r, std::vector<double> ( ell_max+1 ) );
//std::vector< std::vector<double> > E_tilde_2_L     ( n_points_r, std::vector<double> ( ell_max+1 ) );
//std::vector< std::vector<double> > E_corrected_r_L ( n_points_r, std::vector<double> ( ell_max+1 ) );
//std::vector< std::vector<double> > E_corrected_1_L ( n_points_r, std::vector<double> ( ell_max+1 ) );
//std::vector< std::vector<double> > E_corrected_2_L ( n_points_r, std::vector<double> ( ell_max+1 ) );

std::ofstream output_file;

double time_start_seconds = 0;

double T = 0;

int n_gridpoints_inside_light_cylinder = 0;
double percentage_E_points_points_changed_for_second_force_free_condition = 0;	// Output to CSV to keep track of how active the FF condition application is.




//----- Functions -----


void calculate_gridpoints(){
	
	for( int i=0; i<n_points_r; i++ ){
		r[i] = r_min + i * delta_r;
	}
	for( int j=0; j<n_points_t; j++ ){
		t[j] = j * delta_t;
		st[j] = sin( t[j] );
		ct[j] = cos( t[j] );
	}
	
	for( int i=0; i<n_points_r; i++ ){
		for( int j=0; j<n_points_t; j++ ){
			x[i][j] = r[i] * st[j];
			z[i][j] = r[i] * ct[j];
		}
	}
	
}




void calculate_associated_legendre_functions(){
	
	for( int j=0; j<n_points_t; j++ ){
		
		double x = cos( t[j] );
		
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
		
				double sqrt_factor = sqrt( (2.0*ell+1.0) / ( 4.0 * pi ) );
				
				B_r[i][j] += sqrt_factor * P0[j][ell] * B_r_L[i][ell];
				B_t[i][j] += sqrt_factor * P1[j][ell] * B_1_L[i][ell];
				B_t[i][j] += sqrt_factor * P1[j][ell] * B_2_L[i][ell];
				
				E_r[i][j] += sqrt_factor * P0[j][ell] * E_r_L[i][ell];
				E_t[i][j] += sqrt_factor * P1[j][ell] * E_1_L[i][ell];
				E_t[i][j] += sqrt_factor * P1[j][ell] * E_2_L[i][ell];
				
			}
		}
	}
}

void evaluate_radial_derivatives(){
	// Not used by the code, but useful for CSV output to check for strange behaviour in dB/dr and dE/dr.
	// Do second derivatives too.
	
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




void apply_initial_field_values_and_VSH_decompositions(){
	
	for( int i=0; i<n_points_r; i++ ){
		
		for( int j=0; j<n_points_t; j++ ){
			B_r[i][j] = B_r_function( r[i], t[j] );
			B_t[i][j] = B_t_function( r[i], t[j] );
			B_p[i][j] = B_p_function( r[i], t[j] );
			E_r[i][j] = E_r_function( r[i], t[j] );
			E_t[i][j] = E_t_function( r[i], t[j] );
			E_p[i][j] = E_p_function( r[i], t[j] );
		}
		
		for( int L=0; L<=ell_max; L++ ){
			B_r_L[i][L] = B_r_L_function( r[i], L );
			B_1_L[i][L] = B_1_L_function( r[i], L );
			B_2_L[i][L] = B_2_L_function( r[i], L );
			E_r_L[i][L] = E_r_L_function( r[i], L );
			E_1_L[i][L] = E_1_L_function( r[i], L );
			E_2_L[i][L] = E_2_L_function( r[i], L );
		}
		
	}
}




void calculate_n_gridpoints_inside_light_cylinder(){
	for( int i=0; i<n_points_r; i++ ){
		for( int j=0; j<n_points_t; j++ ){
			if( r[i]*st[j] < R_LC ){
				n_gridpoints_inside_light_cylinder++;
			}
		}
	}
}



void apply_second_force_free_condition_within_light_cylinder(){
	// Petri 2012, Section 3.2.
	
	//--- Update the electric field at all points within the light cylinder if the second FF condition B^2-E^2<0 is violated ---
	int num_E_points_changed_for_second_force_free_condition = 0;
	
	for( int i=0; i<n_points_r; i++ ){
		for( int j=0; j<n_points_t; j++ ){
			if( r[i]*st[j] < R_LC ){
				
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
	percentage_E_points_points_changed_for_second_force_free_condition = round( (double)1000 * num_E_points_changed_for_second_force_free_condition / n_gridpoints_inside_light_cylinder ) / 10;
}




void calculate_time_derivatives(){
	for( int i=0; i<n_points_r; i++ ){
		for( int j=0; j<n_points_t; j++ ){
			
			//--- Spatial derivatives ---
			double div_E    = 0;
			double curl_E_r = 0;
			double curl_E_t = 0;
			double curl_E_p = 0;
			double curl_B_r = 0;
			double curl_B_t = 0;
			double curl_B_p = 0;
			
			//- L=0 -
			div_E += sqrt_2Lplus1_over_4pi[0] * ( E_r_L_dr[i][0] + 2.0*pow(r[i],-1) * E_r_L[i][0] );
			
			//- L>0 -
			for( int L=1; L<=ell_max; L++ ){
				
				div_E += sqrt_2Lplus1_over_4pi[L] * ( E_r_L_dr[i][L] + 2.0*pow(r[i],-1) * E_r_L[i][L] - L*(L+1.0)*pow(r[i],-1) * E_1_L[i][L] );
				
				curl_E_r += sqrt_2Lplus1_over_4pi[L] * P0[j][L] * -L*(L+1) * E_2_L[i][L] / r[i];
				curl_E_t += sqrt_2Lplus1_over_4pi[L] * P1[j][L] * ( -E_2_L_dr[i][L] - E_2_L[i][L] / r[i] );
				curl_E_p += sqrt_2Lplus1_over_4pi[L] * P1[j][L] * ( E_1_L_dr[i][L] + ( E_1_L[i][L] - E_r_L[i][L] ) / r[i] );
				
				curl_B_r += sqrt_2Lplus1_over_4pi[L] * P0[j][L] * -L*(L+1) * B_2_L[i][L] / r[i];
				curl_B_t += sqrt_2Lplus1_over_4pi[L] * P1[j][L] * ( -B_2_L_dr[i][L] - B_2_L[i][L] / r[i] );
				curl_B_p += sqrt_2Lplus1_over_4pi[L] * P1[j][L] * ( r[i] * B_r_L_dr2[i][L] + 4.0 * B_r_L_dr[i][L] - (L+2.0)*(L-1.0) * B_r_L[i][L] / r[i] ) / ( L*(L+1.0) );
				
			}
			
			//--- alpha and beta ---
			double B_dot_curl_B = B_r[i][j]*curl_B_r  + B_t[i][j]*curl_B_t  + B_p[i][j]*curl_B_p;
			double E_dot_curl_E = E_r[i][j]*curl_E_r  + E_t[i][j]*curl_E_t  + E_p[i][j]*curl_E_p;
			double B_squared    = B_r[i][j]*B_r[i][j] + B_t[i][j]*B_t[i][j] + B_p[i][j]*B_p[i][j];
			double E_squared    = E_r[i][j]*E_r[i][j] + E_t[i][j]*E_t[i][j] + E_p[i][j]*E_p[i][j];
			
			alpha[i][j] = ( E_dot_curl_E - B_dot_curl_B ) / B_squared;
			beta [i][j]  = - div_E / B_squared;
			
			//--- E cross B ---
			double E_cross_B_r = E_t[i][j] * B_p[i][j] - E_p[i][j] * E_t[i][j];
			double E_cross_B_t = E_p[i][j] * B_r[i][j] - E_r[i][j] * E_p[i][j];
			double E_cross_B_p = E_r[i][j] * B_t[i][j] - E_t[i][j] * E_r[i][j];
			
			//--- Time derivatives ---
			B_r_dT[i][j] = - curl_E_r;
			B_t_dT[i][j] = - curl_E_t;
			B_p_dT[i][j] = - curl_E_p;
			
			E_r_dT[i][j] = curl_B_r + alpha[i][j] * B_r[i][j] + beta[i][j] * E_cross_B_r;
			E_r_dT[i][j] = curl_B_t + alpha[i][j] * B_t[i][j] + beta[i][j] * E_cross_B_t;
			E_r_dT[i][j] = curl_B_p + alpha[i][j] * B_p[i][j] + beta[i][j] * E_cross_B_p;
			
		}
	}
}
			








void calculate_radial_derivatives_of_VSH_coeffs_inner_boundary(){
	// Use right derivatives at inner boundary.
	// Inner BCs mean fields don't change there, so this only needs to be called once. Hence, move into separate function outside loop over time.
	
	for( int L=0; L<=ell_max; L++ ){		
			E_r_L_dr [0][L] = ( E_r_L[1][L] - E_r_L[0][L] ) / delta_r;
			E_1_L_dr [0][L] = ( E_1_L[1][L] - E_1_L[0][L] ) / delta_r;
			E_2_L_dr [0][L] = ( E_2_L[1][L] - E_2_L[0][L] ) / delta_r;
			B_r_L_dr [0][L] = ( B_r_L[1][L] - B_r_L[0][L] ) / delta_r;
			B_2_L_dr [0][L] = ( B_2_L[1][L] - B_2_L[0][L] ) / delta_r;
			
			B_r_L_dr2[0][L] = ( B_r_L[2][L] - 2.0*B_r_L[1][L] + B_r_L[0][L] ) / ( delta_r*delta_r );
	}
}

void calculate_radial_derivatives_of_VSH_coeffs(){
	// Would use right derivatives at inner boundary, but the inner BCs mean fields don't change there, so no need to recalculate at i=0.
	
	for( int L=0; L<=ell_max; L++ ){
		
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
		E_r_L_dr .back()[L] = ( E_r_L.back()[L] - E_r_L[n_points_r-1][L] ) / delta_r;
		E_1_L_dr .back()[L] = ( E_1_L.back()[L] - E_1_L[n_points_r-1][L] ) / delta_r;
		E_2_L_dr .back()[L] = ( E_2_L.back()[L] - E_2_L[n_points_r-1][L] ) / delta_r;
		B_r_L_dr .back()[L] = ( B_r_L.back()[L] - B_r_L[n_points_r-1][L] ) / delta_r;
		B_2_L_dr .back()[L] = ( B_2_L.back()[L] - B_2_L[n_points_r-1][L] ) / delta_r;
		
		B_r_L_dr2.back()[L] = ( B_r_L.back()[L]  - 2.0*B_r_L[n_points_r-1][L] + B_r_L[n_points_r-2][L] ) / ( delta_r*delta_r );
		
	}

}




void calculate_B_1_L_and_B_1_L_dr_inner_boundary(){
	// Inner BCs mean fields don't change there, so this only needs to be called once (at T=0). Hence, move into separate function outside loop over time.
	for( int L=1; L<=ell_max; L++ ){
		B_1_L_dr[0][L] = ( r[0] * B_r_L_dr2[0][L] + 3.0 * B_r_L[0][L] ) / ( L*(L+1.0) );
	}
	
}


void calculate_B_1_L_and_B_1_L_dr_initial(){
	// This version of the function is for T=0, at which time we DON'T want to update B^{(1),L} because it has been read-in from a header file.
	
	for( int L=1; L<=ell_max; L++ ){
		for( int i=1; i<n_points_r; i++ ){
			B_1_L_dr[i][L] = ( r[i] * B_r_L_dr2[i][L] + 3.0 * B_r_L[i][L] ) / ( L*(L+1.0) );
		}
	}
	
}
void calculate_B_1_L_and_B_1_L_dr(){
	// This version of the function is for T>0, at which time we want to update B^{(1),L}.
	
	for( int L=1; L<=ell_max; L++ ){
		for( int i=1; i<n_points_r; i++ ){
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




void apply_outer_boundary_conditions(){
	// Petri 2012, \S3.3.
	// Verified that both versions give the same output 20230818. Only use the simpler V2 going forward.
	for( int j=0; j<n_points_t; j++ ){
		
		
		//--- V1 ---
		double B_theta_PDE = B_t.back()[j];
		double B_phi_PDE   = B_p.back()[j];
		double E_theta_PDE = E_t.back()[j];
		double E_phi_PDE   = E_p.back()[j];
		B_t.back()[j] = 0.5 * ( B_theta_PDE - E_phi_PDE );
		E_p.back()[j] = - B_t.back()[j];
		B_p.back()[j] = 0.5 * ( E_theta_PDE + B_phi_PDE );
		E_t.back()[j] = B_p.back()[j];
		
		
		/*
		//--- V2 ---
		B_t.back()[j] = 0.5 * ( B_t.back()[j] - E_p.back()[j] );
		B_p.back()[j] = 0.5 * ( E_t.back()[j] + B_p.back()[j] );
		*/
	}
}


















void output_headers_to_csv(){
	// Split over separate lines for readability and ease of checking consistency with output_to_csv().
	
	//output_file << "T_index,T,T_sec,i,j,r,t,x,z,B_x,B_y,B_z,B,div(B),curl(B)_x,curl(B)_y,curl(B)_z,E_x,E_y,E_z,E,div(E),curl(E)_x,curl(E)_y,curl(E)_z,alpha,beta,E dot B,E minus B\n";
	
	output_file << "T_index,T,T_sec,i,j,r,t,x,z";
	output_file << ",B_x,B_y,B_z,B,E_x,E_y,E_z,E";
	output_file << ",alpha,beta";
	output_file << ",E_dot_B,Bsq_minus_Esq,sgn(Bsq_minus_Esq),pct_E_pts_chngd";
	output_file << ",B_r_dr,B_t_dr,B_p_dr,E_r_dr,E_t_dr,E_p_dr";
	output_file << "\n";
}


int sgn( double x ){
	//signum function with sgn(x)=0 for x=0.
	return ( x > 0 ) - ( x < 0 );
}

void output_to_csv( int T_index ){
	
	bool condition_T = ( T_index % csv_write_freq_T == 0 ) or ( T_index == n_timesteps );
	
	if( condition_T ){
	
		for( int i=0; i<n_points_r; i++ ){
			for( int j=0; j<n_points_t; j++ ){
				
				//bool condition_r = i > 960 and ( ( i % csv_write_freq_r == 0 ) or ( i == n_points_r - 1 ) );	// modified 20230818 to explore issues near outer boundary.
				bool condition_r = ( i % csv_write_freq_r == 0 ) or ( i == n_points_r - 1 );
				bool condition_t = ( j % csv_write_freq_t == 0 ) or ( j == n_points_t - 1 );
		
				if( condition_r and condition_t ){
					
					double B = sqrt( pow( B_r[i][j], 2 ) + pow( B_t[i][j], 2 ) + pow( B_p[i][j], 2 ) );
					double E = sqrt( pow( E_r[i][j], 2 ) + pow( E_t[i][j], 2 ) + pow( E_p[i][j], 2 ) );
					
					output_file << T_index <<","<< T <<","<< T*T_factor <<","<< i <<","<< j <<","<< r[i] <<","<< t[j] <<","<< x[i][j] <<","<< z[i][j]
								<<","<< B_x[i][j] <<","<< B_y[i][j] <<","<< B_z[i][j] <<","<< B
								<<","<< E_x[i][j] <<","<< E_y[i][j] <<","<< E_z[i][j] <<","<< E
								<<","<< alpha[i][j] <<","<< beta[i][j]
								<<","<< E_dot_B[i][j] <<","<< Bsq_minus_Esq[i][j] <<","<< sgn( Bsq_minus_Esq[i][j] ) <<","<< percentage_E_points_points_changed_for_second_force_free_condition
								<<","<< B_r_dr[i][j] <<","<< B_t_dr[i][j] <<","<<  B_p_dr[i][j] <<","<< E_r_dr[i][j] <<","<< E_t_dr[i][j] <<","<<  E_p_dr[i][j];
					
					output_file << "\n";
					
				}
				
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