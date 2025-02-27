/*

Header_Time_Evolution_01.h

Code snippets for the Time_Evolution_xx.cpp codes. Keeps them shorter and allows for two separate for-loops with Euler and Adams-Bashforth integration
without copying code.

V04: Do the full time-evolution at a single gridpoint, to remove a dimension from the lists.
This can become significant when there are very many gridpoints (e.g. thousands) because the lists would otherwise be (e.g. thousands) of times larger in size.


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

std::vector< std::vector< std::vector<double> > > Y_L1_cross_Psi_L2_coeff_2_L   ( ell_max+1, std::vector< std::vector<double> > ( ell_max+1, std::vector<double> ( ell_max+1 ) ) );
std::vector< std::vector< std::vector<double> > > Y_L1_cross_Phi_L2_coeff_1_L   ( ell_max+1, std::vector< std::vector<double> > ( ell_max+1, std::vector<double> ( ell_max+1 ) ) );
std::vector< std::vector< std::vector<double> > > Psi_L1_cross_Phi_L2_coeff_r_L ( ell_max+1, std::vector< std::vector<double> > ( ell_max+1, std::vector<double> ( ell_max+1 ) ) );

std::vector< std::vector<double> > B_r_L   ( n_points_r, std::vector<double> ( ell_max+1 ) );
std::vector< std::vector<double> > B_1_L   ( n_points_r, std::vector<double> ( ell_max+1 ) );
std::vector< std::vector<double> > B_2_L   ( n_points_r, std::vector<double> ( ell_max+1 ) );
std::vector< std::vector<double> > E_r_L   ( n_points_r, std::vector<double> ( ell_max+1 ) );
std::vector< std::vector<double> > E_1_L   ( n_points_r, std::vector<double> ( ell_max+1 ) );
std::vector< std::vector<double> > E_2_L   ( n_points_r, std::vector<double> ( ell_max+1 ) );
std::vector< std::vector<double> > ExB_r_L ( n_points_r, std::vector<double> ( ell_max+1 ) );
std::vector< std::vector<double> > ExB_1_L ( n_points_r, std::vector<double> ( ell_max+1 ) );
std::vector< std::vector<double> > ExB_2_L ( n_points_r, std::vector<double> ( ell_max+1 ) );

std::vector< std::vector<double> > B_r_L_dr  ( n_points_r, std::vector<double> ( ell_max+1 ) );
std::vector< std::vector<double> > B_1_L_dr  ( n_points_r, std::vector<double> ( ell_max+1 ) );
std::vector< std::vector<double> > B_2_L_dr  ( n_points_r, std::vector<double> ( ell_max+1 ) );
std::vector< std::vector<double> > E_r_L_dr  ( n_points_r, std::vector<double> ( ell_max+1 ) );
std::vector< std::vector<double> > E_1_L_dr  ( n_points_r, std::vector<double> ( ell_max+1 ) );
std::vector< std::vector<double> > E_2_L_dr  ( n_points_r, std::vector<double> ( ell_max+1 ) );
std::vector< std::vector<double> > B_r_L_dr2 ( n_points_r, std::vector<double> ( ell_max+1 ) );

std::vector< std::vector<double> > B_r ( n_points_r, std::vector<double> ( n_points_t ) );
std::vector< std::vector<double> > B_t ( n_points_r, std::vector<double> ( n_points_t ) );
std::vector< std::vector<double> > B_p ( n_points_r, std::vector<double> ( n_points_t ) );
std::vector< std::vector<double> > E_r ( n_points_r, std::vector<double> ( n_points_t ) );
std::vector< std::vector<double> > E_t ( n_points_r, std::vector<double> ( n_points_t ) );
std::vector< std::vector<double> > E_p ( n_points_r, std::vector<double> ( n_points_t ) );
std::vector< std::vector<double> > B_x ( n_points_r, std::vector<double> ( n_points_t ) );
std::vector< std::vector<double> > B_y ( n_points_r, std::vector<double> ( n_points_t ) );
std::vector< std::vector<double> > B_z ( n_points_r, std::vector<double> ( n_points_t ) );
std::vector< std::vector<double> > E_x ( n_points_r, std::vector<double> ( n_points_t ) );
std::vector< std::vector<double> > E_y ( n_points_r, std::vector<double> ( n_points_t ) );
std::vector< std::vector<double> > E_z ( n_points_r, std::vector<double> ( n_points_t ) );

std::vector< std::vector<double> > E_dot_B      ( n_points_r, std::vector<double> ( n_points_t ) );		// To check that force-free conditions are fulfilled.
std::vector< std::vector<double> > Bsq_minus_Esq( n_points_r, std::vector<double> ( n_points_t ) );		// To check that force-free conditions are fulfilled.

std::vector< std::vector<double> > alpha             ( n_points_r, std::vector<double> ( n_points_t ) );
std::vector< std::vector<double> > beta              ( n_points_r, std::vector<double> ( n_points_t ) );
std::vector< std::vector<double> > E_dot_B_over_Bsq  ( n_points_r, std::vector<double> ( n_points_t ) );

std::vector<double> alpha_integrand_0 ( n_points_t );
std::vector<double> alpha_integrand_1 ( n_points_t );
std::vector<double> beta_integrand_0  ( n_points_t );
std::vector<double> beta_integrand_1  ( n_points_t );
std::vector<double> E_VSH_integrand_r ( n_points_t );
std::vector<double> E_VSH_integrand_1 ( n_points_t );
std::vector<double> E_VSH_integrand_2 ( n_points_t );
std::vector< std::vector<double> > alpha_integrals_0( ell_max+1, std::vector<double> ( ell_max+1 ) );
std::vector< std::vector<double> > alpha_integrals_1( ell_max+1, std::vector<double> ( ell_max+1 ) );
std::vector< std::vector<double> > beta_integrals_0 ( ell_max+1, std::vector<double> ( ell_max+1 ) );
std::vector< std::vector<double> > beta_integrals_1 ( ell_max+1, std::vector<double> ( ell_max+1 ) );

// DELETE
std::vector<double> beta_integrand ( n_points_t );
std::vector< std::vector< std::vector<double> > > beta_integrals( ell_max+1, std::vector< std::vector<double> > (ell_max+1, std::vector<double> ( ell_max+1 ) ) );
// end DELETE

std::vector<double> B_r_L_dT_i ( ell_max+1 );
std::vector<double> B_2_L_dT_i ( ell_max+1 );
std::vector<double> E_r_L_dT_i ( ell_max+1 );
std::vector<double> E_1_L_dT_i ( ell_max+1 );
std::vector<double> E_2_L_dT_i ( ell_max+1 );

std::vector<double> B_tilde_integrand_0     ( n_points_t );
std::vector<double> B_tilde_integrand_1     ( n_points_t );
std::vector<double> E_corrected_integrand_0 ( n_points_t );
std::vector<double> E_corrected_integrand_1 ( n_points_t );
std::vector< std::vector<double> > B_tilde_integrals_0    ( ell_max+1, std::vector<double> ( ell_max+1 ) );
std::vector< std::vector<double> > B_tilde_integrals_1    ( ell_max+1, std::vector<double> ( ell_max+1 ) );
std::vector< std::vector<double> > E_corrected_integrals_0( ell_max+1, std::vector<double> ( ell_max+1 ) );
std::vector< std::vector<double> > E_corrected_integrals_1( ell_max+1, std::vector<double> ( ell_max+1 ) );

std::vector< std::vector<double> > E_tilde_r_L     ( n_points_r, std::vector<double> ( ell_max+1 ) );
std::vector< std::vector<double> > E_tilde_1_L     ( n_points_r, std::vector<double> ( ell_max+1 ) );
std::vector< std::vector<double> > E_tilde_2_L     ( n_points_r, std::vector<double> ( ell_max+1 ) );
std::vector< std::vector<double> > E_corrected_r_L ( n_points_r, std::vector<double> ( ell_max+1 ) );
std::vector< std::vector<double> > E_corrected_1_L ( n_points_r, std::vector<double> ( ell_max+1 ) );
std::vector< std::vector<double> > E_corrected_2_L ( n_points_r, std::vector<double> ( ell_max+1 ) );

std::ofstream output_file;

double time_start_seconds = 0;

double T = 0;




//----- Functions -----

std::vector<double> evaluate_VSH_series_and_FD( std::vector< std::vector<double> > A_r, std::vector< std::vector<double> > A_1, std::vector< std::vector<double> > A_2, int i, int j ){
	// Evaluate vector given its VSH series coefficients, which depend on radial position, at a gridpoint (i,j).
	// First three elements are spherical polar; last three elements are Cartesian.
	std::vector<double> A ( 6 );
	
	A[0] = pow( 4.0 * pi, -0.5 ) * A_r[i][0];
	
	for( int L=1; L<=ell_max; L++ ){
		double sqrt_2Lplus1_over_4pi = sqrt( (2.0*L+1.0)/(4.0*pi) );
		A[0] += sqrt_2Lplus1_over_4pi * A_r[i][L] * P0[j][L];
		A[1] += sqrt_2Lplus1_over_4pi * A_1[i][L] * P1[j][L];
		A[2] += sqrt_2Lplus1_over_4pi * A_2[i][L] * P1[j][L];
	}
	
	A[3] = A[0] * st[j] + A[1] * ct[j];
	A[4] = A[2];									// Extra variable assignment not necessary, but use for now to avoid confusion.
	A[5] = A[0] * ct[j] - A[1] * st[j];
	
	return A;
}

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




double wigner_3j( int L1, int L2, int L3, int M1, int M2, int M3 ){
	// https://en.wikipedia.org/wiki/3-j_symbol#Explicit_expression
	
	bool condition_M = ( M1+M2+M3 != 0 );
	bool condition_triangle = ( ( L1+L2-L3 +1 < 0 ) or ( L1-L2+L3 < 0 ) or ( -L1+L2+L3 < 0 ) or ( L1+L2+L3+1 < 0 ) );
	bool condition_LM = ( ( L1-M1 < 0 ) or ( L1+M1 < 0 ) or ( L2-M2 < 0 ) or ( L2+M2 < 0 ) or ( L3-M3 < 0 ) or ( L3+M3 < 0 ) );
	
	if( condition_M or condition_triangle or condition_LM ){
		return 0;
	}
	
	
	double ret = 0;
	
	int k_min = std::max( std::max( 0, -L3+L2-M1 ), -L3+L1+M2 );
	int k_max = std::min( std::min( L1+L2-L3, L1-M1 ), L2+M2 );
	
	for( int k=k_min; k<=k_max; k++ ){
		ret += pow(-1,k) / ( tgamma( k +1 ) * tgamma( L1+L2-L3-k +1 ) * tgamma( L1-M1-k +1 ) * tgamma( L2+M2-k +1 ) * tgamma( L3-L2+M1+k +1 ) * tgamma( L3-L1-M2+k +1 ) );
	}
	
	ret *= sqrt( tgamma( L1-M1 +1 ) * tgamma( L1+M1 +1 ) * tgamma( L2-M2 +1 ) * tgamma( L2+M2 +1 ) * tgamma( L3-M3 +1 ) * tgamma( L3+M3 + 1 ) );
	ret *= sqrt( tgamma( L1+L2-L3 +1 ) * tgamma( L1-L2+L3 +1 ) * tgamma( -L1+L2+L3 +1 ) / tgamma( L1+L2+L3+1 + 1 ) );
	ret *= pow( -1, L1-L2-M3 );
	
	return ret;
}

double wigner_3j_000( int L1, int L2, int L3 ){
	// https://en.wikipedia.org/wiki/3-j_symbol#Explicit_expression
	
	bool condition_triangle = ( ( L1+L2-L3 +1 < 0 ) or ( L1-L2+L3 < 0 ) or ( -L1+L2+L3 < 0 ) or ( L1+L2+L3+1 < 0 ) );
	bool condition_L = ( ( L1 < 0 ) or ( L2 < 0 ) or ( L3 < 0 ) );
	
	if( condition_triangle or condition_L ){
		return 0;
	}
	
	
	double ret = 0;
	
	int k_min = std::max( std::max( 0, L2-L3 ), L1-L3 );
	int k_max = std::min( std::min( L1+L2-L3, L1 ), L2 );
	
	for( int k=k_min; k<=k_max; k++ ){
		ret += pow(-1,k) / ( tgamma( k +1 ) * tgamma( L1+L2-L3-k +1 ) * tgamma( L1-k +1 ) * tgamma( L2-k +1 ) * tgamma( L3-L2+k +1 ) * tgamma( L3-L1+k +1 ) );
	}
	
	ret *= tgamma( L1 +1 ) * tgamma( L2 +1 ) * tgamma( L3 +1 );
	ret *= sqrt( tgamma( L1+L2-L3 +1 ) * tgamma( L1-L2+L3 +1 ) * tgamma( -L1+L2+L3 +1 ) / tgamma( L1+L2+L3+1 + 1 ) );
	ret *= pow( -1, L1-L2 );
	
	return ret;
}




double I_011( int L1, int L2, int L3 ){
	// 1/(4sqrt(pi)) sqrt(2L1+1) sqrt(L2(L2+1)(2L2+1)) sqrt(L3(L3+1)(2L3+1)) int_0^pi P_{L1}^0 P_{L2}^1 P_{L3}^1 sin(theta) dtheta
	
	if( L3 > L1+L2 ){
		return 0;
	}
	
	double ret = 0;
	
	for( int L12=std::max(abs(L1-L2),1); L12<=L1+L2; L12++ ){
		if( (L1+L2+L12)%2 == 0 ){
		
			double G12 = - ( 2*L12+1 ) * wigner_3j_000( L1, L2, L12 ) * wigner_3j( L1, L2, L12, 0, 1, -1 );
			
			if( (L12+L3)%2 == 0 ){
				for( int L123=std::max(abs(L12-L3),2); L123<=L12+L3; L123++ ){
					if( L123 %2 == 0 ){
						double G123          = ( 2*L123+1 ) * wigner_3j_000( L12, L3, L123 ) * wigner_3j( L12, L3, L123, 1, 1, -2 );
						double sqrt_fraction = pow( (L123+2) * (L123+1) * L123 * (L123-1), -0.5 );
						
						ret += G12 * G123 * sqrt_fraction;
					}
				}
			}
			
		}
	}
	
	return ret * sqrt( (2*L1+1) * L2*(L2+1)*(2*L2+1) * L3*(L3+1)*(2*L3+1) / pi );
}




void calculate_AVSH_series_coeffs_of_VSH_cross_products(){
	
	for( int L1=0; L1<=ell_max; L1++ ){
		for( int L2=0; L2<=ell_max; L2++ ){
			for( int L=1; L<=ell_max; L++ ){
				if ( L <= L1 + L2 ){
					Y_L1_cross_Psi_L2_coeff_2_L  [L1][L2][L] = I_011( L1, L2, L ) / ( (double) L*(L+1) );
					Y_L1_cross_Phi_L2_coeff_1_L  [L1][L2][L] = - Y_L1_cross_Psi_L2_coeff_2_L[L1][L2][L];
					Psi_L1_cross_Phi_L2_coeff_r_L[L1][L2][L] = I_011( L, L1, L2 );
				}
			}
		}
	}
	
}




void apply_initial_VSH_coeff_values(){
	
	for( int i=0; i<n_points_r; i++ ){
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




void calculate_radial_derivatives_of_VSH_coeffs_inner_boundary(){
	// Inner BCs mean fields don't change there, so this only needs to be called once. Hence, move into separate function outside loop over time.
	
	for( int L=0; L<=ell_max; L++ ){		
			E_r_L_dr [0][L] = ( E_r_L[1][L] - E_r_L[0][L] ) / ( 2.0*delta_r );
			E_1_L_dr [0][L] = ( E_1_L[1][L] - E_1_L[0][L] ) / ( 2.0*delta_r );
			E_2_L_dr [0][L] = ( E_2_L[1][L] - E_2_L[0][L] ) / ( 2.0*delta_r );
			B_r_L_dr [0][L] = ( B_r_L[1][L] - B_r_L[0][L] ) / ( 2.0*delta_r );
			B_2_L_dr [0][L] = ( B_2_L[1][L] - B_2_L[0][L] ) / ( 2.0*delta_r );
			
			B_r_L_dr2[0][L] = ( B_r_L[1][L] - B_r_L[0][L] ) / ( delta_r*delta_r );
	}
}

void calculate_radial_derivatives_of_VSH_coeffs(){
	// Use one-sided derivatives at the endpoints and symmetric derivatives for all other points.
	// Inner BCs mean fields don't change there, so no need to recalculate at i=0.
	
	for( int L=0; L<=ell_max; L++ ){
		
		for( int i=1; i<n_points_r-1; i++ ){
			E_r_L_dr [i][L] = ( E_r_L[i+1][L] - E_r_L[i-1][L] ) / ( 2.0*delta_r );
			E_1_L_dr [i][L] = ( E_1_L[i+1][L] - E_1_L[i-1][L] ) / ( 2.0*delta_r );
			E_2_L_dr [i][L] = ( E_2_L[i+1][L] - E_2_L[i-1][L] ) / ( 2.0*delta_r );
			B_r_L_dr [i][L] = ( B_r_L[i+1][L] - B_r_L[i-1][L] ) / ( 2.0*delta_r );
			B_2_L_dr [i][L] = ( B_2_L[i+1][L] - B_2_L[i-1][L] ) / ( 2.0*delta_r );
			
			B_r_L_dr2[i][L] = ( B_r_L[i+1][L] - 2.0*B_r_L[i][L] + B_r_L[i-1][L] ) / ( delta_r*delta_r );
		}
		
		E_r_L_dr .back()[L] = ( E_r_L.back()[L] - E_r_L[n_points_r-2][L] ) / ( 2.0*delta_r );
		E_1_L_dr .back()[L] = ( E_1_L.back()[L] - E_1_L[n_points_r-2][L] ) / ( 2.0*delta_r );
		E_2_L_dr .back()[L] = ( E_2_L.back()[L] - E_2_L[n_points_r-2][L] ) / ( 2.0*delta_r );
		B_r_L_dr .back()[L] = ( B_r_L.back()[L] - B_r_L[n_points_r-2][L] ) / ( 2.0*delta_r );
		B_2_L_dr .back()[L] = ( B_2_L.back()[L] - B_2_L[n_points_r-2][L] ) / ( 2.0*delta_r );
		
		B_r_L_dr2.back()[L] = ( -B_r_L.back()[L] + B_r_L[n_points_r-2][L] ) / ( delta_r*delta_r );
		
	}

}




void calculate_B_1_L_and_B_1_L_dr_inner_boundary( bool change_B_1_L = true ){
	// Inner BCs mean fields don't change there, so this only needs to be called once. Hence, move into separate function outside loop over time.
	for( int L=1; L<=ell_max; L++ ){
		if( change_B_1_L ) { B_1_L   [0][L] = ( r[0] * B_r_L_dr [0][L] + 2.0 * B_r_L[0][L] ) / ( L*(L+1.0) ); }
		B_1_L_dr[0][L] = ( r[0] * B_r_L_dr2[0][L] + 3.0 * B_r_L[0][L] ) / ( L*(L+1.0) );
	}
	
}


void calculate_B_1_L_and_B_1_L_dr( bool change_B_1_L = true ){
	
	for( int L=1; L<=ell_max; L++ ){
		for( int i=1; i<n_points_r; i++ ){
			if( change_B_1_L ) { B_1_L[i][L] = ( r[i] * B_r_L_dr[i][L] + 2.0 * B_r_L[i][L] ) / ( L*(L+1.0) ); }
			B_1_L_dr[i][L] = ( r[i] * B_r_L_dr2[i][L] + 3.0 * B_r_L[i][L] ) / ( L*(L+1.0) );
		}
	}
	
}




void calculate_alpha_and_beta(){
	// Along the way, we need B, E, div(E), curl(E) and curl(B), so calculate these.
	// Keep the spherical polar components stored in a vector because we will want them later when calculating Cartesian components for CSV output.
	
	B_r = std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) );
	B_t = std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) );
	B_p = std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) );
	E_r = std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) );
	E_t = std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) );
	E_p = std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) );
	
	for( int i=0; i<n_points_r; i++ ){
		for( int j=0; j<n_points_t; j++ ){
	
			double div_E = 0;
			double curl_E_r = 0;
			double curl_E_t = 0;
			double curl_E_p = 0;
			double curl_B_r = 0;
			double curl_B_t = 0;
			double curl_B_p = 0;
			
			
			//--- L=0 ---
			double sqrt_2Lplus1_over_4pi = 1.0/(2.0*sqrt(pi));
			div_E += sqrt_2Lplus1_over_4pi * ( E_r_L_dr[i][0] + 2.0*pow(r[i],-1) * E_r_L[i][0] );
			E_r[i][j] += sqrt_2Lplus1_over_4pi * E_r_L[i][0];
			B_r[i][j] += sqrt_2Lplus1_over_4pi * B_r_L[i][0];
			
			//--- L>0 ---
			for( int L=1; L<=ell_max; L++ ){
				double sqrt_2Lplus1_over_4pi = sqrt( (2.0*L+1.0)/(4.0*pi) );
				
				E_r[i][j] += sqrt_2Lplus1_over_4pi * E_r_L[i][L] * P0[j][L];
				E_t[i][j] += sqrt_2Lplus1_over_4pi * E_1_L[i][L] * P1[j][L];
				E_p[i][j] += sqrt_2Lplus1_over_4pi * E_2_L[i][L] * P1[j][L];
				
				B_r[i][j] += sqrt_2Lplus1_over_4pi * B_r_L[i][L] * P0[j][L];
				B_t[i][j] += sqrt_2Lplus1_over_4pi * B_1_L[i][L] * P1[j][L];
				B_p[i][j] += sqrt_2Lplus1_over_4pi * B_2_L[i][L] * P1[j][L];
				
				div_E += sqrt_2Lplus1_over_4pi * ( E_r_L_dr[i][L] + 2.0*pow(r[i],-1) * E_r_L[i][L] - L*(L+1.0)*pow(r[i],-1) * E_1_L[i][L] );
				
				curl_E_r += sqrt_2Lplus1_over_4pi * P0[j][L] * -L*(L+1) * E_2_L[i][L] / r[i];
				curl_E_t += sqrt_2Lplus1_over_4pi * P1[j][L] * ( -E_2_L_dr[i][L] - E_2_L[i][L] / r[i] );
				curl_E_p += sqrt_2Lplus1_over_4pi * P1[j][L] * ( E_1_L_dr[i][L] + ( E_1_L[i][L] - E_r_L[i][L] ) / r[i] );
				
				curl_B_r += sqrt_2Lplus1_over_4pi * P0[j][L] * -L*(L+1) * B_2_L[i][L] / r[i];
				curl_B_t += sqrt_2Lplus1_over_4pi * P1[j][L] * ( -B_2_L_dr[i][L] - B_2_L[i][L] / r[i] );
				curl_B_p += sqrt_2Lplus1_over_4pi * P1[j][L] * ( r[i] * B_r_L_dr2[i][L] + 4.0 * B_r_L_dr[i][L] - (L+2.0)*(L-1.0) * B_r_L[i][L] / r[i] ) / ( L*(L+1.0) );
			}
			
			double B_dot_curl_B = B_r[i][j]*curl_B_r  + B_t[i][j]*curl_B_t  + B_p[i][j]*curl_B_p;
			double E_dot_curl_E = E_r[i][j]*curl_E_r  + E_t[i][j]*curl_E_t  + E_p[i][j]*curl_E_p;
			double B_squared    = B_r[i][j]*B_r[i][j] + B_t[i][j]*B_t[i][j] + B_p[i][j]*B_p[i][j];
			double E_squared    = E_r[i][j]*E_r[i][j] + E_t[i][j]*E_t[i][j] + E_p[i][j]*E_p[i][j];
			
			alpha[i][j] = ( E_dot_curl_E - B_dot_curl_B ) / B_squared;
			beta [i][j] = - div_E / B_squared;
			
			E_dot_B      [i][j] = E_r[i][j]*B_r[i][j] + E_t[i][j]*B_t[i][j] + E_p[i][j]*B_p[i][j];
			Bsq_minus_Esq[i][j] = B_squared - E_squared;
			
			E_dot_B_over_Bsq[i][j] = E_dot_B[i][j] / B_squared;
			
		}
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



void calculate_alpha_and_beta_integrals( int i ){
	
	for( int L1=0; L1<=ell_max; L1++ ){
		
		//--- Diagonal terms ---
		for( int j=0; j<n_points_t; j++ ){
			alpha_integrand_0[j] = alpha[i][j] * P0[j][L1] * P0[j][L1] * st[j];
			alpha_integrand_1[j] = alpha[i][j] * P1[j][L1] * P1[j][L1] * st[j];
			beta_integrand_0 [j] = beta [i][j] * P0[j][L1] * P0[j][L1] * st[j];
			beta_integrand_1 [j] = beta [i][j] * P1[j][L1] * P1[j][L1] * st[j];
		}
		alpha_integrals_0[L1][L1] = integral_trapezium( alpha_integrand_0, t );
		alpha_integrals_0[L1][L1] = integral_trapezium( alpha_integrand_1, t );
		beta_integrals_0 [L1][L1] = integral_trapezium( beta_integrand_0 , t );
		beta_integrals_0 [L1][L1] = integral_trapezium( beta_integrand_1 , t );
		
		//--- Upper- and lower-triangular terms ---
		for( int L2=L1+1; L2<=ell_max; L2++ ){
			for( int j=0; j<n_points_t; j++ ){
				alpha_integrand_0[j] = alpha[i][j] * P0[j][L1] * P0[j][L2] * st[j];
				alpha_integrand_1[j] = alpha[i][j] * P1[j][L1] * P1[j][L2] * st[j];
				beta_integrand_0 [j] = beta [i][j] * P0[j][L1] * P0[j][L2] * st[j];
				beta_integrand_1 [j] = beta [i][j] * P1[j][L1] * P1[j][L2] * st[j];
			}
			alpha_integrals_0[L1][L2] = integral_trapezium( alpha_integrand_0, t );
			alpha_integrals_0[L1][L2] = integral_trapezium( alpha_integrand_1, t );
			beta_integrals_0 [L1][L2] = integral_trapezium( beta_integrand_0 , t );
			beta_integrals_0 [L1][L2] = integral_trapezium( beta_integrand_1 , t );
			alpha_integrals_0[L2][L1] = alpha_integrals_0[L1][L2];
			alpha_integrals_1[L2][L1] = alpha_integrals_1[L1][L2];
			beta_integrals_0 [L2][L1] = beta_integrals_0 [L1][L2];
			beta_integrals_1 [L2][L1] = beta_integrals_1 [L1][L2];
		}
	}
}
			



void calculate_beta_integrals( int i ){
	// DELETE
	// In this faster method, we can only swap the [L1][L2] indices because [L3] corresponds to m=0.
	// So, we move the L3 loop to the first loop, and apply the same algorithm as for alpha over L1,L2.
	
	for( int L3=0; L3<=ell_max; L3++ ){
		
		for( int L1=0; L1<=ell_max; L1++ ){
			
			//--- Diagonal terms ---
			for( int j=0; j<n_points_t; j++ ){
				beta_integrand[j] = beta[i][j] * P1[j][L1] * P1[j][L1] * P0[j][L3] * st[j];
			}				
			beta_integrals[L1][L1][L3] = integral_trapezium( beta_integrand, t );
			
			//--- Upper- and lower-triangular terms ---
			for( int L2=L1+1; L2<=ell_max; L2++ ){
				for( int j=0; j<n_points_t; j++ ){
					beta_integrand[j] = beta[i][j] * P1[j][L1] * P1[j][L2] * P0[j][L3] * st[j];
				}
				beta_integrals[L1][L2][L3] = integral_trapezium( beta_integrand, t );
				beta_integrals[L2][L1][L3] = beta_integrals[L1][L2][L3];
			}
		}
	}
}




void calculate_ExB_VSH_coeffs( int i ){
	
	for( int L=0; L<=ell_max; L++ ){
		for( int L1=0; L1<=ell_max; L1++ ){
			for( int L2=1; L2<=ell_max; L2++ ){
				ExB_r_L[i][L] += ( E_1_L[i][L1]*B_2_L[i][L2] - E_2_L[i][L2]*B_1_L[i][L1] ) * Psi_L1_cross_Phi_L2_coeff_r_L[L1][L2][L];
				ExB_1_L[i][L] += ( E_r_L[i][L1]*B_2_L[i][L2] - E_2_L[i][L2]*B_r_L[i][L1] ) * Y_L1_cross_Phi_L2_coeff_1_L  [L1][L2][L];
				ExB_2_L[i][L] += ( E_r_L[i][L1]*B_1_L[i][L2] - E_1_L[i][L2]*B_r_L[i][L1] ) * Y_L1_cross_Psi_L2_coeff_2_L  [L1][L2][L];
			}
		}
	}
	
}
	




void calculate_time_derivatives_of_VSH_coeffs( int i ){
	
	//--- Loop inc. L=0 for radial term ---
	for( int L=0; L<=ell_max; L++ ){
			
		B_r_L_dT_i[L] =  L*(L+1.0) * E_2_L[i][L] / r[i];
		E_r_L_dT_i[L] = -L*(L+1.0) * B_2_L[i][L] / r[i];
		
		double E_r_L_dT_term = 0;
		
		for( int L1=0; L1<=ell_max; L1++ ){
			E_r_L_dT_term += sqrt( 2.0*L1+1.0 ) * ( B_r_L[i][L1] * alpha_integrals_0[L1][L] + ExB_r_L[i][L1] * beta_integrals_0[L1][L] );
		}
		
		E_r_L_dT_term *= 0.5 * sqrt( 2.0*L+1.0 );
		
		E_r_L_dT_i[L] += E_r_L_dT_term;
	
	}
	
	//--- Loop from L=1 for angular terms ---
	for( int L=1; L<=ell_max; L++ ){
			
		B_2_L_dT_i[L] = -E_1_L_dr[i][L] + ( E_r_L[i][L] - E_1_L[i][L] ) / r[i];
		
		E_1_L_dT_i[L] = -B_2_L_dr[i][L] - B_2_L[i][L] / r[i];
		E_2_L_dT_i[L] = B_1_L_dr[i][L] + ( B_1_L[i][L] - B_r_L[i][L] ) / r[i];
		
		
		double E_1_L_dT_term = 0;
		double E_2_L_dT_term = 0;
		
		for( int L1=1; L1<=ell_max; L1++ ){
			E_1_L_dT_term += sqrt( 2.0*L1+1.0 ) * ( B_1_L[i][L1] * alpha_integrals_1[L1][L] + ExB_1_L[i][L1] * beta_integrals_1[L1][L] );
			E_2_L_dT_term += sqrt( 2.0*L1+1.0 ) * ( B_2_L[i][L1] * alpha_integrals_1[L1][L] + ExB_2_L[i][L1] * beta_integrals_1[L1][L] );
		}
		
		E_1_L_dT_term *= sqrt( 2.0*L+1.0 ) / ( 2.0 * L*(L+1.0) );
		E_2_L_dT_term *= sqrt( 2.0*L+1.0 ) / ( 2.0 * L*(L+1.0) );
		
		E_1_L_dT_i[L] += E_1_L_dT_term;
		E_2_L_dT_i[L] += E_2_L_dT_term;
		
	}
	
}




void integrate_VSH_coeffs_wrt_time( int i ){
	for( int L=0; L<=ell_max; L++ ){
		B_r_L[i][L] += B_r_L_dT_i[L] * delta_T;
		B_2_L[i][L] += B_2_L_dT_i[L] * delta_T;
		E_r_L[i][L] += E_r_L_dT_i[L] * delta_T;
		E_1_L[i][L] += E_1_L_dT_i[L] * delta_T;
		E_2_L[i][L] += E_2_L_dT_i[L] * delta_T;
	}
}




void calculate_B_tilde_integrals( int i ){
	// NOT UPDATED 20230810
	// Copy of calculate_alpha_integrals(int i).
	// Could make generalised function with alpha and output_arrays as arguments, but no real point. For now, I only need these three applications.
	
	for( int L1=0; L1<=ell_max; L1++ ){
		
		//--- Diagonal terms ---
		for( int j=0; j<n_points_t; j++ ){
			B_tilde_integrand_0[j] = E_dot_B_over_Bsq[i][j] * P0[j][L1] * P0[j][L1] * st[j];
			B_tilde_integrand_1[j] = E_dot_B_over_Bsq[i][j] * P1[j][L1] * P1[j][L1] * st[j];
		}
		B_tilde_integrals_0[L1][L1] = integral_trapezium( B_tilde_integrand_0, t );
		B_tilde_integrals_0[L1][L1] = integral_trapezium( B_tilde_integrand_1, t );
		
		//--- Upper- and lower-triangular terms ---
		for( int L2=L1+1; L2<=ell_max; L2++ ){
			for( int j=0; j<n_points_t; j++ ){
				B_tilde_integrand_0[j] = E_dot_B_over_Bsq[i][j] * P0[j][L1] * P0[j][L2] * st[j];
				B_tilde_integrand_1[j] = E_dot_B_over_Bsq[i][j] * P1[j][L1] * P1[j][L2] * st[j];
			}
			B_tilde_integrals_0[L1][L2] = integral_trapezium( B_tilde_integrand_0, t );
			B_tilde_integrals_0[L1][L2] = integral_trapezium( B_tilde_integrand_1, t );
			B_tilde_integrals_0[L2][L1] = B_tilde_integrals_0[L1][L2];
			B_tilde_integrals_1[L2][L1] = B_tilde_integrals_1[L1][L2];
		}
	}
}




void calculate_E_corrected_integrals( int i ){
	// NOT UPDATED 20230810
	// Copy of calculate_alpha_integrals(int i).
	// Could make generalised function with alpha and output_arrays as arguments, but no real point. For now, I only need these three applications.
	
	std::vector<double> sqrt_Bsq_over_E_tilde_sq( n_points_t );
	for( int j=0; j<n_points_t; j++ ){
		double Bsq = B_r[i][j]*B_r[i][j] + B_t[i][j]*B_t[i][j] + B_p[i][j]*B_p[i][j];
		std::vector<double> E_tilde = evaluate_VSH_series_and_FD( E_tilde_r_L, E_tilde_1_L, E_tilde_2_L, i, j );
		double E_tilde_sq = E_tilde[0]*E_tilde[0] + E_tilde[1]*E_tilde[1] + E_tilde[2]*E_tilde[2];
		sqrt_Bsq_over_E_tilde_sq[j] = sqrt( Bsq / E_tilde_sq );
	}
	
	for( int L1=0; L1<=ell_max; L1++ ){
		
		//--- Diagonal terms ---
		for( int j=0; j<n_points_t; j++ ){
			E_corrected_integrand_0[j] = sqrt_Bsq_over_E_tilde_sq[j] * P0[j][L1] * P0[j][L1] * st[j];
			E_corrected_integrand_1[j] = sqrt_Bsq_over_E_tilde_sq[j] * P1[j][L1] * P1[j][L1] * st[j];
		}
		E_corrected_integrals_0[L1][L1] = integral_trapezium( E_corrected_integrand_0, t );
		E_corrected_integrals_0[L1][L1] = integral_trapezium( E_corrected_integrand_1, t );
		
		//--- Upper- and lower-triangular terms ---
		for( int L2=L1+1; L2<=ell_max; L2++ ){
			for( int j=0; j<n_points_t; j++ ){
				E_corrected_integrand_0[j] = sqrt_Bsq_over_E_tilde_sq[j] * P0[j][L1] * P0[j][L2] * st[j];
				E_corrected_integrand_1[j] = sqrt_Bsq_over_E_tilde_sq[j] * P1[j][L1] * P1[j][L2] * st[j];
			}
			E_corrected_integrals_0[L1][L2] = integral_trapezium( B_tilde_integrand_0, t );
			E_corrected_integrals_0[L1][L2] = integral_trapezium( B_tilde_integrand_1, t );
			E_corrected_integrals_0[L2][L1] = E_corrected_integrals_0[L1][L2];
			E_corrected_integrals_1[L2][L1] = E_corrected_integrals_1[L1][L2];
		}
	}
}




void reset_arrays(){
	E_tilde_r_L     = std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( ell_max+1 ) );
	E_tilde_1_L     = std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( ell_max+1 ) );
	E_tilde_2_L     = std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( ell_max+1 ) );
	E_corrected_r_L = std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( ell_max+1 ) );
	E_corrected_1_L = std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( ell_max+1 ) );
	E_corrected_2_L = std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( ell_max+1 ) );
}
	




void calculate_E_tilde_coeffs( int i ){
	// NOT UPDATED 20230810
	
	/*
	for( int L=0; L<=ell_max; L++ ){
		
		double E_tilde_r_L_term = 0;
		double E_tilde_1_L_term = 0;
		double E_tilde_2_L_term = 0;
		
		for( int L1=0; L1<=ell_max; L1++ ){
			E_tilde_r_L_term += sqrt(2.0*L1+1.0) * B_r_L[i][L1] * B_tilde_integrals_0[L1][L];
		}
		
		for( int L1=1; L1<=ell_max; L1++ ){
			E_tilde_1_L_term += sqrt(2.0*L1+1.0) * B_1_L[i][L1] * B_tilde_integrals_1[L1][L];
			E_tilde_2_L_term += sqrt(2.0*L1+1.0) * B_2_L[i][L1] * B_tilde_integrals_1[L1][L];
		}
		
		
		E_tilde_r_L_term *= sqrt(2.0*L+1.0) / ( 2.0*pi );
		E_tilde_1_L_term *= sqrt(2.0*L+1.0) / ( (double)L*(L+1.0) * 2.0*pi );
		E_tilde_2_L_term *= sqrt(2.0*L+1.0) / ( (double)L*(L+1.0) * 2.0*pi );
		
		E_tilde_r_L[i][L] = E_r_L[i][L] - E_tilde_r_L_term;
		E_tilde_1_L[i][L] = E_1_L[i][L] - E_tilde_2_L_term;
		E_tilde_2_L[i][L] = E_2_L[i][L] - E_tilde_1_L_term;
	}
	*/
	
	for( int L=0; L<=ell_max; L++ ){
		
		for( int L1=0; L1<=ell_max; L1++ ){
			E_tilde_r_L[i][L] += sqrt(2.0*L1+1.0) * B_r_L[i][L1] * B_tilde_integrals_0[L1][L];
		}
		
		for( int L1=1; L1<=ell_max; L1++ ){
			E_tilde_1_L[i][L] += sqrt(2.0*L1+1.0) * B_1_L[i][L1] * B_tilde_integrals_1[L1][L];
			E_tilde_2_L[i][L] += sqrt(2.0*L1+1.0) * B_2_L[i][L1] * B_tilde_integrals_1[L1][L];
		}
		
		
		E_tilde_r_L[i][L] *= - sqrt(2.0*L+1.0) / ( 2.0*pi );
		E_tilde_1_L[i][L] *= - sqrt(2.0*L+1.0) / ( (double)L*(L+1.0) * 2.0*pi );
		E_tilde_2_L[i][L] *= - sqrt(2.0*L+1.0) / ( (double)L*(L+1.0) * 2.0*pi );
		
		E_tilde_r_L[i][L] += E_r_L[i][L];
		E_tilde_1_L[i][L] += E_1_L[i][L];
		E_tilde_2_L[i][L] += E_2_L[i][L];
		
	}
	
}




void calculate_E_corrected_coeffs( int i ){
	// NOT UPDATED 20230810
	
	for( int L=0; L<=ell_max; L++ ){
		
		for( int L1=0; L1<=ell_max; L1++ ){
			E_corrected_r_L[i][L] += sqrt(2.0*L1+1.0) * E_tilde_r_L[i][L1] * E_corrected_integrals_0[L1][L];
		}
		
		for( int L1=1; L1<=ell_max; L1++ ){
			E_corrected_1_L[i][L] += sqrt(2.0*L1+1.0) * E_tilde_1_L[i][L1] * E_corrected_integrals_1[L1][L];
			E_corrected_2_L[i][L] += sqrt(2.0*L1+1.0) * E_tilde_2_L[i][L1] * E_corrected_integrals_1[L1][L];
		}
		
		E_corrected_r_L[i][L] *= sqrt(2.0*L+1.0) / ( 2.0*pi );
		E_corrected_1_L[i][L] *= sqrt(2.0*L+1.0) / ( (double)L*(L+1.0) * 2.0*pi );
		E_corrected_2_L[i][L] *= sqrt(2.0*L+1.0) / ( (double)L*(L+1.0) * 2.0*pi );
	}
}




void update_E_within_light_cylinder_if_greater_than_Bsq_and_redo_VSH_decomposition(){
	bool E_updated = false;
	int  points_updated = 0;
	
	for( int i=0; i<n_points_r; i++ ){
		for( int j=0; j<n_points_t; j++ ){
			if( ( Bsq_minus_Esq[i][j] < 0 ) and ( r[i]*st[j]<R_LC ) ){
				E_updated = true;
				points_updated += 1;
				double E_dot_B_over_Bsq = ( E_r[i][j]*B_r[i][j]+E_t[i][j]*B_t[i][j]+E_p[i][j]*B_p[i][j] ) / ( B_r[i][j]*B_r[i][j]+B_t[i][j]*B_t[i][j]+B_p[i][j]*B_p[i][j] );
				std::vector<double> E_prime { E_r[i][j] - E_dot_B_over_Bsq * B_r[i][j], E_t[i][j] - E_dot_B_over_Bsq * B_t[i][j], E_p[i][j] - E_dot_B_over_Bsq * B_p[i][j] };
				double sqrt_Bsq_over_Eprimesq = sqrt( ( B_r[i][j]*B_r[i][j]+B_t[i][j]*B_t[i][j]+B_p[i][j]*B_p[i][j] ) / ( E_prime[0]*E_prime[0] + E_prime[1]*E_prime[1] + E_prime[1]*E_prime[1] ) );
				std::vector<double> E_cor { E_prime[0] * sqrt_Bsq_over_Eprimesq, E_prime[1] * sqrt_Bsq_over_Eprimesq,  E_prime[2] * sqrt_Bsq_over_Eprimesq };
				E_r[i][j] = E_cor[0];
				E_t[i][j] = E_cor[1];
				E_p[i][j] = E_cor[2];
			}
		}
	}
	
	if( E_updated ){
		std::cout << "Updating E coeffs:\t" << points_updated << std::endl;
		for( int i=0; i<n_points_r; i++ ){
			for( int ell=0; ell<=ell_max; ell++ ){
			
				//--- Build lists of integrand values ---
				for( int j=0; j<n_points_t; j++ ){
					
					E_VSH_integrand_r[j] = E_r[i][j] * P0[j][ell] * st[j];
					E_VSH_integrand_1[j] = E_t[i][j] * P1[j][ell] * st[j];
					E_VSH_integrand_2[j] = E_p[i][j] * P1[j][ell] * st[j];
					
				}
				
				//--- Integrate and put values into array ---
				E_r_L[i][ell] = integral_trapezium( E_VSH_integrand_r, t ) * sqrt( (2.0*ell+1.0)*pi );
				E_1_L[i][ell] = integral_trapezium( E_VSH_integrand_1, t ) * sqrt( (2.0*ell+1.0)*pi ) / ( (double)ell*(ell+1.0) );
				E_2_L[i][ell] = integral_trapezium( E_VSH_integrand_2, t ) * sqrt( (2.0*ell+1.0)*pi ) / ( (double)ell*(ell+1.0) );
				
				
				if( ell == 0 ){
					E_1_L[i][ell] = 0;
					E_2_L[i][ell] = 0;
				}
				
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




void output_headers_to_csv(){
	// Split over separate lines for readability and ease of checking consistency with output_to_csv().
	
	//output_file << "T_index,T,T_sec,i,j,r,t,x,z,B_x,B_y,B_z,B,div(B),curl(B)_x,curl(B)_y,curl(B)_z,E_x,E_y,E_z,E,div(E),curl(E)_x,curl(E)_y,curl(E)_z,alpha,beta,E dot B,E minus B\n";
	
	output_file << "T_index,T,T_sec,i,j,r,t,x,z";
	output_file << ",B_r,B_t,B_p";
	output_file << ",B_x,B_y,B_z,B,E_x,E_y,E_z,E";
	output_file << ",alpha,beta";
	output_file << ",Etilde_x,Etilde_y,Etilde_z,Etx_check,Ety_check,Etz_check";
	output_file << ",E_dot_B,Bsq_minus_Esq,sgn(Bsq_minus_Esq)";
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
				
				bool condition_r = ( i % csv_write_freq_r == 0 ) or ( i == n_points_r - 1 );
				bool condition_t = ( j % csv_write_freq_t == 0 ) or ( j == n_points_t - 1 );
		
				if( condition_r and condition_t ){
					
					double B = sqrt( pow( B_r[i][j], 2 ) + pow( B_t[i][j], 2 ) + pow( B_p[i][j], 2 ) );
					double E = sqrt( pow( E_r[i][j], 2 ) + pow( E_t[i][j], 2 ) + pow( E_p[i][j], 2 ) );
					
					output_file << T_index <<","<< T <<","<< T*T_factor <<","<< i <<","<< j <<","<< r[i] <<","<< t[j] <<","<< x[i][j] <<","<< z[i][j]
					            <<","<< B_r[i][j] <<","<< B_t[i][j] <<","<< B_p[i][j]
								<<","<< B_x[i][j] <<","<< B_y[i][j] <<","<< B_z[i][j] <<","<< B
								<<","<< E_x[i][j] <<","<< E_y[i][j] <<","<< E_z[i][j] <<","<< E
								<<","<< alpha[i][j] <<","<< beta[i][j];
					
					//--- Evaluate E_tilde to verify that the application of the force-free conditions is done correctly ----
					// This is in spherical coordinates while E_corrected below is in Cartesian coordinates. Initially didn't trust.
					/*
					std::vector<double> E_tilde = evaluate_VSH_series_and_FD( E_tilde_r_L, E_tilde_1_L, E_tilde_2_L, i, j );
					double E_dot_B_over_Bsq = ( E_x[i][j]*B_x[i][j] + E_y[i][j]*B_y[i][j] + E_z[i][j]*B_z[i][j] ) / ( B_x[i][j]*B_x[i][j] + B_y[i][j]*B_y[i][j] + B_z[i][j]*B_z[i][j] );
					std::vector<double> Etildecheck { E_r[i][j] - E_dot_B_over_Bsq * B_r[i][j], E_t[i][j] - E_dot_B_over_Bsq * B_t[i][j], E_p[i][j] - E_dot_B_over_Bsq * B_p[i][j] };
					std::vector<double> Etildecheck_cart { Etildecheck[0] * st[j] + Etildecheck[1] * ct[j], Etildecheck[2], Etildecheck[0] * ct[j] - Etildecheck[1] * st[j] };
					
					output_file <<","<< E_tilde[3] <<","<< E_tilde[4] <<","<< E_tilde[5];
					output_file <<","<< Etildecheck_cart[0] <<","<< Etildecheck_cart[1] <<","<< Etildecheck_cart[2];
					output_file <<","<< abs(E_tilde[3]-Etildecheck_cart[0]) <<","<< abs(E_tilde[4]-Etildecheck_cart[1]) <<","<< abs(E_tilde[5]-Etildecheck_cart[2]) <<"\n";
					*/
					
					//--- Evaluate E_corrected to verify that the application of the force-free conditions is done correctly ----
					
					std::vector<double> E_tilde = evaluate_VSH_series_and_FD( E_tilde_r_L, E_tilde_1_L, E_tilde_2_L, i, j );
					std::vector<double> E_corrected = evaluate_VSH_series_and_FD( E_corrected_r_L, E_corrected_1_L, E_corrected_2_L, i, j );
					double sqrt_Bsq_over_E_tilde_sq = sqrt( ( B_x[i][j]*B_x[i][j] + B_y[i][j]*B_y[i][j] + B_z[i][j]*B_z[i][j] ) / ( E_tilde[3]*E_tilde[3] + E_tilde[4]*E_tilde[4] + E_tilde[5]*E_tilde[5] ) );
					std::vector<double> Ecorcheck { sqrt_Bsq_over_E_tilde_sq*E_tilde[3], sqrt_Bsq_over_E_tilde_sq*E_tilde[4], sqrt_Bsq_over_E_tilde_sq*E_tilde[5] };
					
					output_file <<","<< E_corrected[3] <<","<< E_corrected[4] <<","<< E_corrected[5];
					output_file <<","<< Ecorcheck[0] <<","<< Ecorcheck[1] <<","<< Ecorcheck[2];
					output_file <<","<< E_dot_B[i][j] <<","<< Bsq_minus_Esq[i][j] <<","<< sgn( Bsq_minus_Esq[i][j] );
					
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
			  <<std::left<<std::setw(w)<< "E dot B" <<std::left<<std::setw(w)<< "B^2 minus E^2" <<"|\t"
			  <<std::left<<std::setw(w)<< "% done" <<std::left<<std::setw(w)<< "Est. time left"
			  << std::endl;
}




void output_to_screen( int T_index ){
	
	if( ( T_index % cout_freq_T == 0 ) or ( T_index == n_timesteps ) ){
	
		double percent_done = round( (double)1000*T_index/n_timesteps )/10; // This format makes it exactly one decimal place.
		
		std::cout <<std::left<<std::setw(w)<< T_index <<std::left<<std::setw(w)<< T <<"|\t"
				  <<std::left<<std::setw(w)<< B_r[cout_i][cout_j] <<std::left<<std::setw(w)<< B_t[cout_i][cout_j] <<std::left<<std::setw(w)<< B_p[cout_i][cout_j] <<"|\t"
				  <<std::left<<std::setw(w)<< E_r[cout_i][cout_j] <<std::left<<std::setw(w)<< E_t[cout_i][cout_j] <<std::left<<std::setw(w)<< E_p[cout_i][cout_j] <<"|\t"
				  <<std::left<<std::setw(w)<< E_dot_B[cout_i][cout_j] <<std::left<<std::setw(w)<< Bsq_minus_Esq[cout_i][cout_j] <<"|\t"
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