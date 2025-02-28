// cd OneDrive\PhD\Codes\20221108 Vector spherical harmonics\VSH series identities
// g++ Cross_product_of_VSH_series_01.cpp -o Cross_product_of_VSH_series_01

/*
Fit vector spherical harmonic series to two vectors A(theta,phi) and B(theta,phi) and calculate their cross product.
Compare this to an analytical expression derived on paper (Book 10, p31).
*/

#include <iostream>
#include <iomanip>
#include <chrono>
#include "../../Generate_Spherical_Harmonics.h"
#include "../../Vector_Operations.h"
#include "../Headers_vector_function/Vector_function_05.h"
#include "../Headers_vector_function/Vector_function_06.h"
#include "../Headers_VSH_series/VSH_series_calculation_101.h"	// Works from _101 inclusive.
#include "../Headers_VSH_series/VSH_series_evaluation_101.h"	// Works from _101 inclusive.

int main(){
	
	auto time_start = std::chrono::high_resolution_clock::now();
	
	//----- Define variables -----
	int ell_max = 3;
	int N_t     = 30;
	int N_p     = 30;
	int i_test  = 18;
	int j_test  = 20;
	
	std::vector<double> t_list ( N_t );
	std::vector<double> p_list ( N_p );
	
	std::vector< std::vector< std::vector< std::complex<double> > > > A( 3, std::vector< std::vector< std::complex<double> > > ( N_t, std::vector< std::complex<double> > ( N_p ) ) );		// Precalculated values of the function that we wish to fit.
	std::vector< std::vector< std::vector< std::complex<double> > > > B( 3, std::vector< std::vector< std::complex<double> > > ( N_t, std::vector< std::complex<double> > ( N_p ) ) );		// Precalculated values of the function that we wish to fit.
	
	// Totals of absolute values of fitted coefficients in each direction. No physical interpretation but useful as a quick verification tool.
	std::vector< std::complex<double> > coeffs_total_A ( 3 );
	std::vector< std::complex<double> > coeffs_total_B ( 3 );
	
	
	//----- Fail if the test indices are greater than the maximum number of points -----
	if( i_test > N_t-1 ){
		std::cout << "Test theta index (" << i_test << ") exceeds number of datapoints (" << N_t << "). Stopping the code to prevent undefined behaviour." << std::endl;
		return 1;
	}
	if( j_test > N_p-1 ){
		std::cout << "Test phi index (" << j_test << ") exceeds number of datapoints (" << N_p << "). Stopping the code to prevent undefined behaviour." << std::endl;
		return 1;
	}
	
	
	//----- Build coordinate lists -----
	for( int i=0; i<N_t; i++ ){
		t_list[i] = i * pi / ( (double) N_t - 1 );
	}
	for( int j=0; j<N_p; j++ ){
		p_list[j] = j * 2*pi / ( (double) N_p - 1 );
	}
	
	
	//----- Generate spherical harmonic coefficients -----
	generate_coeffs_spherical_harmonic( ell_max );
	
	
	//----- Precalculate vector A to be fitted -----
	for( int i=0; i<N_t; i++ ){
		for( int j=0; j<N_p; j++ ){
			
			A[0][i][j] = A_r_function( t_list[i], p_list[j] );
			A[1][i][j] = A_t_function( t_list[i], p_list[j] );
			A[2][i][j] = A_p_function( t_list[i], p_list[j] );
			
			B[0][i][j] = B_r_function( t_list[i], p_list[j] );
			B[1][i][j] = B_t_function( t_list[i], p_list[j] );
			B[2][i][j] = B_p_function( t_list[i], p_list[j] );
		}
	}
	
	
	//----- Calculate VSH series -----
	std::vector< std::vector< std::vector< std::complex<double> > > > coeffs_A = generate_vsh_series( A, t_list, p_list, ell_max );
	std::vector< std::vector< std::vector< std::complex<double> > > > coeffs_B = generate_vsh_series( B, t_list, p_list, ell_max );
	
	
	//----- Output VSH series -----
	for( int ell=0; ell<=ell_max; ell++ ){
		for( int m=-ell; m<=ell; m++ ){
			for( int c=0; c<3; c++ ){
				coeffs_total_A[c] += std::complex<double> { abs( coeffs_A[c][ell][ell+m].real() ), abs( coeffs_A[c][ell][ell+m].imag() ) };
				coeffs_total_B[c] += std::complex<double> { abs( coeffs_B[c][ell][ell+m].real() ), abs( coeffs_B[c][ell][ell+m].imag() ) };
			}
		}
	}
	
	std::cout << "Coeffs total A:\t" << coeffs_total_A[0] <<"\t"<< coeffs_total_A[1] <<"\t"<< coeffs_total_A[2] <<"\t\t"<< coeffs_total_A[0]+coeffs_total_A[1]+coeffs_total_A[2] << std::endl;
	std::cout << "Coeffs total B:\t" << coeffs_total_B[0] <<"\t"<< coeffs_total_B[1] <<"\t"<< coeffs_total_B[2] <<"\t\t"<< coeffs_total_B[0]+coeffs_total_B[1]+coeffs_total_B[2] << std::endl;
	
	
	//----- Evaluate cross product at the test point -----
	std::cout << "\nTest fit for (theta,phi) = (" << t_list[i_test] <<","<< p_list[j_test] << "):" << std::endl;
	
	std::vector< std::complex<double> > A_test_calc  = evaluate_vsh_series( t_list[i_test], p_list[j_test], coeffs_A );
	std::vector< std::complex<double> > B_test_calc  = evaluate_vsh_series( t_list[i_test], p_list[j_test], coeffs_B );
	
	std::vector< std::complex<double> > cross_product_exact;
	std::vector< std::complex<double> > cross_product_from_VSH_series;
	
	cross_product_exact.push_back( A[1][i_test][j_test]*B[2][i_test][j_test] - A[2][i_test][j_test]*B[1][i_test][j_test] );
	cross_product_exact.push_back( A[2][i_test][j_test]*B[0][i_test][j_test] - A[0][i_test][j_test]*B[2][i_test][j_test] );
	cross_product_exact.push_back( A[0][i_test][j_test]*B[1][i_test][j_test] - A[1][i_test][j_test]*B[0][i_test][j_test] );
	
	cross_product_from_VSH_series.push_back( A_test_calc[1]*B_test_calc[2] - A_test_calc[2]*B_test_calc[1] );
	cross_product_from_VSH_series.push_back( A_test_calc[2]*B_test_calc[0] - A_test_calc[0]*B_test_calc[2] );
	cross_product_from_VSH_series.push_back( A_test_calc[0]*B_test_calc[1] - A_test_calc[1]*B_test_calc[0] );
	
	print_vector( cross_product_exact          , "Cross product exact          " );
	print_vector( cross_product_from_VSH_series, "Cross product from VSH series" );
	
	
	//----- Calculate dot product from VSH coefficients -----
	std::vector< std::complex<double> > cross_product_from_VSH_coeffs { 0, 0, 0 };
	
	double t = t_list[i_test];
	double p = p_list[j_test];
	
	for( int ell_1=0; ell_1<=ell_max; ell_1++ ){
		for( int m_1=-ell_1; m_1<=ell_1; m_1++ ){
			for( int ell_2=0; ell_2<=ell_max; ell_2++ ){
				for( int m_2=-ell_2; m_2<=ell_2; m_2++ ){
					
					
					//----- V1 -----
					/*std::vector<std::complex<double> > A_at_this_point (3);
					A_at_this_point = vector_sum( A_at_this_point, scalar_times_vector( coeffs_A[0][ell_1][ell_1+m_1], VSH_Y  (t,p,ell_1,m_1) ) );
					A_at_this_point = vector_sum( A_at_this_point, scalar_times_vector( coeffs_A[1][ell_1][ell_1+m_1], VSH_Psi(t,p,ell_1,m_1) ) );
					A_at_this_point = vector_sum( A_at_this_point, scalar_times_vector( coeffs_A[2][ell_1][ell_1+m_1], VSH_Phi(t,p,ell_1,m_1) ) );
					std::vector<std::complex<double> > B_at_this_point (3);
					B_at_this_point = vector_sum( B_at_this_point, scalar_times_vector( coeffs_B[0][ell_2][ell_2+m_2], VSH_Y  (t,p,ell_2,m_2) ) );
					B_at_this_point = vector_sum( B_at_this_point, scalar_times_vector( coeffs_B[1][ell_2][ell_2+m_2], VSH_Psi(t,p,ell_2,m_2) ) );
					B_at_this_point = vector_sum( B_at_this_point, scalar_times_vector( coeffs_B[2][ell_2][ell_2+m_2], VSH_Phi(t,p,ell_2,m_2) ) );
					cross_product_from_VSH_coeffs = vector_sum( cross_product_from_VSH_coeffs, cross_product( A_at_this_point, B_at_this_point ) );*/
					
					
					//----- V2 -----
					/*cross_product_from_VSH_coeffs = vector_sum( cross_product_from_VSH_coeffs, scalar_times_vector( coeffs_A[0][ell_1][ell_1+m_1] * coeffs_B[0][ell_2][ell_2+m_2], cross_product( VSH_Y  (t,p,ell_1,m_1), VSH_Y  (t,p,ell_2,m_2) ) ) );
					cross_product_from_VSH_coeffs = vector_sum( cross_product_from_VSH_coeffs, scalar_times_vector( coeffs_A[0][ell_1][ell_1+m_1] * coeffs_B[1][ell_2][ell_2+m_2], cross_product( VSH_Y  (t,p,ell_1,m_1), VSH_Psi(t,p,ell_2,m_2) ) ) );
					cross_product_from_VSH_coeffs = vector_sum( cross_product_from_VSH_coeffs, scalar_times_vector( coeffs_A[0][ell_1][ell_1+m_1] * coeffs_B[2][ell_2][ell_2+m_2], cross_product( VSH_Y  (t,p,ell_1,m_1), VSH_Phi(t,p,ell_2,m_2) ) ) );
					
					cross_product_from_VSH_coeffs = vector_sum( cross_product_from_VSH_coeffs, scalar_times_vector( coeffs_A[1][ell_1][ell_1+m_1] * coeffs_B[0][ell_2][ell_2+m_2], cross_product( VSH_Psi(t,p,ell_1,m_1), VSH_Y  (t,p,ell_2,m_2) ) ) );
					cross_product_from_VSH_coeffs = vector_sum( cross_product_from_VSH_coeffs, scalar_times_vector( coeffs_A[1][ell_1][ell_1+m_1] * coeffs_B[1][ell_2][ell_2+m_2], cross_product( VSH_Psi(t,p,ell_1,m_1), VSH_Psi(t,p,ell_2,m_2) ) ) );
					cross_product_from_VSH_coeffs = vector_sum( cross_product_from_VSH_coeffs, scalar_times_vector( coeffs_A[1][ell_1][ell_1+m_1] * coeffs_B[2][ell_2][ell_2+m_2], cross_product( VSH_Psi(t,p,ell_1,m_1), VSH_Phi(t,p,ell_2,m_2) ) ) );
					
					cross_product_from_VSH_coeffs = vector_sum( cross_product_from_VSH_coeffs, scalar_times_vector( coeffs_A[2][ell_1][ell_1+m_1] * coeffs_B[0][ell_2][ell_2+m_2], cross_product( VSH_Phi(t,p,ell_1,m_1), VSH_Y  (t,p,ell_2,m_2) ) ) );
					cross_product_from_VSH_coeffs = vector_sum( cross_product_from_VSH_coeffs, scalar_times_vector( coeffs_A[2][ell_1][ell_1+m_1] * coeffs_B[1][ell_2][ell_2+m_2], cross_product( VSH_Phi(t,p,ell_1,m_1), VSH_Psi(t,p,ell_2,m_2) ) ) );
					cross_product_from_VSH_coeffs = vector_sum( cross_product_from_VSH_coeffs, scalar_times_vector( coeffs_A[2][ell_1][ell_1+m_1] * coeffs_B[2][ell_2][ell_2+m_2], cross_product( VSH_Phi(t,p,ell_1,m_1), VSH_Phi(t,p,ell_2,m_2) ) ) );*/
					
					
					//----- V3: (1) and (2) directions will be kept the same; develop r direction within subversions -----
					//--- Subversions with developing r direction ---
					//- V3a -
					/*cross_product_from_VSH_coeffs[0] += coeffs_A[1][ell_1][ell_1+m_1] * coeffs_B[1][ell_2][ell_2+m_2] *   dot_product( VSH_Phi(t,p,ell_1,m_1), VSH_Psi(t,p,ell_2,m_2) );
					cross_product_from_VSH_coeffs[0] += coeffs_A[1][ell_1][ell_1+m_1] * coeffs_B[2][ell_2][ell_2+m_2] *   dot_product( VSH_Psi(t,p,ell_1,m_1), VSH_Psi(t,p,ell_2,m_2) );
					cross_product_from_VSH_coeffs[0] += coeffs_A[2][ell_1][ell_1+m_1] * coeffs_B[1][ell_2][ell_2+m_2] * - dot_product( VSH_Psi(t,p,ell_1,m_1), VSH_Psi(t,p,ell_2,m_2) );
					cross_product_from_VSH_coeffs[0] += coeffs_A[2][ell_1][ell_1+m_1] * coeffs_B[2][ell_2][ell_2+m_2] *   dot_product( VSH_Phi(t,p,ell_1,m_1), VSH_Psi(t,p,ell_2,m_2) );*/
					
					//- V3b -
					/*cross_product_from_VSH_coeffs[0] += ( coeffs_A[1][ell_1][ell_1+m_1] * coeffs_B[1][ell_2][ell_2+m_2] + coeffs_A[2][ell_1][ell_1+m_1] * coeffs_B[2][ell_2][ell_2+m_2] ) * dot_product( VSH_Phi(t,p,ell_1,m_1), VSH_Psi(t,p,ell_2,m_2) );
					cross_product_from_VSH_coeffs[0] += ( coeffs_A[1][ell_1][ell_1+m_1] * coeffs_B[2][ell_2][ell_2+m_2] - coeffs_A[2][ell_1][ell_1+m_1] * coeffs_B[1][ell_2][ell_2+m_2] ) * dot_product( VSH_Psi(t,p,ell_1,m_1), VSH_Psi(t,p,ell_2,m_2) );*/
					
					//- V3c -
					/*std::vector< std::complex<double> > Phi_11 = { 0, 0, ylm_dtheta(t,p,ell_1,m_1) };
					std::vector< std::complex<double> > Psi_11 = { 0, ylm_dtheta(t,p,ell_1,m_1), 0 };
					if( m_1 != 0 ){
						Phi_11[1] = std::complex<double>{0,(double)-m_1} * ylm_over_sintheta(t,p,ell_1,m_1);
						Psi_11[2] = std::complex<double>{0,(double) m_1} * ylm_over_sintheta(t,p,ell_1,m_1);
					}
					std::complex<double> dot_prod_Phi = dot_product( Phi_11, VSH_Psi(t,p,ell_2,m_2) );
					dot_prod_Phi *= coeffs_A[1][ell_1][ell_1+m_1]*coeffs_B[1][ell_2][ell_2+m_2] + coeffs_A[2][ell_1][ell_1+m_1]*coeffs_B[2][ell_2][ell_2+m_2];
					std::complex<double> dot_prod_Psi = dot_product( Psi_11, VSH_Psi(t,p,ell_2,m_2) );
					dot_prod_Psi *= coeffs_A[1][ell_1][ell_1+m_1]*coeffs_B[2][ell_2][ell_2+m_2] - coeffs_A[2][ell_1][ell_1+m_1]*coeffs_B[1][ell_2][ell_2+m_2];
					cross_product_from_VSH_coeffs[0] += dot_prod_Phi + dot_prod_Psi;*/
					
					//- V3d -
					/*std::vector< std::complex<double> > Phi_11 = { 0, 0, ylm_dtheta(t,p,ell_1,m_1) };
					std::vector< std::complex<double> > Psi_11 = { 0, ylm_dtheta(t,p,ell_1,m_1), 0 };
					if( m_1 != 0 ){
						Phi_11[1] = std::complex<double>{0,(double)-m_1} * ylm_over_sintheta(t,p,ell_1,m_1);
						Psi_11[2] = std::complex<double>{0,(double) m_1} * ylm_over_sintheta(t,p,ell_1,m_1);
					}
					Phi_11 = scalar_times_vector( coeffs_A[1][ell_1][ell_1+m_1]*coeffs_B[1][ell_2][ell_2+m_2] + coeffs_A[2][ell_1][ell_1+m_1]*coeffs_B[2][ell_2][ell_2+m_2], Phi_11 );
					Psi_11 = scalar_times_vector( coeffs_A[1][ell_1][ell_1+m_1]*coeffs_B[2][ell_2][ell_2+m_2] - coeffs_A[2][ell_1][ell_1+m_1]*coeffs_B[1][ell_2][ell_2+m_2], Psi_11 );
					std::complex<double> dot_prod_Phi = dot_product( Phi_11, VSH_Psi(t,p,ell_2,m_2) );
					std::complex<double> dot_prod_Psi = dot_product( Psi_11, VSH_Psi(t,p,ell_2,m_2) );
					cross_product_from_VSH_coeffs[0] += dot_prod_Phi + dot_prod_Psi;*/
					
					//- V3e -
					std::vector< std::complex<double> > Phi_11 = { 0, 0, ylm_dtheta(t,p,ell_1,m_1) };
					std::vector< std::complex<double> > Psi_11 = { 0, ylm_dtheta(t,p,ell_1,m_1), 0 };
					if( m_1 != 0 ){
						Phi_11[1] = std::complex<double>{0,(double)-m_1} * ylm_over_sintheta(t,p,ell_1,m_1);
						Psi_11[2] = std::complex<double>{0,(double) m_1} * ylm_over_sintheta(t,p,ell_1,m_1);
					}
					Phi_11 = scalar_times_vector( coeffs_A[1][ell_1][ell_1+m_1]*coeffs_B[1][ell_2][ell_2+m_2] + coeffs_A[2][ell_1][ell_1+m_1]*coeffs_B[2][ell_2][ell_2+m_2], Phi_11 );
					Psi_11 = scalar_times_vector( coeffs_A[1][ell_1][ell_1+m_1]*coeffs_B[2][ell_2][ell_2+m_2] - coeffs_A[2][ell_1][ell_1+m_1]*coeffs_B[1][ell_2][ell_2+m_2], Psi_11 );
					std::vector< std::complex<double> > Phi_plus_Psi = vector_sum( Phi_11, Psi_11 );
					cross_product_from_VSH_coeffs[0] += dot_product( Phi_plus_Psi, VSH_Psi(t,p,ell_2,m_2) );
					
					
					//- V3f (doesn't work 20221210 ) -
					/*std::vector< std::complex<double> > Phi_11 = { 0, 0, ylm_dtheta(t,p,ell_1,m_1) };
					std::vector< std::complex<double> > Psi_11 = { 0, ylm_dtheta(t,p,ell_1,m_1), 0 };
					if( m_1 != 0 ){
						Phi_11[1] = std::complex<double>{0,(double)-m_1} * ylm_over_sintheta(t,p,ell_1,m_1);
						Psi_11[2] = std::complex<double>{0,(double) m_1} * ylm_over_sintheta(t,p,ell_1,m_1);
					}
					Phi_11 = scalar_times_vector( coeffs_A[1][ell_1][ell_1+m_1]*coeffs_B[1][ell_2][ell_2+m_2] + coeffs_A[2][ell_1][ell_1+m_1]*coeffs_B[2][ell_2][ell_2+m_2], Phi_11 );
					Psi_11 = scalar_times_vector( coeffs_A[1][ell_1][ell_1+m_1]*coeffs_B[2][ell_2][ell_2+m_2] - coeffs_A[2][ell_1][ell_1+m_1]*coeffs_B[1][ell_2][ell_2+m_2], Phi_11 );
					std::vector< std::complex<double> > r_part = vector_sum( Phi_11, Psi_11 );
					cross_product_from_VSH_coeffs[0] += dot_product( r_part, VSH_Psi(t,p,ell_2,m_2) );*/
					
					// V3g (doesn't work 20221210 ) -
					/*std::vector< std::complex<double> > V_Phi = { 0, 0, 0 };
					std::vector< std::complex<double> > V_Psi = { 0, 0, 0 };
					V_Phi[1] = ( coeffs_A[1][ell_1][ell_1+m_1]*coeffs_B[1][ell_2][ell_2+m_2] + coeffs_A[2][ell_1][ell_1+m_1]*coeffs_B[2][ell_2][ell_2+m_2] ) * ylm_dtheta(t,p,ell_1,m_1);
					if( m_1 != 0 ){
						V_Phi[2] = ( coeffs_A[1][ell_1][ell_1+m_1]*coeffs_B[1][ell_2][ell_2+m_2] + coeffs_A[2][ell_1][ell_1+m_1]*coeffs_B[2][ell_2][ell_2+m_2] ) * std::complex<double>{0,(double) m_1}*ylm_over_sintheta(t,p,ell_1,m_1);
					}
					V_Psi[2] = ( coeffs_A[1][ell_1][ell_1+m_1]*coeffs_B[2][ell_2][ell_2+m_2] - coeffs_A[2][ell_1][ell_1+m_1]*coeffs_B[1][ell_2][ell_2+m_2] ) * ylm_dtheta(t,p,ell_1,m_1);
					if( m_2 != 0 ){
						V_Psi[1] = ( coeffs_A[1][ell_1][ell_1+m_1]*coeffs_B[2][ell_2][ell_2+m_2] - coeffs_A[2][ell_1][ell_1+m_1]*coeffs_B[1][ell_2][ell_2+m_2] ) * std::complex<double>{0,(double)-m_1}*ylm_over_sintheta(t,p,ell_1,m_1);
					}
					cross_product_from_VSH_coeffs[0] += dot_product( vector_sum( V_Phi, V_Psi ), VSH_Psi(t,p,ell_2,m_2) );*/
					
					//- V3h (doesn't work 20221210 ) -
					/*std::vector< std::complex<double> > r_part = { 0, 0, 0 };
					r_part[1] = ( coeffs_A[1][ell_1][ell_1+m_1]*coeffs_B[1][ell_2][ell_2+m_2] + coeffs_A[2][ell_1][ell_1+m_1]*coeffs_B[2][ell_2][ell_2+m_2] ) * ylm_dtheta(t,p,ell_1,ell_2);
					r_part[2] = ( coeffs_A[1][ell_1][ell_1+m_1]*coeffs_B[2][ell_2][ell_2+m_2] - coeffs_A[2][ell_1][ell_1+m_1]*coeffs_B[1][ell_2][ell_2+m_2] ) * ylm_dtheta(t,p,ell_1,ell_2);
					if( m_1 != 0 ){
						r_part[1] += ( coeffs_A[1][ell_1][ell_1+m_1]*coeffs_B[2][ell_2][ell_2+m_2] - coeffs_A[2][ell_1][ell_1+m_1]*coeffs_B[1][ell_2][ell_2+m_2] ) * std::complex<double>{0,(double)-m_1}*ylm_over_sintheta(t,p,ell_1,m_1);
						r_part[2] += ( coeffs_A[1][ell_1][ell_1+m_1]*coeffs_B[1][ell_2][ell_2+m_2] + coeffs_A[2][ell_1][ell_1+m_1]*coeffs_B[2][ell_2][ell_2+m_2] ) * std::complex<double>{0,(double) m_1}*ylm_over_sintheta(t,p,ell_1,m_1);
					}
					print_vector( r_part );
					cross_product_from_VSH_coeffs[0] = dot_product( r_part, VSH_Psi(t,p,ell_2,m_2) );*/
					
					//- V3e -
					/*cross_product_from_VSH_coeffs[0] += ( coeffs_A[2][ell_1][ell_1+m_1] * coeffs_B[1][ell_2][ell_2+m_2]  - coeffs_A[1][ell_1][ell_1+m_1] * coeffs_B[2][ell_2][ell_2+m_2] ) * ylm_dtheta(t,p,ell_1,m_1) * ylm_dtheta(t,p,ell_2,m_2);
					if( m_1 != 0 ){
						cross_product_from_VSH_coeffs[0] += ( coeffs_A[1][ell_1][ell_1+m_1] * coeffs_B[1][ell_2][ell_2+m_2] + coeffs_A[2][ell_1][ell_1+m_1] * coeffs_B[2][ell_2][ell_2+m_2] ) * std::complex<double>{0,(double) m_2} * ylm_over_sintheta(t,p,ell_2,m_2) * ylm_dtheta(t,p,ell_1,m_1);
					}
					if( m_2 != 0 ){
						cross_product_from_VSH_coeffs[0] += ( coeffs_A[1][ell_1][ell_1+m_1] * coeffs_B[1][ell_2][ell_2+m_2] + coeffs_A[2][ell_1][ell_1+m_1] * coeffs_B[2][ell_2][ell_2+m_2] ) * std::complex<double>{0,(double)-m_1} * ylm_over_sintheta(t,p,ell_1,m_1) * ylm_dtheta(t,p,ell_2,m_2);
					}
					if( ( m_1 != 0 ) && ( m_2 != 0 ) ){
						cross_product_from_VSH_coeffs[0] += ( coeffs_A[1][ell_1][ell_1+m_1] * coeffs_B[2][ell_2][ell_2+m_2] - coeffs_A[2][ell_1][ell_1+m_1] * coeffs_B[1][ell_2][ell_2+m_2] ) * ((double)m_1*m_2) * ylm_over_sintheta(t,p,ell_1,m_1) * ylm_over_sintheta(t,p,ell_2,m_2);
					}*/
					
					
					//--- V3 (1) and (2) directions common to all subversions ---
					cross_product_from_VSH_coeffs = vector_sum( cross_product_from_VSH_coeffs, scalar_times_vector(  coeffs_A[2][ell_1][ell_1+m_1]*coeffs_B[0][ell_2][ell_2+m_2]*ylm(t,p,ell_2,m_2), VSH_Psi(t,p,ell_1,m_1) ) );
					cross_product_from_VSH_coeffs = vector_sum( cross_product_from_VSH_coeffs, scalar_times_vector( -coeffs_A[0][ell_1][ell_1+m_1]*coeffs_B[2][ell_2][ell_2+m_2]*ylm(t,p,ell_1,m_1), VSH_Psi(t,p,ell_2,m_2) ) );
					
					cross_product_from_VSH_coeffs = vector_sum( cross_product_from_VSH_coeffs, scalar_times_vector( -coeffs_A[1][ell_1][ell_1+m_1]*coeffs_B[0][ell_2][ell_2+m_2]*ylm(t,p,ell_2,m_2), VSH_Phi(t,p,ell_1,m_1) ) );
					cross_product_from_VSH_coeffs = vector_sum( cross_product_from_VSH_coeffs, scalar_times_vector(  coeffs_A[0][ell_1][ell_1+m_1]*coeffs_B[1][ell_2][ell_2+m_2]*ylm(t,p,ell_1,m_1), VSH_Phi(t,p,ell_2,m_2) ) );
				}
			}
		}
	}
	print_vector( cross_product_from_VSH_coeffs, "Cross product from VSH coeffs" );
	
	
	
	//----- Output execution time and finish code -----
	auto time_stop = std::chrono::high_resolution_clock::now();
	double exec_time = std::chrono::duration_cast<std::chrono::nanoseconds>( time_stop - time_start ).count() * 1e-9;
	std::cout << "\nExecution time:\t" << floor(exec_time/60) << " m " << exec_time-60*floor(exec_time/60) << " s." << std::endl;
	
	std::cout << "\nFinished" << std::endl;
}