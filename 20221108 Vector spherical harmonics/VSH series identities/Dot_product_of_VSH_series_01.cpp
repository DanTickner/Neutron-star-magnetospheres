// cd OneDrive\PhD\Codes\20221108 Vector spherical harmonics\VSH series identities
// g++ Dot_product_of_VSH_series_01.cpp -o Dot_product_of_VSH_series_01 -std=c++11

/*
Fit vector spherical harmonic series to two vectors A(theta,phi) and B(theta,phi) and calculate their dot product.
Compare this to an analytical expression derived on paper (Book 9, p59).
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
	
	int test_t_i = 18;
	int test_p_j = 20;
	
	std::vector<double> t_list ( N_t );
	std::vector<double> p_list ( N_p );
	
	std::vector< std::vector< std::vector< std::complex<double> > > > A( 3, std::vector< std::vector< std::complex<double> > > ( N_t, std::vector< std::complex<double> > ( N_p ) ) );		// Precalculated values of the function that we wish to fit.
	std::vector< std::vector< std::vector< std::complex<double> > > > B( 3, std::vector< std::vector< std::complex<double> > > ( N_t, std::vector< std::complex<double> > ( N_p ) ) );		// Precalculated values of the function that we wish to fit.
	
	// Totals of absolute values of fitted coefficients in each direction. No physical interpretation but useful as a quick verification tool.
	std::vector< std::complex<double> > coeffs_total_A ( 3 );
	std::vector< std::complex<double> > coeffs_total_B ( 3 );
	
	
	//----- Fail if the test indices are greater than the maximum number of points -----
	if( test_t_i > N_t-1 ){
		std::cout << "Test theta index (" << test_t_i << ") exceeds number of datapoints (" << N_t << "). Stopping the code to prevent undefined behaviour." << std::endl;
		return 1;
	}
	if( test_p_j > N_p-1 ){
		std::cout << "Test phi index (" << test_p_j << ") exceeds number of datapoints (" << N_p << "). Stopping the code to prevent undefined behaviour." << std::endl;
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
	std::vector< std::vector< std::vector< std::complex<double> > > > coeffs_VSH_series_A = generate_vsh_series( A, t_list, p_list, ell_max );
	std::vector< std::vector< std::vector< std::complex<double> > > > coeffs_VSH_series_B = generate_vsh_series( B, t_list, p_list, ell_max );
	
	
	//----- Output VSH series -----
	for( int ell=0; ell<=ell_max; ell++ ){
		for( int m=-ell; m<=ell; m++ ){
			for( int c=0; c<3; c++ ){
				coeffs_total_A[c] += std::complex<double> { abs( coeffs_VSH_series_A[c][ell][ell+m].real() ), abs( coeffs_VSH_series_A[c][ell][ell+m].imag() ) };
				coeffs_total_B[c] += std::complex<double> { abs( coeffs_VSH_series_B[c][ell][ell+m].real() ), abs( coeffs_VSH_series_B	[c][ell][ell+m].imag() ) };
			}
		}
	}
	
	std::cout << "Coeffs total A:\t" << coeffs_total_A[0] <<"\t"<< coeffs_total_A[1] <<"\t"<< coeffs_total_A[2] <<"\t\t"<< coeffs_total_A[0]+coeffs_total_A[1]+coeffs_total_A[2] << std::endl;
	std::cout << "Coeffs total B:\t" << coeffs_total_B[0] <<"\t"<< coeffs_total_B[1] <<"\t"<< coeffs_total_B[2] <<"\t\t"<< coeffs_total_B[0]+coeffs_total_B[1]+coeffs_total_B[2] << std::endl;
	
	
	//----- Evaluate dot product at the test point -----
	std::cout << "\nTest fit for (theta,phi) = (" << t_list[test_t_i] <<","<< p_list[test_p_j] << "):" << std::endl;
	
	std::vector< std::complex<double> > A_test_calc  = evaluate_vsh_series( t_list[test_t_i], p_list[test_p_j], coeffs_VSH_series_A );
	std::vector< std::complex<double> > B_test_calc  = evaluate_vsh_series( t_list[test_t_i], p_list[test_p_j], coeffs_VSH_series_B );
	
	std::complex<double> dot_product_exact;
	std::complex<double> dot_product_from_VSH_series;
	for( int c=0; c<3; c++ ){
		dot_product_exact += A[c][test_t_i][test_p_j] * B[c][test_t_i][test_p_j];
		dot_product_from_VSH_series += A_test_calc[c] * B_test_calc[c];
	}
	
	std::cout << "Dot product exact          :\t" << dot_product_exact           << std::endl;
	std::cout << "Dot product from VSH series:\t" << dot_product_from_VSH_series << std::endl;
	
	
	//----- Calculate dot product from VSH coefficients -----
	std::complex<double> dot_product_from_VSH_coeffs = 0;
	
	for( int ell_1=0; ell_1<=ell_max; ell_1++ ){
		for( int m_1=-ell_1; m_1<=ell_1; m_1++ ){
			for( int ell_2=0; ell_2<=ell_max; ell_2++ ){
				for( int m_2=-ell_2; m_2<=ell_2; m_2++ ){
					
					
					//----- V1 -----
					/*std::vector<std::complex<double> > A_at_this_point (3);
					A_at_this_point = vector_sum_cart( A_at_this_point, scalar_times_vector_cart( coeffs_VSH_series_A[0][ell_1][ell_1+m_1], VSH_Y  (t_list[test_t_i],p_list[test_p_j],ell_1,m_1) ) );
					A_at_this_point = vector_sum_cart( A_at_this_point, scalar_times_vector_cart( coeffs_VSH_series_A[1][ell_1][ell_1+m_1], VSH_Psi(t_list[test_t_i],p_list[test_p_j],ell_1,m_1) ) );
					A_at_this_point = vector_sum_cart( A_at_this_point, scalar_times_vector_cart( coeffs_VSH_series_A[2][ell_1][ell_1+m_1], VSH_Phi(t_list[test_t_i],p_list[test_p_j],ell_1,m_1) ) );
					std::vector<std::complex<double> > B_at_this_point (3);
					B_at_this_point = vector_sum_cart( B_at_this_point, scalar_times_vector_cart( coeffs_VSH_series_B[0][ell_2][ell_2+m_2], VSH_Y  (t_list[test_t_i],p_list[test_p_j],ell_2,m_2) ) );
					B_at_this_point = vector_sum_cart( B_at_this_point, scalar_times_vector_cart( coeffs_VSH_series_B[1][ell_2][ell_2+m_2], VSH_Psi(t_list[test_t_i],p_list[test_p_j],ell_2,m_2) ) );
					B_at_this_point = vector_sum_cart( B_at_this_point, scalar_times_vector_cart( coeffs_VSH_series_B[2][ell_2][ell_2+m_2], VSH_Phi(t_list[test_t_i],p_list[test_p_j],ell_2,m_2) ) );
					dot_product_from_VSH_coeffs += dot_product_cart( A_at_this_point, B_at_this_point );*/
					
					
					//----- V2 -----
					/*dot_product_from_VSH_coeffs += coeffs_VSH_series_A[0][ell_1][ell_1+m_1] * coeffs_VSH_series_B[0][ell_2][ell_2+m_2] * dot_product_cart( VSH_Y  (t_list[test_t_i],p_list[test_p_j],ell_1,m_1), VSH_Y  (t_list[test_t_i],p_list[test_p_j],ell_2,m_2) );
					dot_product_from_VSH_coeffs += coeffs_VSH_series_A[0][ell_1][ell_1+m_1] * coeffs_VSH_series_B[1][ell_2][ell_2+m_2] * dot_product_cart( VSH_Y  (t_list[test_t_i],p_list[test_p_j],ell_1,m_1), VSH_Psi(t_list[test_t_i],p_list[test_p_j],ell_2,m_2) );
					dot_product_from_VSH_coeffs += coeffs_VSH_series_A[0][ell_1][ell_1+m_1] * coeffs_VSH_series_B[2][ell_2][ell_2+m_2] * dot_product_cart( VSH_Y  (t_list[test_t_i],p_list[test_p_j],ell_1,m_1), VSH_Phi(t_list[test_t_i],p_list[test_p_j],ell_2,m_2) );
					
					dot_product_from_VSH_coeffs += coeffs_VSH_series_A[1][ell_1][ell_1+m_1] * coeffs_VSH_series_B[0][ell_2][ell_2+m_2] * dot_product_cart( VSH_Psi(t_list[test_t_i],p_list[test_p_j],ell_1,m_1), VSH_Y  (t_list[test_t_i],p_list[test_p_j],ell_2,m_2) );
					dot_product_from_VSH_coeffs += coeffs_VSH_series_A[1][ell_1][ell_1+m_1] * coeffs_VSH_series_B[1][ell_2][ell_2+m_2] * dot_product_cart( VSH_Psi(t_list[test_t_i],p_list[test_p_j],ell_1,m_1), VSH_Psi(t_list[test_t_i],p_list[test_p_j],ell_2,m_2) );
					dot_product_from_VSH_coeffs += coeffs_VSH_series_A[1][ell_1][ell_1+m_1] * coeffs_VSH_series_B[2][ell_2][ell_2+m_2] * dot_product_cart( VSH_Psi(t_list[test_t_i],p_list[test_p_j],ell_1,m_1), VSH_Phi(t_list[test_t_i],p_list[test_p_j],ell_2,m_2) );
					
					dot_product_from_VSH_coeffs += coeffs_VSH_series_A[2][ell_1][ell_1+m_1] * coeffs_VSH_series_B[0][ell_2][ell_2+m_2] * dot_product_cart( VSH_Phi(t_list[test_t_i],p_list[test_p_j],ell_1,m_1), VSH_Y  (t_list[test_t_i],p_list[test_p_j],ell_2,m_2) );
					dot_product_from_VSH_coeffs += coeffs_VSH_series_A[2][ell_1][ell_1+m_1] * coeffs_VSH_series_B[1][ell_2][ell_2+m_2] * dot_product_cart( VSH_Phi(t_list[test_t_i],p_list[test_p_j],ell_1,m_1), VSH_Psi(t_list[test_t_i],p_list[test_p_j],ell_2,m_2) );
					dot_product_from_VSH_coeffs += coeffs_VSH_series_A[2][ell_1][ell_1+m_1] * coeffs_VSH_series_B[2][ell_2][ell_2+m_2] * dot_product_cart( VSH_Phi(t_list[test_t_i],p_list[test_p_j],ell_1,m_1), VSH_Phi(t_list[test_t_i],p_list[test_p_j],ell_2,m_2) );*/
					
					
					//----- V3 -----
					/*dot_product_from_VSH_coeffs += coeffs_VSH_series_A[0][ell_1][ell_1+m_1] * coeffs_VSH_series_B[0][ell_2][ell_2+m_2] * dot_product_cart( VSH_Y  (t_list[test_t_i],p_list[test_p_j],ell_1,m_1), VSH_Y  (t_list[test_t_i],p_list[test_p_j],ell_2,m_2) );
					dot_product_from_VSH_coeffs += 0;
					dot_product_from_VSH_coeffs += 0;
					
					dot_product_from_VSH_coeffs += 0;
					dot_product_from_VSH_coeffs += coeffs_VSH_series_A[1][ell_1][ell_1+m_1] * coeffs_VSH_series_B[1][ell_2][ell_2+m_2] * dot_product_cart( VSH_Psi(t_list[test_t_i],p_list[test_p_j],ell_1,m_1), VSH_Psi(t_list[test_t_i],p_list[test_p_j],ell_2,m_2) );
					dot_product_from_VSH_coeffs += coeffs_VSH_series_A[1][ell_1][ell_1+m_1] * coeffs_VSH_series_B[2][ell_2][ell_2+m_2] * dot_product_cart( VSH_Psi(t_list[test_t_i],p_list[test_p_j],ell_1,m_1), VSH_Phi(t_list[test_t_i],p_list[test_p_j],ell_2,m_2) );
					
					dot_product_from_VSH_coeffs += 0;
					dot_product_from_VSH_coeffs += coeffs_VSH_series_A[2][ell_1][ell_1+m_1] * coeffs_VSH_series_B[1][ell_2][ell_2+m_2] * dot_product_cart( VSH_Phi(t_list[test_t_i],p_list[test_p_j],ell_1,m_1), VSH_Psi(t_list[test_t_i],p_list[test_p_j],ell_2,m_2) );
					dot_product_from_VSH_coeffs += coeffs_VSH_series_A[2][ell_1][ell_1+m_1] * coeffs_VSH_series_B[2][ell_2][ell_2+m_2] * dot_product_cart( VSH_Phi(t_list[test_t_i],p_list[test_p_j],ell_1,m_1), VSH_Phi(t_list[test_t_i],p_list[test_p_j],ell_2,m_2) );*/
					
					
					//----- V4 -----
					dot_product_from_VSH_coeffs += coeffs_VSH_series_A[0][ell_1][ell_1+m_1] * coeffs_VSH_series_B[0][ell_2][ell_2+m_2] * ylm(t_list[test_t_i],p_list[test_p_j],ell_1,m_1) * ylm(t_list[test_t_i],p_list[test_p_j],ell_2,m_2);
					
					dot_product_from_VSH_coeffs += ( coeffs_VSH_series_A[1][ell_1][ell_1+m_1] * coeffs_VSH_series_B[1][ell_2][ell_2+m_2] + coeffs_VSH_series_A[2][ell_1][ell_1+m_1] * coeffs_VSH_series_B[2][ell_2][ell_2+m_2] ) * (
					                               ylm_dtheta(t_list[test_t_i],p_list[test_p_j],ell_1,m_1) * ylm_dtheta(t_list[test_t_i],p_list[test_p_j],ell_2,m_2) + pow(sin(t_list[test_t_i]),-2) * ylm_dphi(t_list[test_t_i],p_list[test_p_j],ell_1,m_1) * ylm_dphi(t_list[test_t_i],p_list[test_p_j],ell_2,m_2) );
												   
					dot_product_from_VSH_coeffs += ( coeffs_VSH_series_A[2][ell_1][ell_1+m_1] * coeffs_VSH_series_B[1][ell_2][ell_2+m_2] - coeffs_VSH_series_A[1][ell_1][ell_1+m_1] * coeffs_VSH_series_B[2][ell_2][ell_2+m_2] ) * pow(sin(t_list[test_t_i]),-1) * (
					                               ylm_dtheta(t_list[test_t_i],p_list[test_p_j],ell_1,m_1) * ylm_dphi(t_list[test_t_i],p_list[test_p_j],ell_2,m_2) - ylm_dphi(t_list[test_t_i],p_list[test_p_j],ell_1,m_1) * ylm_dtheta(t_list[test_t_i],p_list[test_p_j],ell_2,m_2) );
				}
			}
		}
	}
	std::cout << "Dot product from VSH series:\t" << dot_product_from_VSH_coeffs << std::endl;
	
	
	//----- Output execution time and finish code -----
	auto time_stop = std::chrono::high_resolution_clock::now();
	double exec_time = std::chrono::duration_cast<std::chrono::nanoseconds>( time_stop - time_start ).count() * 1e-9;
	std::cout << "\nExecution time:\t" << floor(exec_time/60) << " m " << exec_time-60*floor(exec_time/60) << " s." << std::endl;
	
	std::cout << "\nFinished" << std::endl;
}