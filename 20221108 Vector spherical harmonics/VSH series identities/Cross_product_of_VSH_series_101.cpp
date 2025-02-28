// cd OneDrive\PhD\Codes\20221108 Vector spherical harmonics\VSH series identities
// g++ Cross_product_of_VSH_series_101.cpp -o Cross_product_of_VSH_series_101

/*
Fit vector spherical harmonic series to two vectors A(theta,phi) and B(theta,phi) and calculate their cross product.
Compare this to an analytical expression derived on paper (Book 10, p56).
Different derivation to the previous code.
*/

#include <iostream>
#include <iomanip>
#include <chrono>
#include "../../Generate_Spherical_Harmonics.h"
//#include "../../Vector_Operations.h"
#include "../Headers_vector_function/Vector_function_05.h"
#include "../Headers_vector_function/Vector_function_06.h"
#include "../Headers_VSH_series/VSH_series_calculation_101.h"
#include "../Headers_VSH_series/VSH_series_evaluation_101.h"

int main(){
	
	auto time_start = std::chrono::high_resolution_clock::now();
	
	//----- Define variables -----
	int ell_max = 8;
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
	double t = t_list[i_test];
	double p = p_list[j_test];
	
	std::cout << "\nTest fit for (theta,phi) = (" << t <<","<< p << "):" << std::endl;
	
	std::vector< std::complex<double> > A_test_calc  = evaluate_vsh_series( t, p, coeffs_A );
	std::vector< std::complex<double> > B_test_calc  = evaluate_vsh_series( t, p, coeffs_B );
	
	std::vector< std::complex<double> > A_test_exact{ A[0][i_test][j_test], A[1][i_test][j_test], A[2][i_test][j_test] };
	std::vector< std::complex<double> > B_test_exact{ B[0][i_test][j_test], B[1][i_test][j_test], B[2][i_test][j_test] };
	
	std::vector< std::complex<double> > cross_product_exact{ A_test_exact[1]*B_test_exact[2]-A_test_exact[2]*B_test_exact[1],
	                                                         A_test_exact[2]*B_test_exact[0]-A_test_exact[0]*B_test_exact[2],
															 A_test_exact[0]*B_test_exact[1]-A_test_exact[1]*B_test_exact[0] };
	
	std::vector< std::complex<double> > cross_product_from_VSH_series{ A_test_calc[1]*B_test_calc[2]-A_test_calc[2]*B_test_calc[1],
	                                                                   A_test_calc[2]*B_test_calc[0]-A_test_calc[0]*B_test_calc[2],
															           A_test_calc[0]*B_test_calc[1]-A_test_calc[1]*B_test_calc[0] };
	
	print_vector( cross_product_exact          , "Cross product exact          " );
	print_vector( cross_product_from_VSH_series, "Cross product from VSH series" );
	
	
	//----- Calculate cross product from VSH coefficients -----
	std::vector< std::complex<double> > cross_product_from_VSH_coeffs { 0, 0, 0 };
	
	std::vector< std::complex<double> > cross_product_from_VSH_coeffs_00 { 0, 0, 0 };
	std::vector< std::complex<double> > cross_product_from_VSH_coeffs_01 { 0, 0, 0 };
	std::vector< std::complex<double> > cross_product_from_VSH_coeffs_02 { 0, 0, 0 };
	std::vector< std::complex<double> > cross_product_from_VSH_coeffs_10 { 0, 0, 0 };
	std::vector< std::complex<double> > cross_product_from_VSH_coeffs_11 { 0, 0, 0 };
	std::vector< std::complex<double> > cross_product_from_VSH_coeffs_12 { 0, 0, 0 };
	std::vector< std::complex<double> > cross_product_from_VSH_coeffs_20 { 0, 0, 0 };
	std::vector< std::complex<double> > cross_product_from_VSH_coeffs_21 { 0, 0, 0 };
	std::vector< std::complex<double> > cross_product_from_VSH_coeffs_22 { 0, 0, 0 };
	
	for( int L1=0; L1<=ell_max; L1++ ){
		for( int M1=-L1; M1<=L1; M1++ ){
			for( int L2=0; L2<=ell_max; L2++ ){
				for( int M2=-L2; M2<=L2; M2++ ){
					
					//--- Precalculate vectors and VSHs at the test point for these indices ---
					std::vector< std::complex<double> > Ar12 { coeffs_A[0][L1][L1+M1], coeffs_A[1][L1][L1+M1], coeffs_A[2][L1][L1+M1] };
					std::vector< std::complex<double> > Br12 { coeffs_B[0][L2][L2+M2], coeffs_B[1][L2][L2+M2], coeffs_B[2][L2][L2+M2] };
					std::vector< std::complex<double> > Y1   = VSH_Y  (t,p,L1,M1);
					std::vector< std::complex<double> > Y2   = VSH_Y  (t,p,L2,M2);
					std::vector< std::complex<double> > Psi1 = VSH_Psi(t,p,L1,M1);
					std::vector< std::complex<double> > Psi2 = VSH_Psi(t,p,L2,M2);
					std::vector< std::complex<double> > Phi1 = VSH_Phi(t,p,L1,M1);
					std::vector< std::complex<double> > Phi2 = VSH_Phi(t,p,L2,M2);
					
					//--- Add to each vector in the sum ---
					cross_product_from_VSH_coeffs_00 = vector_sum( cross_product_from_VSH_coeffs_00, scalar_times_vector( Ar12[0] * Br12[0], cross_product( Y1  , Y2   ) ) );
					cross_product_from_VSH_coeffs_01 = vector_sum( cross_product_from_VSH_coeffs_01, scalar_times_vector( Ar12[0] * Br12[1], cross_product( Y1  , Psi2 ) ) );
					cross_product_from_VSH_coeffs_02 = vector_sum( cross_product_from_VSH_coeffs_02, scalar_times_vector( Ar12[0] * Br12[2], cross_product( Y1  , Phi2 ) ) );
					cross_product_from_VSH_coeffs_10 = vector_sum( cross_product_from_VSH_coeffs_10, scalar_times_vector( Ar12[1] * Br12[0], cross_product( Psi1, Y2   ) ) );
					cross_product_from_VSH_coeffs_11 = vector_sum( cross_product_from_VSH_coeffs_11, scalar_times_vector( Ar12[1] * Br12[1], cross_product( Psi1, Psi2 ) ) );
					cross_product_from_VSH_coeffs_12 = vector_sum( cross_product_from_VSH_coeffs_12, scalar_times_vector( Ar12[1] * Br12[2], cross_product( Psi1, Phi2 ) ) );
					cross_product_from_VSH_coeffs_20 = vector_sum( cross_product_from_VSH_coeffs_20, scalar_times_vector( Ar12[2] * Br12[0], cross_product( Phi1, Y2   ) ) );
					cross_product_from_VSH_coeffs_21 = vector_sum( cross_product_from_VSH_coeffs_21, scalar_times_vector( Ar12[2] * Br12[1], cross_product( Phi1, Psi2 ) ) );
					cross_product_from_VSH_coeffs_22 = vector_sum( cross_product_from_VSH_coeffs_22, scalar_times_vector( Ar12[2] * Br12[2], cross_product( Phi1, Phi2 ) ) );
				}
			}
		}
	}
	
	cross_product_from_VSH_coeffs = vector_sum( cross_product_from_VSH_coeffs, cross_product_from_VSH_coeffs_00 );
	cross_product_from_VSH_coeffs = vector_sum( cross_product_from_VSH_coeffs, cross_product_from_VSH_coeffs_01 );
	cross_product_from_VSH_coeffs = vector_sum( cross_product_from_VSH_coeffs, cross_product_from_VSH_coeffs_02 );
	cross_product_from_VSH_coeffs = vector_sum( cross_product_from_VSH_coeffs, cross_product_from_VSH_coeffs_10 );
	cross_product_from_VSH_coeffs = vector_sum( cross_product_from_VSH_coeffs, cross_product_from_VSH_coeffs_11 );
	cross_product_from_VSH_coeffs = vector_sum( cross_product_from_VSH_coeffs, cross_product_from_VSH_coeffs_12 );
	cross_product_from_VSH_coeffs = vector_sum( cross_product_from_VSH_coeffs, cross_product_from_VSH_coeffs_20 );
	cross_product_from_VSH_coeffs = vector_sum( cross_product_from_VSH_coeffs, cross_product_from_VSH_coeffs_21 );
	cross_product_from_VSH_coeffs = vector_sum( cross_product_from_VSH_coeffs, cross_product_from_VSH_coeffs_22 );
	
	print_vector( cross_product_from_VSH_coeffs, "Cross product from VSH coeffs" );
	
	std::cout << "\n\n" << std::endl;
	
	int w = 20;
	print_vector_setw( cross_product_from_VSH_coeffs_00, w, "Y_1   x Y_2  " );
	print_vector_setw( cross_product_from_VSH_coeffs_01, w, "Y_1   x Psi_2" );
	print_vector_setw( cross_product_from_VSH_coeffs_02, w, "Y_1   x Phi_2" );
	print_vector_setw( cross_product_from_VSH_coeffs_10, w, "Psi_1 x Y_2  " );
	print_vector_setw( cross_product_from_VSH_coeffs_11, w, "Psi_1 x Psi_2" );
	print_vector_setw( cross_product_from_VSH_coeffs_12, w, "Psi_1 x Phi_2" );
	print_vector_setw( cross_product_from_VSH_coeffs_20, w, "Phi_1 x Y_2  " );
	print_vector_setw( cross_product_from_VSH_coeffs_21, w, "Phi_1 x Psi_2" );
	print_vector_setw( cross_product_from_VSH_coeffs_22, w, "Phi_1 x Phi_2" );
	
	
	
	//----- Output execution time and finish code -----
	auto time_stop = std::chrono::high_resolution_clock::now();
	double exec_time = std::chrono::duration_cast<std::chrono::nanoseconds>( time_stop - time_start ).count() * 1e-9;
	std::cout << "\nExecution time:\t" << floor(exec_time/60) << " m " << exec_time-60*floor(exec_time/60) << " s." << std::endl;
	
	std::cout << "\nFinished" << std::endl;
}