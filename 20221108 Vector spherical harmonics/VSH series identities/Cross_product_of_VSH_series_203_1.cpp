// cd OneDrive\PhD\Codes\20221108 Vector spherical harmonics\VSH series identities
// g++ Cross_product_of_VSH_series_203_1.cpp -o Cross_product_of_VSH_series_203_1

/*
Fit vector spherical harmonic series to two vectors A(theta,phi) and B(theta,phi) and calculate their cross product.
Compare this to an analytical expression derived on paper (Book 10, p58).
Calculate _01 sum separately.
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
	int k_max = pow(ell_max+1,2) - 1;
	std::vector< std::complex<double> > cross_product_from_VSH_coeffs { 0, 0, 0 };
	
	std::vector< std::complex<double> > cross_product_from_VSH_coeffs_01 { 0, 0, 0 };
	std::vector< std::complex<double> > cross_product_from_VSH_coeffs_02 { 0, 0, 0 };
	std::vector< std::complex<double> > cross_product_from_VSH_coeffs_10 { 0, 0, 0 };
	std::vector< std::complex<double> > cross_product_from_VSH_coeffs_11 { 0, 0, 0 };
	std::vector< std::complex<double> > cross_product_from_VSH_coeffs_12 { 0, 0, 0 };
	std::vector< std::complex<double> > cross_product_from_VSH_coeffs_20 { 0, 0, 0 };
	std::vector< std::complex<double> > cross_product_from_VSH_coeffs_21 { 0, 0, 0 };
	std::vector< std::complex<double> > cross_product_from_VSH_coeffs_22 { 0, 0, 0 };
	
	std::vector< std::complex<double> > cross_product_from_VSH_coeffs_01a { 0, 0, 0 };	// until calculation correct
	
	
	//----- Calculate _01 sum -----
	
	//--- V0 (same as if it were in the main loop) ---
	/*for( int k1=0; k1<=k_max; k1++ ){
		int L1 = sqrt(k1);
		int M1 = k1 - L1*(L1+1);
		for( int k2=0; k2<=k_max; k2++ ){	
			int L2 = sqrt(k2);
			int M2 = k2 - L2*(L2+1);
			//--- Precalculate vectors and VSHs at the test point for these indices ---
			std::vector< std::complex<double> > Ar12 { coeffs_A[0][L1][L1+M1], coeffs_A[1][L1][L1+M1], coeffs_A[2][L1][L1+M1] };
			std::vector< std::complex<double> > Br12 { coeffs_B[0][L2][L2+M2], coeffs_B[1][L2][L2+M2], coeffs_B[2][L2][L2+M2] };
			std::vector< std::complex<double> > Phi2 = VSH_Phi(t,p,L2,M2);
			std::complex<double> SH1 = ylm(t,p,L1,M1);
			
			//--- Add to each vector in the sum ---
			cross_product_from_VSH_coeffs_01a = vector_sum( cross_product_from_VSH_coeffs_01a, scalar_times_vector( Ar12[0] * Br12[1] * SH1, Phi2 ) );
		}
	}*/
	
	
	//--- V1 ---
	/*std::complex<double> sum_k1 = 0;
	for( int k1=0; k1<=k_max; k1++ ){
		int L1 = sqrt(k1);
		int M1 = k1 - L1*(L1+1);
		sum_k1 += coeffs_A[0][L1][L1+M1] * ylm(t,p,L1,M1);
	}
	
	std::vector< std::complex<double> > sum_k2 { 0, 0, 0 };
	for( int k2=0; k2<=k_max; k2++ ){	
		int L2 = sqrt(k2);
		int M2 = k2 - L2*(L2+1);
		sum_k2 = vector_sum( sum_k2, scalar_times_vector( coeffs_B[1][L2][L2+M2], VSH_Phi(t,p,L2,M2) ) );
	}
	
	cross_product_from_VSH_coeffs_01a = scalar_times_vector( sum_k1, sum_k2 );*/
	
	
	//--- V2 ---
	/*for( int k1=0; k1<=k_max; k1++ ){
		int L1 = sqrt(k1);
		int M1 = k1 - L1*(L1+1);
		
		std::vector< std::complex<double> > sum_over_Phi_k2 { 0, 0, 0 };
		for( int k2=k1+1; k2<=k_max; k2++ ){	
			int L2 = sqrt(k2);
			int M2 = k2 - L2*(L2+1);
			sum_over_Phi_k2 = vector_sum( sum_over_Phi_k2, scalar_times_vector( coeffs_B[1][L2][L2+M2], VSH_Phi(t,p,L2,M2) ) );
		}
		sum_over_Phi_k2 = scalar_times_vector( coeffs_A[0][L1][L1+M1] * ylm(t,p,L1,M1), sum_over_Phi_k2 );
		
		std::complex<double> sum_over_ylm_k2 = 0;
		for( int k2=k1; k2<=k_max; k2++ ){	
			int L2 = sqrt(k2);
			int M2 = k2 - L2*(L2+1);
			sum_over_ylm_k2 += coeffs_A[0][L2][L2+M2] * ylm(t,p,L2,M2);
		}
		std::vector< std::complex<double> > coeff_times_Phi_k1 = scalar_times_vector( coeffs_B[1][L1][L1+M1] * sum_over_ylm_k2, VSH_Phi(t,p,L1,M1) );
		
		std::vector< std::complex<double> > vector_for_this_k1 = vector_sum( sum_over_Phi_k2, coeff_times_Phi_k1 );
		cross_product_from_VSH_coeffs_01a = vector_sum( cross_product_from_VSH_coeffs_01a, vector_for_this_k1 );
	}*/
	
	
	//--- V3 ---
	std::vector< std::complex<double> > cross_product_from_VSH_coeffs_01a_Phi_k1_part { 0, 0, 0 };
	for( int k1=0; k1<=k_max; k1++ ){
		int L1 = sqrt(k1);
		int M1 = k1 - L1*(L1+1);
		std::complex<double> sum_over_k2 = 0;
		for( int k2=k1; k2<=k_max; k2++ ){	
			int L2 = sqrt(k2);
			int M2 = k2 - L2*(L2+1);
			sum_over_k2 += coeffs_A[0][L2][L2+M2] * ylm(t,p,L2,M2);
		}
		sum_over_k2 *= coeffs_B[1][L1][L1+M1];
		cross_product_from_VSH_coeffs_01a_Phi_k1_part = vector_sum( cross_product_from_VSH_coeffs_01a_Phi_k1_part, scalar_times_vector( sum_over_k2, VSH_Phi(t,p,L1,M1) ) );
	}
	
	std::vector< std::complex<double> > cross_product_from_VSH_coeffs_01a_Phi_k2_part { 0, 0, 0 };
	for( int k1=0; k1<=k_max; k1++ ){
		int L1 = sqrt(k1);
		int M1 = k1 - L1*(L1+1);
		std::vector< std::complex<double> > sum_over_k2 { 0, 0, 0 };
		for( int k2=k1+1; k2<=k_max; k2++ ){	
			int L2 = sqrt(k2);
			int M2 = k2 - L2*(L2+1);
			sum_over_k2 = vector_sum( sum_over_k2, scalar_times_vector( coeffs_B[1][L2][L2+M2], VSH_Phi(t,p,L2,M2) ) );
		}
		sum_over_k2 = scalar_times_vector( coeffs_A[0][L1][L1+M1] * ylm(t,p,L1,M1), sum_over_k2 );
		cross_product_from_VSH_coeffs_01a_Phi_k2_part = vector_sum( cross_product_from_VSH_coeffs_01a_Phi_k2_part, sum_over_k2 );
	}
	
	cross_product_from_VSH_coeffs_01a = vector_sum( cross_product_from_VSH_coeffs_01a, cross_product_from_VSH_coeffs_01a_Phi_k1_part );
	cross_product_from_VSH_coeffs_01a = vector_sum( cross_product_from_VSH_coeffs_01a, cross_product_from_VSH_coeffs_01a_Phi_k2_part );
		
	
	//----- Calculate remaining sums -----
	for( int k1=0; k1<=k_max; k1++ ){
		int L1 = sqrt(k1);
		int M1 = k1 - L1*(L1+1);
		for( int k2=0; k2<=k_max; k2++ ){	
			int L2 = sqrt(k2);
			int M2 = k2 - L2*(L2+1);
			//std::cout << k1 <<" ("<<L1<<","<<M1<<")\t" << k2 <<" ("<<L2<<","<<M2<<")" << std::endl;
			
			//--- Precalculate vectors and VSHs at the test point for these indices ---
			std::vector< std::complex<double> > Ar12 { coeffs_A[0][L1][L1+M1], coeffs_A[1][L1][L1+M1], coeffs_A[2][L1][L1+M1] };
			std::vector< std::complex<double> > Br12 { coeffs_B[0][L2][L2+M2], coeffs_B[1][L2][L2+M2], coeffs_B[2][L2][L2+M2] };
			std::vector< std::complex<double> > Y1   = VSH_Y  (t,p,L1,M1);
			std::vector< std::complex<double> > Y2   = VSH_Y  (t,p,L2,M2);
			std::vector< std::complex<double> > Psi1 = VSH_Psi(t,p,L1,M1);
			std::vector< std::complex<double> > Psi2 = VSH_Psi(t,p,L2,M2);
			std::vector< std::complex<double> > Phi1 = VSH_Phi(t,p,L1,M1);
			std::vector< std::complex<double> > Phi2 = VSH_Phi(t,p,L2,M2);
			std::complex<double> SH1 = ylm(t,p,L1,M1);
			std::complex<double> SH2 = ylm(t,p,L2,M2);
			
			//--- Add to each vector in the sum ---
			cross_product_from_VSH_coeffs_01 = vector_sum( cross_product_from_VSH_coeffs_01, scalar_times_vector( Ar12[0] * Br12[1] * SH1, Phi2 ) );
			cross_product_from_VSH_coeffs_02 = vector_sum( cross_product_from_VSH_coeffs_02, scalar_times_vector( Ar12[0] * Br12[2] *-SH1, Psi2 ) );
			cross_product_from_VSH_coeffs_10 = vector_sum( cross_product_from_VSH_coeffs_10, scalar_times_vector( Ar12[1] * Br12[0] *-SH2, Phi1 ) );
			cross_product_from_VSH_coeffs_11 = vector_sum( cross_product_from_VSH_coeffs_11, scalar_times_vector( Ar12[1] * Br12[1], cross_product( Psi1, Psi2 ) ) );
			cross_product_from_VSH_coeffs_12 = vector_sum( cross_product_from_VSH_coeffs_12, scalar_times_vector( Ar12[1] * Br12[2], cross_product( Psi1, Phi2 ) ) );
			cross_product_from_VSH_coeffs_20 = vector_sum( cross_product_from_VSH_coeffs_20, scalar_times_vector( Ar12[2] * Br12[0] * SH2, Psi1 ) );
			cross_product_from_VSH_coeffs_21 = vector_sum( cross_product_from_VSH_coeffs_21, scalar_times_vector( Ar12[2] * Br12[1], cross_product( Phi1, Psi2 ) ) );
			cross_product_from_VSH_coeffs_22 = vector_sum( cross_product_from_VSH_coeffs_22, scalar_times_vector( Ar12[2] * Br12[2], cross_product( Phi1, Phi2 ) ) );
		}
	}
	
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
	
	print_vector_setw( cross_product_from_VSH_coeffs_01a    , w, "Y_1   x Psi_2 (added)" );
	
	print_vector_setw( cross_product_from_VSH_coeffs_01     , w, "Y_1   x Psi_2        " );
	/*print_vector_setw( cross_product_from_VSH_coeffs_02, w, "Y_1   x Phi_2" );
	print_vector_setw( cross_product_from_VSH_coeffs_10, w, "Psi_1 x Y_2  " );
	print_vector_setw( cross_product_from_VSH_coeffs_11, w, "Psi_1 x Psi_2" );
	print_vector_setw( cross_product_from_VSH_coeffs_12, w, "Psi_1 x Phi_2" );
	print_vector_setw( cross_product_from_VSH_coeffs_20, w, "Phi_1 x Y_2  " );
	print_vector_setw( cross_product_from_VSH_coeffs_21, w, "Phi_1 x Psi_2" );
	print_vector_setw( cross_product_from_VSH_coeffs_22, w, "Phi_1 x Phi_2" );*/
	
	
	
	//----- Output execution time and finish code -----
	auto time_stop = std::chrono::high_resolution_clock::now();
	double exec_time = std::chrono::duration_cast<std::chrono::nanoseconds>( time_stop - time_start ).count() * 1e-9;
	std::cout << "\nExecution time:\t" << floor(exec_time/60) << " m " << exec_time-60*floor(exec_time/60) << " s." << std::endl;
	
	std::cout << "\nFinished" << std::endl;
}