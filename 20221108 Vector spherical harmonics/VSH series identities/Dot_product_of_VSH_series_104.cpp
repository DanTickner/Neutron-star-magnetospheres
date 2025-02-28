// cd OneDrive\PhD\Codes\20221108 Vector spherical harmonics\VSH series identities
// g++ Dot_product_of_VSH_series_104.cpp -o Dot_product_of_VSH_series_104

/*
Fit vector spherical harmonic series to two vectors A(theta,phi) and B(theta,phi) and calculate their dot product.
Compare this to an analytical expression derived on paper (Book 11, p4).
Now that only the nonzero terms are considered, no need to store them in separate variables.
No longer precalculate the VSHs or rely on functions for dot products.
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
	int ell_max = 12;
	int N_t     = 100;
	int N_p     = 100;
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
			
			B[0][i][j] = A_r_function( t_list[i], p_list[j] );
			B[1][i][j] = A_t_function( t_list[i], p_list[j] );
			B[2][i][j] = A_p_function( t_list[i], p_list[j] );
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
	
	
	//----- Evaluate dot product at the test point -----
	double t = t_list[i_test];
	double p = p_list[j_test];
	
	std::cout << "\nTest fit for (theta,phi) = (" << t <<","<< p << "):" << std::endl;
	
	std::vector< std::complex<double> > A_test_calc  = evaluate_vsh_series( t, p, coeffs_A );
	std::vector< std::complex<double> > B_test_calc  = evaluate_vsh_series( t, p, coeffs_B );
	
	std::vector< std::complex<double> > A_test_exact{ A[0][i_test][j_test], A[1][i_test][j_test], A[2][i_test][j_test] };
	std::vector< std::complex<double> > B_test_exact{ B[0][i_test][j_test], B[1][i_test][j_test], B[2][i_test][j_test] };
	
	std::complex<double> dot_product_exact1           = A_test_exact[0]*B_test_exact[0] + A_test_exact[1]*B_test_exact[1] + A_test_exact[2]*B_test_exact[2];
	std::complex<double> dot_product_from_VSH_series1 = A_test_calc [0]*B_test_calc [0] + A_test_calc [1]*B_test_calc [1] + A_test_calc [2]*B_test_calc [2];
	
	print_vector( A_test_exact, "Vector A (exact)" );
	print_vector( B_test_exact, "Vector B (exact)" );
	std::cout << "Dot product exact          :\t" << dot_product_exact1           << std::endl;
	std::cout << "Dot product from VSH series:\t" << dot_product_from_VSH_series1 << std::endl;
	
	
	//----- Calculate dot product from VSH coefficients -----
	int k_max = pow(ell_max+1,2) - 1;
	
	std::complex<double> dot_product_from_VSH_coeffs = 0;
	
	std::complex<double> dot_product_from_VSH_coeffs_a = 0;
	std::complex<double> dot_product_from_VSH_coeffs_b = 0;
	std::complex<double> dot_product_from_VSH_coeffs_c = 0;
	
	//--- V1 (keep precalculated vectors Ar12 Br12) ---
	/*for( int k1=0; k1<=k_max; k1++ ){
		int L1 = sqrt(k1);
		int M1 = k1 - L1*(L1+1);
		std::complex<double> SH_1               = ylm              (t,p,L1,M1);
		std::complex<double> SH_1_over_sintheta = ylm_over_sintheta(t,p,L1,M1);
		std::complex<double> SH_1_dtheta        = ylm_dtheta       (t,p,L1,M1);
		
		for( int k2=0; k2<=k_max; k2++ ){	
			int L2 = sqrt(k2);
			int M2 = k2 - L2*(L2+1);
			std::complex<double> SH_2 = ylm(t,p,L2,M2);
			std::complex<double> SH_2_over_sintheta  = ylm_over_sintheta(t,p,L2,M2);
			std::complex<double> SH_2_dtheta         = ylm_dtheta       (t,p,L2,M2);
			
			//--- Precalculate vectors and VSHs at the test point for these indices ---
			std::vector< std::complex<double> > Ar12 { coeffs_A[0][L1][L1+M1], coeffs_A[1][L1][L1+M1], coeffs_A[2][L1][L1+M1] };
			std::vector< std::complex<double> > Br12 { coeffs_B[0][L2][L2+M2], coeffs_B[1][L2][L2+M2], coeffs_B[2][L2][L2+M2] };
			
			//--- Add to each vector in the sum ---
			std::complex<double> Y1_dot_Y1 = SH_1 * SH_2;
			
			std::complex<double> Psi1_dot_Psi2 = SH_1_dtheta * SH_2_dtheta;
			if( ( M1 != 0 ) && ( M2 != 0 ) ){
				Psi1_dot_Psi2 -= (double)M1*M2 * SH_1_over_sintheta * SH_2_over_sintheta;
			}
			
			std::complex<double> Psi1_dot_Phi2 = 0;
			if( M1 != 0 ){
				Psi1_dot_Phi2 += std::complex<double>{0,(double)M1} * SH_1_over_sintheta * SH_2_dtheta;
			}
			if( M2 != 0 ){
				Psi1_dot_Phi2 -= std::complex<double>{0,(double)M2} * SH_2_over_sintheta * SH_1_dtheta;
			}
			
			dot_product_from_VSH_coeffs += Ar12[0] * Br12[0] * Y1_dot_Y1;
			dot_product_from_VSH_coeffs += ( Ar12[1] * Br12[1] + Ar12[2] * Br12[2] ) * Psi1_dot_Psi2;
			dot_product_from_VSH_coeffs += ( Ar12[1] * Br12[2] - Ar12[2] * Br12[1] ) * Psi1_dot_Phi2;
			
		}
	}*/
	
	for( int k1=0; k1<=k_max; k1++ ){
		int L1 = sqrt(k1);
		int M1 = k1 - L1*(L1+1);
		std::complex<double> SH_1               = ylm              (t,p,L1,M1);
		std::complex<double> SH_1_over_sintheta = ylm_over_sintheta(t,p,L1,M1);
		std::complex<double> SH_1_dtheta        = ylm_dtheta       (t,p,L1,M1);
		
		for( int k2=0; k2<=k_max; k2++ ){	
			int L2 = sqrt(k2);
			int M2 = k2 - L2*(L2+1);
			std::complex<double> SH_2 = ylm(t,p,L2,M2);
			std::complex<double> SH_2_over_sintheta  = ylm_over_sintheta(t,p,L2,M2);
			std::complex<double> SH_2_dtheta         = ylm_dtheta       (t,p,L2,M2);
			
			//--- Precalculate vectors and VSHs at the test point for these indices ---
			std::vector< std::complex<double> > Ar12 { coeffs_A[0][L1][L1+M1], coeffs_A[1][L1][L1+M1], coeffs_A[2][L1][L1+M1] };
			std::vector< std::complex<double> > Br12 { coeffs_B[0][L2][L2+M2], coeffs_B[1][L2][L2+M2], coeffs_B[2][L2][L2+M2] };
			
			//--- Add to each vector in the sum ---
			std::complex<double> Y1_dot_Y1 = SH_1 * SH_2;
			
			std::complex<double> Psi1_dot_Psi2 = SH_1_dtheta * SH_2_dtheta;
			if( ( M1 != 0 ) && ( M2 != 0 ) ){
				Psi1_dot_Psi2 -= (double)M1*M2 * SH_1_over_sintheta * SH_2_over_sintheta;
			}
			
			std::complex<double> Psi1_dot_Phi2 = 0;
			if( M1 != 0 ){
				Psi1_dot_Phi2 += std::complex<double>{0,(double)M1} * SH_1_over_sintheta * SH_2_dtheta;
			}
			if( M2 != 0 ){
				Psi1_dot_Phi2 -= std::complex<double>{0,(double)M2} * SH_2_over_sintheta * SH_1_dtheta;
			}
			
			dot_product_from_VSH_coeffs += Y1_dot_Y1     *   coeffs_A[0][L1][L1+M1]*coeffs_B[0][L2][L2+M2];
			dot_product_from_VSH_coeffs += Psi1_dot_Psi2 * ( coeffs_A[1][L1][L1+M1]*coeffs_B[1][L2][L2+M2] + coeffs_A[2][L1][L1+M1]*coeffs_B[2][L2][L2+M2] );
			dot_product_from_VSH_coeffs += Psi1_dot_Phi2 * ( coeffs_A[1][L1][L1+M1]*coeffs_B[2][L2][L2+M2] - coeffs_A[2][L1][L1+M1]*coeffs_B[1][L2][L2+M2] );
			
			//--- Extra part that adds to the three terms separately for future reference. Not needed for the calculation. ---
			dot_product_from_VSH_coeffs_a += Y1_dot_Y1     *   coeffs_A[0][L1][L1+M1]*coeffs_B[0][L2][L2+M2];
			dot_product_from_VSH_coeffs_b += Psi1_dot_Psi2 * ( coeffs_A[1][L1][L1+M1]*coeffs_B[1][L2][L2+M2] + coeffs_A[2][L1][L1+M1]*coeffs_B[2][L2][L2+M2] );
			dot_product_from_VSH_coeffs_c += Psi1_dot_Phi2 * ( coeffs_A[1][L1][L1+M1]*coeffs_B[2][L2][L2+M2] - coeffs_A[2][L1][L1+M1]*coeffs_B[1][L2][L2+M2] );
			
		}
	}
	
	std::cout << "Dot product from VSH coeffs:\t" << dot_product_from_VSH_coeffs << "\n\n" << std::endl;
	std::cout << "Contribution from each term:\t" << dot_product_from_VSH_coeffs_a <<"\t"<< dot_product_from_VSH_coeffs_b <<"\t"<< dot_product_from_VSH_coeffs_c << std::endl;
	
	
	
	//----- Output execution time and finish code -----
	auto time_stop = std::chrono::high_resolution_clock::now();
	double exec_time = std::chrono::duration_cast<std::chrono::nanoseconds>( time_stop - time_start ).count() * 1e-9;
	std::cout << "\nExecution time:\t" << floor(exec_time/60) << " m " << exec_time-60*floor(exec_time/60) << " s." << std::endl;
	
	std::cout << "\nFinished" << std::endl;
}