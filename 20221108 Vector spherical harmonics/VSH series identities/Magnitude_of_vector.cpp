// cd OneDrive\PhD\Codes\20221108 Vector spherical harmonics\VSH series identities
// g++ Magnitude_of_vector.cpp -o Magnitude_of_vector

/*
Fit vector spherical harmonic series to vector A(theta,phi) and calculate its magnitude.
Compare this to an analytical expression derived on paper (Book 10, p34).
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
	
	// Totals of absolute values of fitted coefficients in each direction. No physical interpretation but useful as a quick verification tool.
	std::vector< std::complex<double> > coeffs_total ( 3 );
	
	
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
		}
	}
	
	
	//----- Calculate VSH series -----
	std::vector< std::vector< std::vector< std::complex<double> > > > coeffs_A = generate_vsh_series( A, t_list, p_list, ell_max );
	
	
	//----- Output VSH series -----
	for( int ell=0; ell<=ell_max; ell++ ){
		for( int m=-ell; m<=ell; m++ ){
			for( int c=0; c<3; c++ ){
				coeffs_total[c] += std::complex<double> { abs( coeffs_A[c][ell][ell+m].real() ), abs( coeffs_A[c][ell][ell+m].imag() ) };
			}
		}
	}
	
	std::cout << "Coeffs total:\t" << coeffs_total[0] <<"\t"<< coeffs_total[1] <<"\t"<< coeffs_total[2] <<"\t\t"<< coeffs_total[0]+coeffs_total[1]+coeffs_total[2] << std::endl;
	
	
	//----- Evaluate cross product at the test point -----
	std::cout << "\nTest fit for (theta,phi) = (" << t_list[i_test] <<","<< p_list[j_test] << "):" << std::endl;
	
	std::vector< std::complex<double> > A_test_calc  = evaluate_vsh_series( t_list[i_test], p_list[j_test], coeffs_A );
	
	std::complex<double> magnitude_exact = A_test_calc[0];
	
	std::complex<double> magnitude_test = 0;
	
	for( int ell=0; ell<=ell_max; ell++ ){
		for( int m=-ell; m<=ell; m++ ){
			magnitude_test += coeffs_A[0][ell][ell+m] * ylm( t_list[i_test], p_list[j_test], ell, m );
		}
	}
	
	std::cout << "Magnitude exact expression:\t" << A[0][i_test][j_test] << std::endl;
	std::cout << "Magnitude from VSH series :\t" << magnitude_exact      << std::endl;
	std::cout << "Magnitude from VSH coeffs :\t" << magnitude_test       << std::endl;
	
	
	
	return 0;
}