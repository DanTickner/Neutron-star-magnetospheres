// cd OneDrive\PhD\Codes\20221108 Vector spherical harmonics\VSH series axisymmetric
// g++ Fit_VSH_series_axisymmetric_no_r_dependence_01.cpp -o Fit_VSH_series_axisymmetric_no_r_dependence_01

/*
Fit a vector spherical harmonic series to an axisymmetric (no phi-dependence) vector A(theta).
*/

#include <iostream>
#include <iomanip>
#include <chrono>
#include "../../Generate_Associated_Legendre_Functions.h"
#include "../../Generate_spherical_harmonics.h"								// Only needed to test that result agrees with VSH functions.
#include "../Headers_vector_function/Vector_function_07_axisymmetric.h"
#include "../../Numerical_Integration.h"
#include "../../Vector_Operations.h"
#include "../Headers_VSH/VSH.h"

int main(){
	
	auto time_start = std::chrono::high_resolution_clock::now();
	
	//----- Define variables -----
	int ell_max = 10;
	int N_t     = 50;
	int N_p     = 50;
	
	int i_test = 18;
	int j_test = 20;
	
	std::vector<double> t_list ( N_t );
	std::vector<double> p_list ( N_p );
	
	std::vector< std::vector< std::vector<double> > > A( 3, std::vector< std::vector<double> > ( N_t, std::vector<double> ( N_p ) ) );		// Precalculated values of the function that we wish to fit.
	std::vector< std::vector< std::vector< std::vector< std::complex<double> > > > > SH;				// Precalculated spherical harmonics
	std::vector< std::vector<double> > coeffs_VSH_series ( 3, ( std::vector<double> ( ell_max+1 ) ) );	// Coefficients of VSH series to be calculated.
	std::vector< std::vector<double> > P0( ell_max+1, ( std::vector<double> ( N_t ) ) );
	std::vector< std::vector<double> > P1( ell_max+1, ( std::vector<double> ( N_t ) ) );
	
	// Totals of absolute values of fitted coefficients in each direction. No physical interpretation but useful as a quick verification tool.
	double coeffs_total_r = 0;
	double coeffs_total_1 = 0;
	double coeffs_total_2 = 0;
	
	
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
	
	
	//----- Generate associated Legendre function coefficients -----
	generate_coeffs_associated_legendre( ell_max );
	generate_coeffs_spherical_harmonic ( ell_max );		// Only needed to test that result agrees with VSH functions.
	
	
	//----- Precalculate vector A to be fitted -----
	for( int i=0; i<N_t; i++ ){
		for( int j=0; j<N_p; j++ ){
			A[0][i][j] = A_r_function( t_list[i], p_list[j] );
			A[1][i][j] = A_t_function( t_list[i], p_list[j] );
			A[2][i][j] = A_p_function( t_list[i], p_list[j] );
		}
	}
	
	
	//----- Precalculate spherical harmonics -----
	for( int ell=0; ell<=ell_max; ell++ ){
		for( int i=0; i<N_t; i++ ){
			double x = cos( t_list[i] );
			P0[ell][i] = associated_legendre_function( x, ell, 0 );
			P1[ell][i] = associated_legendre_function( x, ell, 1 );
		}
	}
	
	
	//----- Calculate VSH series -----
	for( int ell=0; ell<=ell_max; ell++ ){
			
		//--- Define arrays to store integrand values ---
		std::vector< std::vector<double> > integrand_r( N_t, std::vector<double> ( N_p ) );
		std::vector< std::vector<double> > integrand_1( N_t, std::vector<double> ( N_p ) );
		std::vector< std::vector<double> > integrand_2( N_t, std::vector<double> ( N_p ) );
		
		//--- Calculate the VSH series integrands by their definitions ---
		for( int i=0; i<N_t; i++ ){
			for( int j=0; j<N_p; j++ ){
				integrand_r[i][j] += A[0][i][j] * P0[ell][i];
				integrand_1[i][j] += A[1][i][j] * P1[ell][i];
				integrand_2[i][j] += A[2][i][j] * P1[ell][i];
			}
		}
		
		//--- Calculate the integral. Must be done in separate loop because we look ahead to (i+1,j+1). ---
		for( int i=0; i<N_t-1; i++ ){
			for( int j=0; j<N_p-1; j++ ){
				double dA = sqrt( ( 2.0*ell+1.0 ) / ( 4.0*pi ) ) * 0.25 * ( cos(t_list[i]) - cos(t_list[i+1]) ) * ( p_list[j+1]-p_list[j] );
				coeffs_VSH_series[0][ell] += ( integrand_r[i][j] + integrand_r[i][j+1] + integrand_r[i+1][j] + integrand_r[i+1][j+1] ) * dA;
				coeffs_VSH_series[1][ell] += ( integrand_1[i][j] + integrand_1[i][j+1] + integrand_1[i+1][j] + integrand_1[i+1][j+1] ) * dA;
				coeffs_VSH_series[2][ell] += ( integrand_2[i][j] + integrand_2[i][j+1] + integrand_2[i+1][j] + integrand_2[i+1][j+1] ) * dA;
			}
		}
		
		if( ell != 0 ){
			coeffs_VSH_series[1][ell] /= (double) ell*(ell+1);
			coeffs_VSH_series[2][ell] /= (double) ell*(ell+1);
		}
	}
	
	//----- Output VSH series -----
	std::cout <<std::left<<std::setw(4)<< "ell" <<":\t"<< std::left<<std::setw(15)<< "A^r" <<"\t"<< std::left<<std::setw(15)<< "A^1" <<"\t"<< std::left<<std::setw(15)<< "A^2" << std::endl;
			
	for( int ell=0; ell<=ell_max; ell++ ){
		std::cout <<std::left<<std::setw(4)<< ell <<":\t"<< std::left<<std::setw(15)<< coeffs_VSH_series[0][ell] <<"\t"<< std::left<<std::setw(15)<< coeffs_VSH_series[1][ell] <<"\t"<< std::left<<std::setw(15)<< coeffs_VSH_series[2][ell] << std::endl;

		coeffs_total_r += coeffs_VSH_series[0][ell];
		coeffs_total_1 += coeffs_VSH_series[1][ell];
		coeffs_total_2 += coeffs_VSH_series[2][ell];
	}
	
	std::cout << "\nCoeffs total:\t" << coeffs_total_r <<"\t"<< coeffs_total_1 <<"\t"<< coeffs_total_2 <<"\t\t"<< coeffs_total_r+coeffs_total_1+coeffs_total_2 << std::endl;
	
	
	//----- Evaluate VSH series at the test point -----
	double t = t_list[i_test];
	double p = p_list[j_test];

	std::cout << "\nTest fit for (theta,phi) = (" << t <<","<< p << "):" << std::endl;
	
	/*
	std::vector< std::complex<double> > A_test_calc (3);	// Needs to be complex in order to use the predefined VSHs for verification.
	
	for( int ell=0; ell<=ell_max; ell++ ){
		A_test_calc = vector_sum( A_test_calc, scalar_times_vector( coeffs_VSH_series[0][ell], VSH_Y  (t,p,ell,0) ) );
		A_test_calc = vector_sum( A_test_calc, scalar_times_vector( coeffs_VSH_series[1][ell], VSH_Psi(t,p,ell,0) ) );
		A_test_calc = vector_sum( A_test_calc, scalar_times_vector( coeffs_VSH_series[2][ell], VSH_Phi(t,p,ell,0) ) );
	}
	
	std::vector< std::complex<double> > A_test_exact { A[0][i_test][j_test], A[1][i_test][j_test], A[2][i_test][j_test] };
	*/
	
	std::vector<double> A_test_calc (3);
	
	for( int ell=0; ell<=ell_max; ell++ ){
		double factor = sqrt( ( 2.0*ell+1.0 ) / ( 4.0*pi ) );
		A_test_calc[0] += coeffs_VSH_series[0][ell] * factor * P0[ell][i_test];
		A_test_calc[1] += coeffs_VSH_series[1][ell] * factor * P1[ell][i_test];
		A_test_calc[2] += coeffs_VSH_series[2][ell] * factor * P1[ell][i_test];
	}
	
	std::vector<double> A_test_exact { A[0][i_test][j_test], A[1][i_test][j_test], A[2][i_test][j_test] };
	
	std::cout <<std::left<<std::setw(7)<< "Cmpnt." <<std::left<<std::setw(28)<< "Calc."        <<std::left<<std::setw(28)<< "Exact"         << std::endl;
	std::cout <<std::left<<std::setw(7)<< "r"      <<std::left<<std::setw(28)<< A_test_calc[0] <<std::left<<std::setw(28)<< A_test_exact[0] << std::endl;
	std::cout <<std::left<<std::setw(7)<< "theta"  <<std::left<<std::setw(28)<< A_test_calc[1] <<std::left<<std::setw(28)<< A_test_exact[1] << std::endl;
	std::cout <<std::left<<std::setw(7)<< "phi"    <<std::left<<std::setw(28)<< A_test_calc[2] <<std::left<<std::setw(28)<< A_test_exact[2] << std::endl;
	
	return 0;
}