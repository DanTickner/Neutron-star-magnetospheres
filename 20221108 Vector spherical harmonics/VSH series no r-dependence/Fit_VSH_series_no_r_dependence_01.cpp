// cd OneDrive\PhD\Codes\20221108 Vector spherical harmonics\VSH series no r-dependence
// g++ Fit_VSH_series_no_r_dependence_01.cpp -o Fit_VSH_series_no_r_dependence_01

/*
Fit a vector spherical harmonic series to a vector A(theta,phi).
*/

#include <iostream>
#include <iomanip>
#include <chrono>
#include "../../Generate_Spherical_Harmonics.h"
#include "../Headers_vector_function/Vector_function_01.h"
#include "../Headers_VSH_series/VSH_series_calculation_01.h"
#include "../Headers_VSH_series/VSH_series_evaluation_03.h"

int main(){
	
	auto time_start = std::chrono::high_resolution_clock::now();
	
	//----- Define variables -----
	int ell_max = 5;
	int N_t     = 300;
	int N_p     = 300;
	
	int test_t_i = 28;
	int test_p_j = 42;
	
	std::vector<double> t_list ( N_t );
	std::vector<double> p_list ( N_p );
	
	std::vector< std::vector< std::vector<double> > > A( 3, std::vector< std::vector<double> > ( N_t, std::vector<double> ( N_p ) ) );		// Precalculated values of the function that we wish to fit.
	
	
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
	std::vector< std::vector< std::complex<double> > > coeffs_VSH_series_r = generate_vsh_series_Cr( A, t_list, p_list, ell_max );
	std::vector< std::vector< std::complex<double> > > coeffs_VSH_series_1 = generate_vsh_series_C1( A, t_list, p_list, ell_max );
	std::vector< std::vector< std::complex<double> > > coeffs_VSH_series_2 = generate_vsh_series_C2( A, t_list, p_list, ell_max );
	
	
	//----- Output VSH series -----
	/* Keep this layout for when I develop a way of predicting what the exact coeffs should be.
	for( int ell=0; ell<=ell_max; ell++ ){
		for( int m=-ell; m<=ell; m++ ){
			std::cout << ell <<"\t"<< m <<"\t"<< coeffs_VSH_series_r[ell][ell+m] << std::endl;
		}
	}
	
	std::cout << std::endl;
	for( int ell=0; ell<=ell_max; ell++ ){
		for( int m=-ell; m<=ell; m++ ){
			std::cout << ell <<"\t"<< m <<"\t"<< coeffs_VSH_series_1[ell][ell+m] << std::endl;
		}
	}
	
	std::cout << std::endl;
	for( int ell=0; ell<=ell_max; ell++ ){
		for( int m=-ell; m<=ell; m++ ){
			std::cout << ell <<"\t"<< m <<"\t"<< coeffs_VSH_series_2[ell][ell+m] << std::endl;
		}
	}*/
	
	std::cout <<std::left<<std::setw(4)<< "ell" <<std::left<<std::setw(4)<< "m"   <<":\t"
	<<std::left<<std::setw(15)<< "Re(C^r)" <<std::left<<std::setw(15)<< "Im(C^r)" <<"\t/\t"
	<<std::left<<std::setw(15)<< "Re(C^1)" <<std::left<<std::setw(15)<< "Im(C^1)" <<"\t/\t"
	<<std::left<<std::setw(15)<< "Re(C^2)" <<std::left<<std::setw(15)<< "Im(C^2)" << std::endl;
			
	for( int ell=0; ell<=ell_max; ell++ ){
		for( int m=-ell; m<=ell; m++ ){
			std::cout <<std::left<<std::setw(4)<< ell <<std::left<<std::setw(4)<< m <<":\t"
			<<std::left<<std::setw(15)<< coeffs_VSH_series_r[ell][ell+m].real() <<std::left<<std::setw(15)<< coeffs_VSH_series_r[ell][ell+m].imag() <<"\t/\t"
			<<std::left<<std::setw(15)<< coeffs_VSH_series_1[ell][ell+m].real() <<std::left<<std::setw(15)<< coeffs_VSH_series_1[ell][ell+m].imag() <<"\t/\t"
			<<std::left<<std::setw(15)<< coeffs_VSH_series_2[ell][ell+m].real() <<std::left<<std::setw(15)<< coeffs_VSH_series_2[ell][ell+m].imag()	<< std::endl;
		}
		std::cout << std::endl;
	}
	
	
	//----- Evaluate VSH series at the test point -----
	std::cout << "\nTest fit for (theta,phi) = (" << t_list[test_t_i] <<","<< p_list[test_p_j] << "):" << std::endl;
	
	std::vector< std::complex<double> > A_test_calc  = evaluate_vsh_series( t_list[test_t_i], p_list[test_p_j], coeffs_VSH_series_r, coeffs_VSH_series_1, coeffs_VSH_series_2 );
	std::vector< std::complex<double> > A_test_exact { A[0][test_t_i][test_p_j], A[1][test_t_i][test_p_j], A[2][test_t_i][test_p_j] };
	
	std::cout <<std::left<<std::setw(7)<< "Cmpnt." <<std::left<<std::setw(28)<< "Calc."        <<std::left<<std::setw(28)<< "Exact"         << std::endl;
	std::cout <<std::left<<std::setw(7)<< "r"      <<std::left<<std::setw(28)<< A_test_calc[0] <<std::left<<std::setw(28)<< A_test_exact[0] << std::endl;
	std::cout <<std::left<<std::setw(7)<< "theta"  <<std::left<<std::setw(28)<< A_test_calc[1] <<std::left<<std::setw(28)<< A_test_exact[1] << std::endl;
	std::cout <<std::left<<std::setw(7)<< "phi"    <<std::left<<std::setw(28)<< A_test_calc[2] <<std::left<<std::setw(28)<< A_test_exact[2] << std::endl;
	
	std::complex<double> A_test_r = 0;
	for( int ell=0; ell<=ell_max; ell++ ){
		for( int m=-ell; m<=ell; m++ ){
			A_test_r += coeffs_VSH_series_r[ell][ell+m] * ylm( t_list[test_t_i], p_list[test_p_j], ell, m );
		}
	}
	std::cout << "\nCheck different way to calculate A_r:\t"<< A_test_r << std::endl;
	
	
	//----- Calculate divergence at the test point -----
	
	
	//----- Output execution time and finish code -----
	auto time_stop = std::chrono::high_resolution_clock::now();
	double exec_time = std::chrono::duration_cast<std::chrono::nanoseconds>( time_stop - time_start ).count() * 1e-9;
	std::cout << "\nExecution time:\t" << floor(exec_time/60) << " m " << exec_time-60*floor(exec_time/60) << " s." << std::endl;
	
	std::cout << "\nFinished" << std::endl;
}