/*

cd OneDrive\PhD\Codes\VSH series
g++ AVSH_series_01.cpp -o AVSH_series_01 --std=c++11
AVSH_series_01

Read a vector function 
A(r,theta) = A_r(r,theta) e_r + A_theta(r,theta) e_theta + A_phi(r,theta) e_phi
and calculate its AVSH (axisymmetric vector spherical harmonic) series coefficients at a given value of r.

Based on VSH_series_01.cpp.

*/


#include <iostream>
#include <iomanip>
#include <math.h>
#include <vector>
#include <chrono>

#include "../Generate_Associated_Legendre_Functions.h"	// Includes "../Vector_Operations.h", which we also need here.
#include "../Numerical_Integration.h"

const double pi = acos(-1);


//----- Call header file with target function definitions -----
#include "Headers/AVSH_series_dipole.h"


int main(){
	
	//----- Define variables -----
	int ell_max = 1;
	int n_steps_integration = 1000;
	int w = 14;							// Fixed width for screen output.
	
	double test_r     = 1.00900900900901;//0.321;				// Value of x at which to test the series.
	double test_theta = 0.432;
	
	
	std::vector<double> integrand_values_r ( n_steps_integration );
	std::vector<double> integrand_values_1 ( n_steps_integration );
	std::vector<double> integrand_values_2 ( n_steps_integration );
	
	std::vector<double> coeffs_r ( ell_max+1 );
	std::vector<double> coeffs_1 ( ell_max+1 );
	std::vector<double> coeffs_2 ( ell_max+1 );
	
	std::vector<double> A_r_values     ( n_steps_integration );
	std::vector<double> A_theta_values ( n_steps_integration );
	std::vector<double> A_phi_values   ( n_steps_integration );
	
	

	
	//----- Start timing -----
	auto time_start = std::chrono::high_resolution_clock::now();	
	
	//----- Build lists of coordinates and precalculate vector values -----
	double delta_theta = pi / ( (double) n_steps_integration - 1.0 );
	
	std::vector<double> theta_list ( n_steps_integration );
	
	for( int i=0; i<n_steps_integration; i++ ){
		
		theta_list[i] = (double) i * delta_theta;
		
		A_r_values    [i] = A_r_function    ( test_r, theta_list[i] );
		A_theta_values[i] = A_theta_function( test_r, theta_list[i] );
		A_phi_values  [i] = A_phi_function  ( test_r, theta_list[i] );
		
	}
	
	
	//----- Calculate VSH series coefficients -----
	
	std::cout << "ell" << "\t"<< "0" <<"\t|\t"<< std::left<<std::setw(w)<< "A^{r,ell}_m" <<std::left<<std::setw(w)<< "" <<"\t|\t"<< std::left<<std::setw(w)<< "A^{(1),ell}_m" <<std::left<<std::setw(w)<< "" <<"\t|\t"<< std::left<<std::setw(w)<< "A^{(2),ell}_m" <<std::left<<std::setw(w)<< "" << std::endl;
	
	for( int ell=0; ell<=ell_max; ell++ ){
		
		//--- Build lists of integrand values ---
		for( int i=0; i<n_steps_integration; i++ ){
				
			double theta = theta_list[i];
			
			double P_ell_0 = p_ell  ( cos(theta), ell );
			double P_ell_1 = p_ell_1( cos(theta), ell );
			
			integrand_values_r[i] = A_r_values    [i] * P_ell_0 * sin( theta );
			integrand_values_1[i] = A_theta_values[i] * P_ell_1 * sin( theta );
			integrand_values_2[i] = A_phi_values  [i] * P_ell_1 * sin( theta );
			
		}
		
		//--- Integrate and put values into array ---
		coeffs_r[ell] = trapezium( integrand_values_r, theta_list ) * sqrt( (2.0*ell+1.0)*pi );
		coeffs_1[ell] = trapezium( integrand_values_1, theta_list ) * sqrt( (2.0*ell+1.0)*pi ) / ( (double)ell*(ell+1.0) );
		coeffs_2[ell] = trapezium( integrand_values_2, theta_list ) * sqrt( (2.0*ell+1.0)*pi ) / ( (double)ell*(ell+1.0) );
		
		
		if( ell == 0 ){
			coeffs_1[0] = 0;
			coeffs_2[0] = 0;
		}
		
		std::cout << ell << "\t"<< 0 <<"\t|\t"<< std::left<<std::setw(w)<< coeffs_r[ell] <<std::left<<std::setw(w)<< A_r_ell_guess(test_r,ell) <<"\t|\t"<<std::left<<std::setw(w)<< coeffs_1[ell] <<std::left<<std::setw(w)<< A_1_ell_guess(test_r,ell) <<"\t|\t"<< std::left<<std::setw(w)<< coeffs_2[ell] <<std::left<<std::setw(w)<< A_2_ell_guess(test_r,ell) << std::endl;
	}
	
	
	
	
	//----- Evaluate VSH series at the chosen gridpoint -----
	
	std::vector< std::complex<double> > A_exact { A_r_function(test_r,test_theta), A_theta_function(test_r,test_theta), A_phi_function(test_r,test_theta) };
	
	std::vector< std::complex<double> > A_series { 0, 0, 0 };
	
	for( int ell=0; ell<=ell_max; ell++ ){
		
		double sqrt_factor = sqrt( (2.0*ell+1.0) / ( 4.0 * pi ) );
		
		double P_ell_0 = p_ell  ( cos(test_theta), ell );
		double P_ell_1 = p_ell_1( cos(test_theta), ell );
		
		A_series[0] += sqrt_factor * P_ell_0 * coeffs_r[ell];
		A_series[1] += sqrt_factor * P_ell_1 * coeffs_1[ell];
		A_series[2] += sqrt_factor * P_ell_1 * coeffs_2[ell];
		
	}
	
	std::cout << std::endl;
	cout_vector_setw( A_exact , w, "A_exact " );
	cout_vector_setw( A_series, w, "A_series" );
	
	
	
	
	//----- Output execution time and finish code -----
	auto   time_stop   = std::chrono::high_resolution_clock::now();
	double exec_time   = std::chrono::duration_cast<std::chrono::nanoseconds>( time_stop - time_start ).count() * 1e-9;
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
	
	
	std::cout << "done" << std::endl;
	return 0;
}