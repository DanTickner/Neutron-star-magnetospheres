/*

cd OneDrive\PhD\Codes\VSH series
g++ VSH_series_01.cpp -o VSH_series_01 --std=c++11
VSH_series_01

Read a vector function 
A(r,theta,phi) = A_r(r,theta,phi) e_r + A_theta(r,theta,phi) e_theta + A_phi(r,theta,phi) e_phi
and calculate its VSH series coefficients at a given value of r.

Inefficient due to not precalculating any of the VSHs or storing the vector A in an array, but good enough for prototyping.
Using 10^2  terms takes on the order of seconds and is accurate within around 5%.
Using 100^2 terms takes on the order of minutes and is accurate within around 0.5%.

*/


#include <iostream>
#include <iomanip>
#include <math.h>
#include <vector>
#include <chrono>

#include "../Vector_Operations.h"
#include "../Generate_Spherical_Harmonics.h"
#include "../Numerical_Integration.h"


//----- Call header file with target function definitions -----
#include "Headers/20230726_VSH_series_test_02.h"


int main(){
	
	//----- Define variables -----
	int ell_max = 3;
	int n_steps_integration_theta = 3e2;
	int n_steps_integration_phi   = 3e2;
	int w = 35;							// Fixed width for screen output.
	
	double test_r     = 0.321;				// Value of x at which to test the series.
	double test_theta = 0.432;
	double test_phi   = 0.567;
	
	
	std::vector< std::vector< std::complex<double> > > integrand_values_Y   ( n_steps_integration_theta, std::vector< std::complex<double> > ( n_steps_integration_phi ) );
	std::vector< std::vector< std::complex<double> > > integrand_values_Psi ( n_steps_integration_theta, std::vector< std::complex<double> > ( n_steps_integration_phi ) );
	std::vector< std::vector< std::complex<double> > > integrand_values_Phi ( n_steps_integration_theta, std::vector< std::complex<double> > ( n_steps_integration_phi ) );
	
	std::vector< std::vector< std::complex<double> > > coeffs_r ( ell_max+1,std::vector< std::complex<double> > ( 2*ell_max+1 ) );
	std::vector< std::vector< std::complex<double> > > coeffs_1 ( ell_max+1,std::vector< std::complex<double> > ( 2*ell_max+1 ) );
	std::vector< std::vector< std::complex<double> > > coeffs_2 ( ell_max+1,std::vector< std::complex<double> > ( 2*ell_max+1 ) );
	
	
	//----- Generate coefficients for spherical harmonics -----
	// Need ell_max+1 because dYlm/dtheta * 1/sin(theta) is defined in terms of Y_{ell+1}^{m\pm1}.
	generate_coeffs_spherical_harmonic( ell_max + 1 );
	
	//----- Build lists of coordinates -----
	double delta_theta =       pi / ( (double) n_steps_integration_theta - 1.0 );
	double delta_phi   = 2.0 * pi / ( (double) n_steps_integration_phi   - 1.0 );
	
	std::vector<double> theta_list ( n_steps_integration_theta );
	std::vector<double> phi_list   ( n_steps_integration_phi   );
	
	for( int i=0; i<n_steps_integration_theta; i++ ){
		theta_list[i] = (double) i * delta_theta;
	}
	
	for( int j=0; j<n_steps_integration_phi; j++ ){
		phi_list[j] = (double) j * delta_phi;
	}
	
	
	//----- Calculate VSH series coefficients -----
	auto time_start = std::chrono::high_resolution_clock::now();
	
	//----- ell=0 case separately, as it only applies to the r-coefficients -----
	std::cout << "ell" << "\t"<< "m" <<"\t|\t"<< std::left<<std::setw(w)<< "A^{r,ell}_m" <<std::left<<std::setw(w)<< "A^{(1),ell}_m" <<std::left<<std::setw(w)<< "A^{(2),ell}_m" << std::endl;
	
	for( int i=0; i<n_steps_integration_theta; i++ ){
		for( int j=0; j<n_steps_integration_phi; j++ ){
			
			double theta = theta_list[i];
			double phi   = phi_list  [j];
			
			std::vector< std::complex<double> > A { A_r_function(test_r,theta,phi), A_theta_function(test_r,theta,phi), A_phi_function(test_r,theta,phi) };
			
			integrand_values_Y  [i][j] = dot_product( A, vector_conj( VSH_Y  (theta,phi,0,0) ) );
			
		}
	}
	
	coeffs_r[0][0+0] = integral_spherical_surface( integrand_values_Y, theta_list, phi_list );
	
	std::cout << 0 << "\t"<< 0 <<"\t|\t"<< std::left<<std::setw(w)<< coeffs_r[0][0+0] <<std::left<<std::setw(w)<< coeffs_1[0][0+0] <<std::left<<std::setw(w)<< coeffs_2[0][0+0] <<"\n"<< std::endl;
	
	
	//----- ell>0 cases -----
	for( int ell=1; ell<=ell_max; ell++ ){
		for( int m=-ell; m<=ell; m++ ){
		
			//--- Build lists of integrand values ---
			for( int i=0; i<n_steps_integration_theta; i++ ){
				for( int j=0; j<n_steps_integration_phi; j++ ){
					
					double theta = theta_list[i];
					double phi   = phi_list  [j];
					
					std::vector< std::complex<double> > A { A_r_function(test_r,theta,phi), A_theta_function(test_r,theta,phi), A_phi_function(test_r,theta,phi) };
					
					integrand_values_Y  [i][j] = dot_product( A, vector_conj( VSH_Y  (theta,phi,ell,m) ) );
					integrand_values_Psi[i][j] = dot_product( A, vector_conj( VSH_Psi(theta,phi,ell,m) ) ) / ( (double)ell*(ell+1.0) );
					integrand_values_Phi[i][j] = dot_product( A, vector_conj( VSH_Phi(theta,phi,ell,m) ) ) / ( (double)ell*(ell+1.0) );
					
				}
			}
			
			//--- Integrate and put values into array ---
			// Note that we offset the m-index by ell to accommodate for negative values.
			coeffs_r[ell][ell+m] = integral_spherical_surface( integrand_values_Y  , theta_list, phi_list );
			coeffs_1[ell][ell+m] = integral_spherical_surface( integrand_values_Psi, theta_list, phi_list );
			coeffs_2[ell][ell+m] = integral_spherical_surface( integrand_values_Phi, theta_list, phi_list );
			
			std::cout << ell << "\t"<< m <<"\t|\t"<< std::left<<std::setw(w)<< coeffs_r[ell][ell+m] <<std::left<<std::setw(w)<< coeffs_1[ell][ell+m] <<std::left<<std::setw(w)<< coeffs_2[ell][ell+m] << std::endl;
		}
		std::cout << std::endl;
	}
	
	
	
	
	//----- Evaluate VSH series at the chosen gridpoint -----
	std::vector< std::complex<double> > A_exact { A_r_function(test_r,test_theta,test_phi), A_theta_function(test_r,test_theta,test_phi), A_phi_function(test_r,test_theta,test_phi) };
	
	std::vector< std::complex<double> > A_series { 0, 0, 0 };
	
	for( int ell=0; ell<=ell_max; ell++ ){
		for( int m=-ell; m<=ell; m++ ){
			A_series = vector_sum( A_series, scalar_times_vector( coeffs_r[ell][ell+m], VSH_Y  (test_theta,test_phi,ell,m) ) );
			A_series = vector_sum( A_series, scalar_times_vector( coeffs_1[ell][ell+m], VSH_Psi(test_theta,test_phi,ell,m) ) );
			A_series = vector_sum( A_series, scalar_times_vector( coeffs_2[ell][ell+m], VSH_Phi(test_theta,test_phi,ell,m) ) );
		}
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