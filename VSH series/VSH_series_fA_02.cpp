/*

cd OneDrive\PhD\Codes\VSH series
g++ VSH_series_fA_02.cpp -o VSH_series_fA_02 --std=c++11
VSH_series_fA_02

Read a vector function 
A(r,theta,phi) = A_r(r,theta,phi) e_r + A_theta(r,theta,phi) e_theta + A_phi(r,theta,phi) e_phi
and a scalar function
f(r,theta,phi)
and calculate the VSH series coefficients of fA at a given value of r.

V02: Get rid of separate ell=0 consideration, similarly to VSH_series_03.cpp vs VSH_series_02.cpp.

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
#include "Headers/20230727_VSH_series_test_fA_01.h"


int main(){
	
	//----- Define variables -----
	int ell_max = 7;						// Remember to allow for higher ell_max than you would need for just A alone, because f effectively makes the vector different and so it requires a different amount of ells.
	int n_steps_integration_theta = 5e1;
	int n_steps_integration_phi   = 5e1;
	int w = 35;							// Fixed width for screen output.
	
	double test_r     = 0.321;				// Value of x at which to test the series.
	double test_theta = 0.432;
	double test_phi   = 0.567;
	
	
	std::vector< std::vector< std::complex<double> > > integrand_values_A_r  ( n_steps_integration_theta, std::vector< std::complex<double> > ( n_steps_integration_phi ) );
	std::vector< std::vector< std::complex<double> > > integrand_values_A_1  ( n_steps_integration_theta, std::vector< std::complex<double> > ( n_steps_integration_phi ) );
	std::vector< std::vector< std::complex<double> > > integrand_values_A_2  ( n_steps_integration_theta, std::vector< std::complex<double> > ( n_steps_integration_phi ) );
	std::vector< std::vector< std::complex<double> > > integrand_values_fA_r ( n_steps_integration_theta, std::vector< std::complex<double> > ( n_steps_integration_phi ) );	// Dual purpose for numerical and series calculations
	std::vector< std::vector< std::complex<double> > > integrand_values_fA_1 ( n_steps_integration_theta, std::vector< std::complex<double> > ( n_steps_integration_phi ) );
	std::vector< std::vector< std::complex<double> > > integrand_values_fA_2 ( n_steps_integration_theta, std::vector< std::complex<double> > ( n_steps_integration_phi ) );
	
	std::vector< std::vector< std::complex<double> > > coeffs_A_r ( ell_max+1,std::vector< std::complex<double> > ( 2*ell_max+1 ) );
	std::vector< std::vector< std::complex<double> > > coeffs_A_1 ( ell_max+1,std::vector< std::complex<double> > ( 2*ell_max+1 ) );
	std::vector< std::vector< std::complex<double> > > coeffs_A_2 ( ell_max+1,std::vector< std::complex<double> > ( 2*ell_max+1 ) );
	
	std::vector< std::vector< std::complex<double> > > coeffs_fA_r_numerical ( ell_max+1,std::vector< std::complex<double> > ( 2*ell_max+1 ) );	// For calculation of the VSH coeffs of fA the "normal numerical way"
	std::vector< std::vector< std::complex<double> > > coeffs_fA_1_numerical ( ell_max+1,std::vector< std::complex<double> > ( 2*ell_max+1 ) );
	std::vector< std::vector< std::complex<double> > > coeffs_fA_2_numerical ( ell_max+1,std::vector< std::complex<double> > ( 2*ell_max+1 ) );
	
	std::vector< std::vector< std::complex<double> > > coeffs_fA_r_series ( ell_max+1,std::vector< std::complex<double> > ( 2*ell_max+1 ) );	// For calculation of the VSH coeffs of fA in terms of the series coeffs of A times some integrals.
	std::vector< std::vector< std::complex<double> > > coeffs_fA_1_series ( ell_max+1,std::vector< std::complex<double> > ( 2*ell_max+1 ) );
	std::vector< std::vector< std::complex<double> > > coeffs_fA_2_series ( ell_max+1,std::vector< std::complex<double> > ( 2*ell_max+1 ) );
	
	std::vector< std::vector< std::complex<double> > > A_r_values      ( n_steps_integration_theta, std::vector< std::complex<double> > ( n_steps_integration_phi ) );
	std::vector< std::vector< std::complex<double> > > A_theta_values  ( n_steps_integration_theta, std::vector< std::complex<double> > ( n_steps_integration_phi ) );
	std::vector< std::vector< std::complex<double> > > A_phi_values    ( n_steps_integration_theta, std::vector< std::complex<double> > ( n_steps_integration_phi ) );
	std::vector< std::vector< std::complex<double> > > f_values        ( n_steps_integration_theta, std::vector< std::complex<double> > ( n_steps_integration_phi ) );
	std::vector< std::vector< std::complex<double> > > fA_r_values     ( n_steps_integration_theta, std::vector< std::complex<double> > ( n_steps_integration_phi ) );
	std::vector< std::vector< std::complex<double> > > fA_theta_values ( n_steps_integration_theta, std::vector< std::complex<double> > ( n_steps_integration_phi ) );
	std::vector< std::vector< std::complex<double> > > fA_phi_values   ( n_steps_integration_theta, std::vector< std::complex<double> > ( n_steps_integration_phi ) );
	
	std::vector< std::vector< std::vector< std::vector< std::complex<double> > > > > ylm_values                    ( n_steps_integration_theta, std::vector< std::vector< std::vector< std::complex<double> > > > ( n_steps_integration_phi, std::vector< std::vector< std::complex<double> > > ( ell_max+1, std::vector< std::complex<double> > ( 2*ell_max+1 ) ) ) );
	std::vector< std::vector< std::vector< std::vector< std::complex<double> > > > > ylm_dtheta_values             ( n_steps_integration_theta, std::vector< std::vector< std::vector< std::complex<double> > > > ( n_steps_integration_phi, std::vector< std::vector< std::complex<double> > > ( ell_max+1, std::vector< std::complex<double> > ( 2*ell_max+1 ) ) ) );
	std::vector< std::vector< std::vector< std::vector< std::complex<double> > > > > ylm_dphi_over_sintheta_values ( n_steps_integration_theta, std::vector< std::vector< std::vector< std::complex<double> > > > ( n_steps_integration_phi, std::vector< std::vector< std::complex<double> > > ( ell_max+1, std::vector< std::complex<double> > ( 2*ell_max+1 ) ) ) );

	
	//----- Start timing -----
	auto time_start = std::chrono::high_resolution_clock::now();
	
	//----- Generate coefficients for spherical harmonics -----
	// Need ell_max+1 because dYlm/dtheta * 1/sin(theta) is defined in terms of Y_{ell+1}^{m\pm1}.
	generate_coeffs_spherical_harmonic( ell_max + 1 );
	
	//----- Build lists of coordinates and precalculate vector values -----
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
	
	for( int i=0; i<n_steps_integration_theta; i++ ){
		for( int j=0; j<n_steps_integration_phi; j++ ){
			
			double theta = theta_list[i];
			double phi   = phi_list  [j];
			
			A_r_values    [i][j] = A_r_function    ( test_r, theta, phi );
			A_theta_values[i][j] = A_theta_function( test_r, theta, phi );
			A_phi_values  [i][j] = A_phi_function  ( test_r, theta, phi );
			
			f_values[i][j] = f_function( test_r, theta, phi );
			
			fA_r_values    [i][j] = f_values[i][j] * A_r_values    [i][j];
			fA_theta_values[i][j] = f_values[i][j] * A_theta_values[i][j];
			fA_phi_values  [i][j] = f_values[i][j] * A_phi_values  [i][j];
			
		}
	}
			
		
	
	
	//----- Calculate VSH series coefficients -----
	
	std::cout << "ell" << "\t"<< "m" <<"\t|\t"<< std::left<<std::setw(w)<< "A^{r,ell}_m" <<std::left<<std::setw(w)<< "A^{(1),ell}_m" <<std::left<<std::setw(w)<< "A^{(2),ell}_m" << std::endl;
	
	for( int ell=0; ell<=ell_max; ell++ ){
		for( int m=-ell; m<=ell; m++ ){
		
			//--- Build lists of integrand values ---
			for( int i=0; i<n_steps_integration_theta; i++ ){
				for( int j=0; j<n_steps_integration_phi; j++ ){
					
					double theta = theta_list[i];
					double phi   = phi_list  [j];
					
					std::complex<double> ylm_conj                    = std::conj( ylm                   ( theta, phi, ell, m ) );
					std::complex<double> ylm_dtheta_conj             = std::conj( ylm_dtheta            ( theta, phi, ell, m ) );
					std::complex<double> ylm_dphi_over_sintheta_conj = std::conj( ylm_dphi_over_sintheta( theta, phi, ell, m ) );
					
					integrand_values_A_r[i][j] = A_r_values[i][j] * ylm_conj;
					integrand_values_A_1[i][j] = ( A_theta_values[i][j] *  ylm_dtheta_conj             + A_phi_values[i][j] * ylm_dphi_over_sintheta_conj ) / ( (double)ell*(ell+1.0) );
					integrand_values_A_2[i][j] = ( A_theta_values[i][j] * -ylm_dphi_over_sintheta_conj + A_phi_values[i][j] * ylm_dtheta_conj             ) / ( (double)ell*(ell+1.0) );
					
					integrand_values_fA_r[i][j] = fA_r_values[i][j] * ylm_conj;
					integrand_values_fA_1[i][j] = ( fA_theta_values[i][j] *  ylm_dtheta_conj             + fA_phi_values[i][j] * ylm_dphi_over_sintheta_conj ) / ( (double)ell*(ell+1.0) );
					integrand_values_fA_2[i][j] = ( fA_theta_values[i][j] * -ylm_dphi_over_sintheta_conj + fA_phi_values[i][j] * ylm_dtheta_conj             ) / ( (double)ell*(ell+1.0) );
					
				}
			}
			
			//--- Integrate and put values into array ---
			// Note that we offset the m-index by ell to accommodate for negative values.
			coeffs_A_r[ell][ell+m] = integral_spherical_surface( integrand_values_A_r, theta_list, phi_list );
			coeffs_A_1[ell][ell+m] = integral_spherical_surface( integrand_values_A_1, theta_list, phi_list );
			coeffs_A_2[ell][ell+m] = integral_spherical_surface( integrand_values_A_2, theta_list, phi_list );
			
			coeffs_fA_r_numerical[ell][ell+m] = integral_spherical_surface( integrand_values_fA_r, theta_list, phi_list );
			coeffs_fA_1_numerical[ell][ell+m] = integral_spherical_surface( integrand_values_fA_1, theta_list, phi_list );
			coeffs_fA_2_numerical[ell][ell+m] = integral_spherical_surface( integrand_values_fA_2, theta_list, phi_list );
			
			if( ell == 0 ){
				coeffs_A_1           [0][0] = 0;
				coeffs_A_2           [0][0] = 0;
				coeffs_fA_1_numerical[0][0] = 0;
				coeffs_fA_2_numerical[0][0] = 0;
			}
			
			std::cout << ell << "\t"<< m <<"\t|\t"<< std::left<<std::setw(w)<< coeffs_A_r[ell][ell+m] <<std::left<<std::setw(w)<< coeffs_A_1[ell][ell+m] <<std::left<<std::setw(w)<< coeffs_A_2[ell][ell+m] << std::endl;
		}
		std::cout << std::endl;
	}
	
	
	
	
	//----- Calculate fA coefficients by integration -----
	
	//--- Precalculate spherical harmonics and function values ---
	for( int i=0; i<n_steps_integration_theta; i++ ){
		for( int j=0; j<n_steps_integration_phi; j++ ){
			
			double theta = theta_list[i];
			double phi   = phi_list  [j];
			
			for( int ell=0; ell<=ell_max; ell++ ){
				for( int m=-ell; m<=ell; m++ ){
					
					ylm_values                   [i][j][ell][ell+m] = ylm                   ( theta, phi, ell, m );
					ylm_dtheta_values            [i][j][ell][ell+m] = ylm_dtheta            ( theta, phi, ell, m );
					ylm_dphi_over_sintheta_values[i][j][ell][ell+m] = ylm_dphi_over_sintheta( theta, phi, ell, m );
					
				}
			}
		}
	}
	
	
	//--- Perform integration ---
	for( int ell=0; ell<=ell_max; ell++ ){
		for( int m=-ell; m<=ell; m++ ){
			for( int L1=0; L1<=ell_max; L1++ ){
				for( int M1=-L1; M1<=L1; M1++ ){
				
					//--- Build integrand lists ---
					for( int i=0; i<n_steps_integration_theta; i++ ){
						for( int j=0; j<n_steps_integration_phi; j++ ){
							
							double theta = theta_list[i];
							double phi   = phi_list  [j];
							
							std::complex<double> f        = f_values[i][j];
							std::complex<double> ylm_conj = std::conj( ylm_values[i][j][ell][ell+m] );
							std::complex<double> ylm_1    = ylm_values[i][j][L1][L1+M1];
							
							std::complex<double> ylm_dtheta_conj             = std::conj( ylm_dtheta_values            [i][j][ell][ell+m] );
							std::complex<double> ylm_dphi_over_sintheta_conj = std::conj( ylm_dphi_over_sintheta_values[i][j][ell][ell+m] );
							
							std::complex<double> ylm_1_dtheta             = ylm_dtheta_values            [i][j][L1][L1+M1];
							std::complex<double> ylm_1_dphi_over_sintheta = ylm_dphi_over_sintheta_values[i][j][L1][L1+M1];
							
							integrand_values_fA_r[i][j] = coeffs_A_r[L1][L1+M1] * f * ylm_conj * ylm_1;
							
							integrand_values_fA_1[i][j] =   coeffs_A_1[L1][L1+M1] * f * ( ylm_1_dphi_over_sintheta * ylm_dphi_over_sintheta_conj + ylm_1_dtheta * ylm_dtheta_conj )
							                              + coeffs_A_2[L1][L1+M1] * f * ( ylm_1_dtheta * ylm_dphi_over_sintheta_conj - ylm_1_dphi_over_sintheta * ylm_dtheta_conj );
							
							integrand_values_fA_2[i][j] =   coeffs_A_2[L1][L1+M1] * f * ( ylm_1_dphi_over_sintheta * ylm_dphi_over_sintheta_conj + ylm_1_dtheta * ylm_dtheta_conj )
							                              + coeffs_A_1[L1][L1+M1] * f * ( ylm_1_dphi_over_sintheta * ylm_dtheta_conj - ylm_1_dtheta * ylm_dphi_over_sintheta_conj );
						}
					}
					
					//--- Evaluate integrals ---
					coeffs_fA_r_series[ell][ell+m] += integral_spherical_surface( integrand_values_fA_r, theta_list, phi_list );
					if( ell > 0 ){
						coeffs_fA_1_series[ell][ell+m] += integral_spherical_surface( integrand_values_fA_1, theta_list, phi_list ) / ( (double)ell*(ell+1.0) );
						coeffs_fA_2_series[ell][ell+m] += integral_spherical_surface( integrand_values_fA_2, theta_list, phi_list ) / ( (double)ell*(ell+1.0) );
					}
				}
			}
			
			//std::cout << coeffs_fA_r_numerical[ell][ell+m] <<"\t"<< coeffs_fA_r_series[ell][ell+m] << std::endl;
		}
	}
	
	
	
	
	
	
	//----- Evaluate VSH series at the chosen gridpoint -----
	std::vector< std::complex<double> > A_exact  { A_r_function(test_r,test_theta,test_phi), A_theta_function(test_r,test_theta,test_phi), A_phi_function(test_r,test_theta,test_phi) };
	std::vector< std::complex<double> > fA_exact { f_function(test_r,test_theta,test_phi)*A_exact[0], f_function(test_r,test_theta,test_phi)*A_exact[1], f_function(test_r,test_theta,test_phi)*A_exact[2] };
	
	std::vector< std::complex<double> > A_series     { 0, 0, 0 };
	std::vector< std::complex<double> > fA_numerical { 0, 0, 0 };
	std::vector< std::complex<double> > fA_series    { 0, 0, 0 };
	
	for( int ell=0; ell<=ell_max; ell++ ){
		for( int m=-ell; m<=ell; m++ ){
			
			A_series = vector_sum( A_series, scalar_times_vector( coeffs_A_r[ell][ell+m], VSH_Y  (test_theta,test_phi,ell,m) ) );
			A_series = vector_sum( A_series, scalar_times_vector( coeffs_A_1[ell][ell+m], VSH_Psi(test_theta,test_phi,ell,m) ) );
			A_series = vector_sum( A_series, scalar_times_vector( coeffs_A_2[ell][ell+m], VSH_Phi(test_theta,test_phi,ell,m) ) );
			
			fA_numerical = vector_sum( fA_numerical, scalar_times_vector( coeffs_fA_r_numerical[ell][ell+m], VSH_Y  (test_theta,test_phi,ell,m) ) );
			fA_numerical = vector_sum( fA_numerical, scalar_times_vector( coeffs_fA_1_numerical[ell][ell+m], VSH_Psi(test_theta,test_phi,ell,m) ) );
			fA_numerical = vector_sum( fA_numerical, scalar_times_vector( coeffs_fA_2_numerical[ell][ell+m], VSH_Phi(test_theta,test_phi,ell,m) ) );
			
			fA_series = vector_sum( fA_series, scalar_times_vector( coeffs_fA_r_series[ell][ell+m], VSH_Y  (test_theta,test_phi,ell,m) ) );
			fA_series = vector_sum( fA_series, scalar_times_vector( coeffs_fA_1_series[ell][ell+m], VSH_Psi(test_theta,test_phi,ell,m) ) );
			fA_series = vector_sum( fA_series, scalar_times_vector( coeffs_fA_2_series[ell][ell+m], VSH_Phi(test_theta,test_phi,ell,m) ) );
			
		}
	}
	
	std::cout << std::endl;
	cout_vector_setw( A_exact , w, "A_exact " );
	cout_vector_setw( A_series, w, "A_series" );
	std::cout << std::endl;
	cout_vector_setw( fA_exact     , w, "fA_exact    " );
	cout_vector_setw( fA_numerical , w, "fA_numerical" );
	cout_vector_setw( fA_series    , w, "fA_series   " );
	
	
	
	
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