// cd OneDrive\PhD\Codes\20221108 Vector spherical harmonics\Dot products
// g++ Check_dot_product_Psi_Psi.cpp -o Check_dot_product_Psi_Psi -std=c++11

/*
Psi_ell^m cdot Psi_ell'^m' = d/dtheta(Y_ell^m) d/dtheta(Y_ell'^m') - m m'/sin^2(theta) * Y_ell^m Y_ell'^m'
*/

#include <iostream>
#include "../../Generate_Spherical_Harmonics.h"
#include "../Headers_VSH/VSH_Psi_01.h"
#include "../../Vector_Operations.h"

int main(){
	
	//----- Define variables -----
	int ell_max = 3;
	
	int N_theta = 100;
	int N_phi   = 100;
	
	std::vector<double> theta_list ( N_theta );
	std::vector<double> phi_list   ( N_phi   );
	
	
	//----- Build coordinate lists -----
	double theta_min = 0.1 * pi;	// 1/sin^2(theta) term appears
	double theta_max = 0.9 * pi;
	for( int i=0; i<N_theta; i++ ){
		theta_list[i] = theta_min + i * ( theta_max - theta_min ) / ( (double) N_theta - 1 );
	}
	for( int j=0; j<N_phi; j++ ){
		phi_list[j] = j * 2*pi / ( (double) N_phi - 1 );
	}
	
	
	//----- Generate spherical harmonic coefficients -----
	generate_coeffs_spherical_harmonic( ell_max );
	
	
	//----- Test for single values -----
	/*
	double theta_test = theta_list[2];
	double phi_test   = phi_list  [2];
	
	for( int ell=0; ell<=ell_max; ell++ ){
		for( int m=-ell; m<=ell; m++ ){
			for( int ell_prime=0; ell_prime<=ell_max; ell_prime++ ){
				for( int m_prime=-ell_prime; m_prime<=ell_prime; m_prime++ ){
					
					std::complex<double> dot_product_exact = dot_product_cart( VSH_Psi( theta_test, phi_test, ell, m ), VSH_Psi( theta_test, phi_test, ell_prime, m_prime ) );
					std::complex<double> dot_product_guess = ylm_dtheta( theta_test, phi_test, ell, m ) * ylm_dtheta( theta_test, phi_test, ell_prime, m_prime ) - ( (double) m*m_prime ) * pow(sin(theta_test),-2) * ylm( theta_test, phi_test, ell, m ) * ylm( theta_test, phi_test, ell_prime, m_prime );
					
					std::cout << dot_product_exact <<"\t"<< dot_product_guess <<"\t"<< dot_product_exact - dot_product_guess << std::endl;
				}
			}
		}
	}
	*/
	
	
	//----- Calculate stdev over entire range of theta,phi -----
	
	std::cout << "ell\tm\tell'\tm'\tstdev" << std::endl;
	
	for( int ell=0; ell<=ell_max; ell++ ){
		for( int m=-ell; m<=ell; m++ ){
			for( int ell_prime=0; ell_prime<=ell_max; ell_prime++ ){
				for( int m_prime=-ell_prime; m_prime<=ell_prime; m_prime++ ){
					
					std::complex<double> stdev = 0;
					
					for( int i=0; i<N_theta; i++ ){
						for( int j=0; j<N_phi; j++ ){
					
							std::complex<double> dot_product_exact = dot_product_cart( VSH_Psi( theta_list[i], phi_list[j], ell, m ), VSH_Psi( theta_list[i], phi_list[j], ell_prime, m_prime ) );
							std::complex<double> dot_product_guess = ylm_dtheta( theta_list[i], phi_list[j], ell, m ) * ylm_dtheta( theta_list[i], phi_list[j], ell_prime, m_prime ) - ( (double) m*m_prime ) * pow(sin(theta_list[i]),-2) * ylm( theta_list[i], phi_list[j], ell, m ) * ylm( theta_list[i], phi_list[j], ell_prime, m_prime );
							
							stdev += pow( dot_product_exact - dot_product_guess, 2 );
						}
					}
					stdev = sqrt( stdev / ( (double) N_theta * N_phi ) );
					std::cout << ell <<"\t"<< m <<"\t"<< ell_prime <<"\t"<< m_prime <<"\t"<< stdev << std::endl;
					
				}
			}
		}
	}
	
	std::cout << "ell\tm\tell'\tm'\tstdev" << std::endl;
	
	
	return 0;
}