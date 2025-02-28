// cd OneDrive\PhD\Codes\20221108 Vector spherical harmonics\Dot products
// g++ Check_dot_product_Psi_Phi_and_Phi_Psi.cpp -o Check_dot_product_Psi_Phi_and_Phi_Psi

/*
Psi_ell^m cdot Phi_ell'^m' = - Psi_ell'^m' cdot Phi_ell^m = i/sin(theta) ( m Y_ell^m d/dtheta(Y_ell'^m') - m' Y_ell'^m' d/dtheta(Y_ell^m)
*/

#include <iostream>
#include "../../Generate_Spherical_Harmonics.h"
#include "../Headers_VSH/VSH_Psi_02.h"
#include "../Headers_VSH/VSH_Phi_01.h"
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
	
	
	//----- Test for single values ----- NOT CHANGED FROM COPIED CODE 20221219
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
	
	for( int L1=0; L1<=ell_max; L1++ ){
		for( int M1=-L1; M1<=L1; M1++ ){
			for( int L2=0; L2<=ell_max; L2++ ){
				for( int M2=-L2; M2<=L2; M2++ ){
					
					std::complex<double> stdev_1 = 0;
					std::complex<double> stdev_2 = 0;
					
					for( int i=0; i<N_theta; i++ ){
						for( int j=0; j<N_phi; j++ ){
							
							double t = theta_list[i];
							double p = phi_list  [j];
					
							std::complex<double> dot_product_exact_1 = dot_product( VSH_Psi(t,p,L1,M1), VSH_Phi(t,p,L2,M2) );
							std::complex<double> dot_product_exact_2 = dot_product( VSH_Psi(t,p,L2,M2), VSH_Phi(t,p,L1,M1) );
							std::complex<double> dot_product_guess_1 = std::complex<double>{0,(double)M1} * ylm_over_sintheta(t,p,L1,M1) * ylm_dtheta(t,p,L2,M2) - std::complex<double>{0,(double)M2} * ylm_over_sintheta(t,p,L2,M2) * ylm_dtheta(t,p,L1,M1);
							std::complex<double> dot_product_guess_2 = - dot_product_guess_1;
							
							stdev_1 += pow( dot_product_exact_1 - dot_product_guess_1, 2 );
							stdev_2 += pow( dot_product_exact_2 - dot_product_guess_2, 2 );
						}
					}
					stdev_1 = sqrt( stdev_1 / ( (double) N_theta * N_phi ) );
					stdev_2 = sqrt( stdev_2 / ( (double) N_theta * N_phi ) );
					std::cout << L1 <<"\t"<< M1 <<"\t"<< L2 <<"\t"<< M2 <<"\t"<< stdev_1 <<"\t"<< stdev_2 << std::endl;
					
				}
			}
		}
	}
	
	std::cout << "ell\tm\tell'\tm'\tstdev" << std::endl;
	
	
	return 0;
}