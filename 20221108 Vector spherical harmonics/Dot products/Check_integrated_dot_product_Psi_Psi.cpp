// cd OneDrive\PhD\Codes\20221108 Vector spherical harmonics\Dot products
// g++ Check_integrated_dot_product_Psi_Psi.cpp -o Check_integrated_dot_product_Psi_Psi

/*
\int Psi_ell^m cdot ( Psi_ell'^m' )^* dOmega = ell (ell+1) delta_{ell,ell'} delta^{m,m'}
*/

#include <iostream>
#include <iomanip>
#include "../../Generate_Spherical_Harmonics.h"
#include "../Headers_VSH/VSH_Psi_02.h"
#include "../../Vector_Operations.h"
#include "../../Numerical_Integration.h"

int main(){
	
	//----- Define variables -----
	int ell_max = 3;
	int N_theta = 500;
	int N_phi   = 500;
	
	std::vector<double> theta_list ( N_theta );
	std::vector<double> phi_list   ( N_phi   );
	
	
	//----- Build coordinate lists -----
	for( int i=0; i<N_theta; i++ ){
		theta_list[i] = i * pi / ( (double) N_theta - 1 );
	}
	for( int j=0; j<N_phi; j++ ){
		phi_list[j] = j * 2*pi / ( (double) N_phi - 1 );
	}
	
	
	//----- Generate spherical harmonic coefficients -----
	generate_coeffs_spherical_harmonic( ell_max );
	
	
	//----- Calculate surface integral over entire range of theta,phi -----
	int N_integrals = 0;	// Could calculate from m,ell relations but may as well increment as we go.
	std::complex<double> stdev_total = 0;
	std::cout << "ell\tm\tell'\tm'\t" << std::left << std::setw(30) << "integral_exact" << std::left << std::setw(30) <<  "integral_guess" << std::left << std::setw(30) << "difference" << std::endl;
	
	for( int ell=0; ell<=ell_max; ell++ ){
		for( int m=-ell; m<=ell; m++ ){
			for( int ell_prime=0; ell_prime<=ell_max; ell_prime++ ){
				for( int m_prime=-ell_prime; m_prime<=ell_prime; m_prime++ ){
					
					std::vector< std::vector< std::complex<double> > > integrand_values ( N_theta, std::vector< std::complex<double> > ( N_phi ) );
					for( int i=0; i<N_theta; i++ ){
						for( int j=0; j<N_phi; j++ ){
							integrand_values[i][j] = dot_product_cart( VSH_Psi( theta_list[i], phi_list[j], ell, m ), vector_conj( VSH_Psi( theta_list[i], phi_list[j], ell_prime, m_prime ) ) );
						}
					}
					
					std::complex<double> integral_exact = integral_spherical_surface( integrand_values, theta_list, phi_list );
					
					std::complex<double> integral_guess = 0;
					if( ( ell == ell_prime ) && ( m == m_prime ) ){
						integral_guess = (std::complex<double>) ( ell * (ell+1) );
					}
					
					std::cout << ell <<"\t"<< m <<"\t"<< ell_prime <<"\t"<< m_prime <<"\t";
					std::cout << std::left << std::setw(30);
					std::cout << integral_exact;
					std::cout << std::left << std::setw(30);
					std::cout << integral_guess;
					std::cout << std::left << std::setw(30);
					std::cout << integral_exact - integral_guess << std::endl;
					
					stdev_total += pow( integral_exact - integral_guess, 2 );
					N_integrals++;
					
				}
			}
		}
	}
	
	std::cout << "ell\tm\tell'\tm'\t" << std::left << std::setw(30) << "integral_exact" << std::left << std::setw(30) <<  "integral_guess" << std::left << std::setw(30) << "difference" << std::endl;
	
	std::cout << "\nTotal stdev:\t" << sqrt( stdev_total / ( (double) N_integrals ) ) << std::endl;
	
	return 0;
}