// cd OneDrive\PhD\Codes\Test mathematical expressions
// g++ 20221206_Check_dphi_Y_ell_m_star.cpp -o 20221206_Check_dphi_Y_ell_m_star

/*
Check the following relations:
d/dphi [Y_ell^m]^* = - e^{-i2mphi} d/dphi Y_ell^m
                   = - e^{-i2mphi} im Y_ell^m.
*/

#include <vector>
#include <iostream>
#include <iomanip>										// std::setw
#include "../Generate_Spherical_Harmonics.h"

int main(){
	
	//----- Define variables -----
	int sw = 12;		// setw() text value.
	
	int ell_max = 3;
	
	int N_theta = 100;
	int N_phi   = 100;
	
	int i_test = 22;
	int j_test = 33;
	
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
	
	
	//----- Test for single values -----
	std::cout << "Testing for (theta,phi) = ( " << theta_list[i_test]/pi << " pi , " << phi_list[j_test]/pi << " pi )" << std::endl;
	
	for( int ell=0; ell<=ell_max; ell++ ){
		for( int m=-ell; m<=ell; m++ ){
			std::complex<double> SH_exact = std::conj( ylm_dphi( theta_list[i_test], phi_list[j_test], ell, m ) );
			std::complex<double> SH_test1 = - exp(std::complex<double>{0,-2*m*phi_list[j_test]}) * ylm_dphi( theta_list[i_test], phi_list[j_test], ell, m );
			std::complex<double> SH_test2 = - exp(std::complex<double>{0,-2*m*phi_list[j_test]}) * std::complex<double>{0,(double)m} * ylm( theta_list[i_test], phi_list[j_test], ell, m );
			
			std::cout << ell << "\t" << m << "\t";
			std::cout << std::left<<std::setw(sw) << SH_exact.real() <<"\t"<< std::left<<std::setw(sw) << SH_exact.imag();
			std::cout << "\t/\t";
			std::cout << std::left<<std::setw(sw) << SH_test1.real() <<"\t"<< std::left<<std::setw(sw) << SH_test1.imag();
			std::cout << "\t/\t";
			std::cout << std::left<<std::setw(sw) << SH_test2.real() <<"\t"<< std::left<<std::setw(sw) << SH_test2.imag();
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}
	
	
	//----- Calculate standard deviation for each set of indices -----
	std::cout << "\nStandard deviation for each (ell,m)" << std::endl;
	
	for( int ell=0; ell<=ell_max; ell++ ){
		for( int m=-ell; m<=ell; m++ ){
			std::complex<double> stdev1 = 0;
			std::complex<double> stdev2 = 0;
			
			for( int i=0; i<N_theta; i++ ){
				for( int j=0; j<N_phi; j++ ){
					std::complex<double> SH_exact = std::conj( ylm_dphi( theta_list[i], phi_list[j], ell, m ) );
					std::complex<double> SH_test1 = - exp(std::complex<double>{0,-2*m*phi_list[j]}) * ylm_dphi( theta_list[i], phi_list[j], ell, m );
					std::complex<double> SH_test2 = - exp(std::complex<double>{0,-2*m*phi_list[j]}) * std::complex<double>{0,(double)m} * ylm( theta_list[i], phi_list[j], ell, m );
					stdev1 += pow( SH_exact - SH_test1, 2 );
					stdev2 += pow( SH_exact - SH_test2, 2 );
				}
			}
			
			stdev1 = sqrt( stdev1 / ( (double) N_theta * N_phi ) );
			stdev2 = sqrt( stdev2 / ( (double) N_theta * N_phi ) );
			std::cout << ell <<","<< std::right << std::setw(2) << m <<"\t"<< stdev1 <<"\t"<< stdev2 << std::endl;
		}
		std::cout << std::endl;
	}
	
	return 0;
}