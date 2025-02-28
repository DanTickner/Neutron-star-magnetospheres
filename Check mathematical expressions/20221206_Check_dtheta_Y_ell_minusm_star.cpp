// cd OneDrive\PhD\Codes\Test mathematical expressions
// g++ 20221206_Check_dtheta_Y_ell_minusm_star.cpp -o 20221206_Check_dtheta_Y_ell_minusm_star

/*
Check the following relation:
d/dtheta [Y_ell^{-m}]^* = (-1)^m d/dtheta Y_ell^m.
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
		for( int m=0; m<=ell; m++ ){
			std::complex<double> SH_exact = std::conj( ylm_dtheta( theta_list[i_test], phi_list[j_test], ell, -m ) );
			std::complex<double> SH_test  = pow(-1,m) * ylm_dtheta( theta_list[i_test], phi_list[j_test], ell, m );
			
			std::cout << ell << "\t" << m << "\t";
			std::cout << std::left<<std::setw(sw) << SH_exact.real() <<"\t"<< std::left<<std::setw(sw) << SH_exact.imag();
			std::cout << "\t/\t";
			std::cout << std::left<<std::setw(sw) << SH_test .real() <<"\t"<< std::left<<std::setw(sw) << SH_test .imag();
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}
	
	
	//----- Calculate standard deviation for each set of indices -----
	std::cout << "\nStandard deviation for each (ell,m)" << std::endl;
	
	for( int ell=0; ell<=ell_max; ell++ ){
		for( int m=0; m<=ell; m++ ){
			std::complex<double> stdev = 0;
			
			for( int i=0; i<N_theta; i++ ){
				for( int j=0; j<N_phi; j++ ){
					std::complex<double> SH_exact = std::conj( ylm_dtheta( theta_list[i], phi_list[j], ell, -m ) );
					std::complex<double> SH_test  = pow(-1,m) * ylm_dtheta( theta_list[i], phi_list[j], ell, m );
					stdev += pow( SH_exact - SH_test, 2 );
				}
			}
			
			stdev = sqrt( stdev / ( (double) N_theta * N_phi ) );
			std::cout << ell <<","<< std::right << std::setw(2) << m <<"\t"<< stdev << std::endl;
		}
		std::cout << std::endl;
	}
	
	return 0;
}