// cd OneDrive\PhD\Codes\20221108 Vector spherical harmonics
// g++ Check_VSH_negative_m_Y.cpp -o Check_VSH_negative_m_Y

/*
Check the following relation:
Y_ell^{-m} = (-1)^m e^{-i2mphi} Y_ell^m
where Y_ell^m are the first vector spherical harmonics.
*/

#include <iostream>
#include <iomanip>										// std::setw
#include "../Generate_Spherical_Harmonics.h"
#include "Headers_VSH/VSH_Y_01.h"
#include "../Vector_Operations.h"						// scalar_times_vector_cart()

int main(){
	
	//----- Define variables -----
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
			std::vector< std::complex<double> > VSH_test  = scalar_times_vector_cart( pow(-1,m)*exp(std::complex<double>{0,-2*m*phi_list[j_test]}), VSH_Y( theta_list[i_test], phi_list[j_test], ell, m ) );
			std::vector< std::complex<double> > VSH_exact = VSH_Y( theta_list[i_test], phi_list[j_test], ell, -m );
			std::cout << ell << "\t" << m << "\t" << std::left << std::setw(18);
			std::cout << VSH_test[0] <<"\t"<< VSH_test[1] <<"\t"<< VSH_test[2];
			std::cout << "\t/\t" << std::left << std::setw(18);
			std::cout << VSH_exact[0] <<"\t"<< VSH_exact[1] <<"\t"<< VSH_exact[2];
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
					std::vector< std::complex<double> > VSH_test  = scalar_times_vector_cart( pow(-1,m)*exp(std::complex<double>{0,-2*m*phi_list[j_test]}), VSH_Y( theta_list[i_test], phi_list[j_test], ell, m ) );
					std::vector< std::complex<double> > VSH_exact = VSH_Y( theta_list[i_test], phi_list[j_test], ell, -m );
					for( int c=0; c<3; c++ ){
						stdev += pow( VSH_test[c] - VSH_exact[c], 2 );
					}
				}
			}
			
			stdev = sqrt( stdev / ( (double) N_theta * N_phi ) );
			std::cout << ell <<","<< std::right << std::setw(2) << m <<"\t"<< stdev << std::endl;
		}
		std::cout << std::endl;
	}
	
	return 0;
}