// cd OneDrive\PhD\Codes\20221108 Vector spherical harmonics
// g++ Check_VSH_Phi_01.cpp -o Check_VSH_Phi_01

/*
Check that the automatically calculated VSHs Phi agree with the hardcoded versions.
Initial version requires the scalar Y_ell^m to have been pregenerated. Might lose this requirement in future versions if it provides a speedup.
*/

#include <iostream>
#include <iomanip>										// std::setw
#include "../Generate_Spherical_Harmonics.h"
#include "Headers_VSH/VSH_Phi_01.h"
#include "../Hardcoded_Vector_Spherical_Harmonics.h"
#include "../Vector_Operations.h"						// print_vector()

//----- Global variables and function definitions -----
void test_VSH_Phi( std::vector< std::complex<double> > (*VSH_hard_function)(double,double), double theta, double phi, int ell, int m, bool new_line=false );
void print_stdev_VSH_Phi( std::vector< std::complex<double> > (*VSH_hard_function)(double,double), std::vector<double> theta_list, std::vector<double> phi_list, int ell, int m, bool new_line=false );

int main(){
	
	//----- Define variables -----
	int ell_max = 3;
	
	int N_theta = 100;
	int N_phi   = 100;
	
	int i_test = 0;
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
	
	test_VSH_Phi( VSH_Phi_00 , theta_list[i_test], phi_list[j_test], 0,  0, true );
	
	test_VSH_Phi( VSH_Phi_11 , theta_list[i_test], phi_list[j_test], 1,  1, true );
	test_VSH_Phi( VSH_Phi_10 , theta_list[i_test], phi_list[j_test], 1,  0       );
	test_VSH_Phi( VSH_Phi_1m1, theta_list[i_test], phi_list[j_test], 1, -1       );
	
	test_VSH_Phi( VSH_Phi_22 , theta_list[i_test], phi_list[j_test], 2,  2, true );
	test_VSH_Phi( VSH_Phi_21 , theta_list[i_test], phi_list[j_test], 2,  1       );
	test_VSH_Phi( VSH_Phi_20 , theta_list[i_test], phi_list[j_test], 2,  0       );
	test_VSH_Phi( VSH_Phi_2m1, theta_list[i_test], phi_list[j_test], 2, -1       );
	test_VSH_Phi( VSH_Phi_2m2, theta_list[i_test], phi_list[j_test], 2, -2       );
	
	test_VSH_Phi( VSH_Phi_33 , theta_list[i_test], phi_list[j_test], 3,  3, true );
	test_VSH_Phi( VSH_Phi_32 , theta_list[i_test], phi_list[j_test], 3,  2       );
	test_VSH_Phi( VSH_Phi_31 , theta_list[i_test], phi_list[j_test], 3,  1       );
	test_VSH_Phi( VSH_Phi_30 , theta_list[i_test], phi_list[j_test], 3,  0       );
	test_VSH_Phi( VSH_Phi_3m1, theta_list[i_test], phi_list[j_test], 3, -1       );
	test_VSH_Phi( VSH_Phi_3m2, theta_list[i_test], phi_list[j_test], 3, -2       );
	test_VSH_Phi( VSH_Phi_3m3, theta_list[i_test], phi_list[j_test], 3, -3       );
	
	
	//----- Calculate standard deviation for each set of indices -----
	std::cout << "\nStandard deviation for each (ell,m)" << std::endl;
	
	print_stdev_VSH_Phi( VSH_Phi_00 , theta_list, phi_list, 0,  0, true );
	
	print_stdev_VSH_Phi( VSH_Phi_11 , theta_list, phi_list, 1,  1, true );
	print_stdev_VSH_Phi( VSH_Phi_10 , theta_list, phi_list, 1,  0       );
	print_stdev_VSH_Phi( VSH_Phi_1m1, theta_list, phi_list, 1, -1       );
	
	print_stdev_VSH_Phi( VSH_Phi_22 , theta_list, phi_list, 2,  2, true );
	print_stdev_VSH_Phi( VSH_Phi_21 , theta_list, phi_list, 2,  1       );
	print_stdev_VSH_Phi( VSH_Phi_20 , theta_list, phi_list, 2,  0       );
	print_stdev_VSH_Phi( VSH_Phi_2m1, theta_list, phi_list, 2, -1       );
	print_stdev_VSH_Phi( VSH_Phi_2m2, theta_list, phi_list, 2, -2       );
	
	print_stdev_VSH_Phi( VSH_Phi_33 , theta_list, phi_list, 3,  3, true );
	print_stdev_VSH_Phi( VSH_Phi_32 , theta_list, phi_list, 3,  2       );
	print_stdev_VSH_Phi( VSH_Phi_31 , theta_list, phi_list, 3,  1       );
	print_stdev_VSH_Phi( VSH_Phi_30 , theta_list, phi_list, 3,  0       );
	print_stdev_VSH_Phi( VSH_Phi_3m1, theta_list, phi_list, 3, -1       );
	print_stdev_VSH_Phi( VSH_Phi_3m2, theta_list, phi_list, 3, -2       );
	print_stdev_VSH_Phi( VSH_Phi_3m3, theta_list, phi_list, 3, -3       );
	
	return 0;
}




void test_VSH_Phi( std::vector< std::complex<double> > (*VSH_hard_function)(double,double), double theta, double phi, int ell, int m, bool new_line ){
	// Use setw() for fixed widths. We expect the r components to be zero, so no need to fix a width on this components.
	if( new_line ){
		std::cout << std::endl;
	}
	std::vector< std::complex<double> > VSH_calc = VSH_Phi( theta, phi, ell, m );
	std::vector< std::complex<double> > VSH_hard = VSH_hard_function( theta, phi );
	std::cout << VSH_calc[0] << "\t";
	std::cout << std::left << std::setw(16);
	std::cout << VSH_calc[1] <<"\t";
	std::cout << std::left << std::setw(16);
	std::cout << VSH_calc[2] <<"\t";
	std::cout << " / " << VSH_hard[0] <<"\t";
	std::cout << std::left << std::setw(16);
	std::cout << VSH_hard[1] <<"\t";
	std::cout << std::left << std::setw(16);
	std::cout << VSH_hard[2] << std::endl;
}




void print_stdev_VSH_Phi( std::vector< std::complex<double> > (*VSH_hard_function)(double,double), std::vector<double> theta_list, std::vector<double> phi_list, int ell, int m, bool new_line ){
	double N_theta = theta_list.size();
	double N_phi   = phi_list.size();
	std::complex<double> ret = 0;
	for( int i=0; i<N_theta; i++ ){
		for( int j=0; j<N_phi; j++ ){
			std::vector< std::complex<double> > VSH_calc = VSH_Phi( theta_list[i], phi_list[j], ell, m );
			std::vector< std::complex<double> > VSH_hard = VSH_hard_function( theta_list[i], phi_list[j] );
			for( int c=0; c<3; c++ ){
				ret += pow( VSH_calc[c] - VSH_hard[c], 2 );
			}
		}
	}
	
	ret = sqrt( ret / ( (double) N_theta * N_phi ) );
	
	if( new_line ){
		std::cout << std::endl;
	}
	std::cout << ell <<","<< std::right << std::setw(2) << m <<"\t"<< ret << std::endl;
}