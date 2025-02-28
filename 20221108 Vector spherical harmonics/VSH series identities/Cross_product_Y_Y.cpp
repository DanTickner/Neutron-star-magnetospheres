// cd OneDrive\PhD\Codes\20221108 Vector spherical harmonics\VSH series identities
// g++ Cross_product_Y_Y.cpp -o Cross_product_Y_Y

/*
Cross products of vector spherical harmonics.
Y_{ell_1}^{m_1} x Y_{ell_2}^{m_2} = 0.
*/

#include <iostream>
#include <iomanip>
#include "../../Generate_Spherical_Harmonics.h"
#include "../../Vector_Operations.h"
#include "../Headers_VSH/VSH_Y_01.h"

//----- Function definitions -----
std::vector< std::complex<double> > exact( double theta, double phi, int ell_1, int m_1, int ell_2, int m_2 );
std::vector< std::complex<double> > test ( double theta, double phi, int ell_1, int m_1, int ell_2, int m_2 );

int main(){
	
	//----- Define variables -----
	int ell_max = 3;
	int N_t     = 30;
	int N_p     = 30;
	
	int i_test = 18;
	int j_test = 20;
	
	std::vector<double> t_list ( N_t );
	std::vector<double> p_list ( N_p );
	
	int sw = 18;	// std::setw( sw ) for consistent spacing.
	
	
	//----- Fail if the test indices are greater than the maximum number of points -----
	if( i_test > N_t-1 ){
		std::cout << "Test theta index (" << i_test << ") exceeds number of datapoints (" << N_t << "). Stopping the code to prevent undefined behaviour." << std::endl;
		return 1;
	}
	if( j_test > N_p-1 ){
		std::cout << "Test phi index (" << j_test << ") exceeds number of datapoints (" << N_p << "). Stopping the code to prevent undefined behaviour." << std::endl;
		return 1;
	}
	
	
	//----- Build coordinate lists -----
	for( int i=0; i<N_t; i++ ){
		t_list[i] = i * pi / ( (double) N_t - 1 );
	}
	for( int j=0; j<N_p; j++ ){
		p_list[j] = j * 2*pi / ( (double) N_p - 1 );
	}
	
	
	//----- Generate spherical harmonic coefficients -----
	generate_coeffs_spherical_harmonic( ell_max );
	
	
	//----- Calculate cross products at the test points -----
	for( int ell_1=0; ell_1<=ell_max; ell_1++ ){
		for( int m_1=-ell_1; m_1<=ell_1; m_1++ ){
			for( int ell_2=0; ell_2<=ell_max; ell_2++ ){
				for( int m_2=-ell_2; m_2<=ell_2; m_2++ ){
					
					std::vector< std::complex<double> > LHS = exact( t_list[i_test], p_list[j_test], ell_1, m_1, ell_2, m_2 );
					std::vector< std::complex<double> > RHS = test ( t_list[i_test], p_list[j_test], ell_1, m_1, ell_2, m_2 );
					
					std::cout << "(" << ell_1 <<","<< std::right<<std::setw(2) << m_1 << "), (" << ell_2 <<","<< std::right<<std::setw(2) << m_2 << "):\t\t"
					<< std::left<<std::setw(sw) << LHS[0] << std::left<<std::setw(sw) << LHS[1] << std::left<<std::setw(sw) << LHS[2] << "\t//\t"
					<< std::left<<std::setw(sw) << RHS[0] << std::left<<std::setw(sw) << RHS[1] << std::left<<std::setw(sw) << RHS[2] << std::endl;
					
				}
			}
		}
	}
	std::cout << std::endl;
	
	
	//----- Calculate stdev of expression over all points -----
	std::vector< std::complex<double> > stdev_total { 0, 0, 0 };
	
	for( int ell_1=0; ell_1<=ell_max; ell_1++ ){
		for( int m_1=-ell_1; m_1<=ell_1; m_1++ ){
			for( int ell_2=0; ell_2<=ell_max; ell_2++ ){
				for( int m_2=-ell_2; m_2<=ell_2; m_2++ ){
					
					std::vector< std::complex<double> > stdev { 0, 0, 0 };
					for( int i=0; i<N_t; i++ ){
						for( int j=0; j<N_p; j++ ){
							
							std::vector< std::complex<double> > LHS = exact( t_list[i], p_list[j], ell_1, m_1, ell_2, m_2 );
							std::vector< std::complex<double> > RHS = test ( t_list[i], p_list[j], ell_1, m_1, ell_2, m_2 );
							for( int c=0; c<3; c++ ){
								stdev[c] += pow( LHS[c] - RHS[c], 2 );
							}
						}
					}
					
					for( int c=0; c<3; c++ ){
						stdev_total[c] += stdev[c];
						stdev[c] = sqrt( stdev[c] ) / ( (double) N_t * N_p );
					}
					
					std::cout << "(" << ell_1 <<","<< std::right<<std::setw(2) << m_1 << "), (" << ell_2 <<","<< std::right<<std::setw(2) << m_2 << "):\t\t"
					<< std::left<<std::setw(sw) << stdev[0] << std::left<<std::setw(sw) << stdev[1] << std::left<<std::setw(sw) << stdev[2] << std::endl;
				}
			}
		}
	}
	
	for( int c=0; c<3; c++ ){
		stdev_total[c] = sqrt( stdev_total[c] ) / ( (double) N_t * N_p );
	}
	std::cout << "\nstdev total:\t" << std::left<<std::setw(sw) << stdev_total[0] << std::left<<std::setw(sw) << stdev_total[1] << std::left<<std::setw(sw) << stdev_total[2] << std::endl;
	
	
	return 0;
}




std::vector< std::complex<double> > exact( double theta, double phi, int ell_1, int m_1, int ell_2, int m_2 ){
	return cross_product( VSH_Y( theta, phi, ell_1, m_1 ), VSH_Y( theta, phi, ell_2, m_2 ) );
}


std::vector< std::complex<double> > test ( double theta, double phi, int ell_1, int m_1, int ell_2, int m_2 ){
	return std::vector< std::complex<double> > { 0, 0, 0 };
}