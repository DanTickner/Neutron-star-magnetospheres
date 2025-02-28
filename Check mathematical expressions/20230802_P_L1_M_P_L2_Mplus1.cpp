/*
cd OneDrive\PhD\Codes\Check mathematical expressions
g++ 20230802_P_L1_M_P_L2_Mplus1.cpp -o 20230802_P_L1_M_P_L2_Mplus1 --std=c++11
20230802_P_L1_M_P_L2_Mplus1


P_{L1}^M P_{L2}^{M+1} = 

OneNote Force-Free Documentation Proofs -> Spherical harmonics
*/

#include <iostream>
#include <iomanip>
#include <math.h>
#include <vector>
#include <chrono>

#include "../Generate_Spherical_Harmonics.h"
#include "../Generate_Associated_Legendre_Functions.h"


int main(){
	
	//----- Define variables -----
	
	int L1 = 3;
	int L2 = 5;
	int M = 2;
	double x = 0.345;
	int w = 40;			// Fixed width.
	
	
	generate_coeffs_spherical_harmonic( L1+L2 );
	generate_coeffs_associated_legendre( L1+L2 );
	
	double P1 = associated_legendre_function( x, L1, M   );
	double P2 = associated_legendre_function( x, L2, M+1 );
	
	std::cout << "P1 = P_" << L1 <<"^"<< M   <<" = "<< P1 << std::endl;
	std::cout << "P2 = P_" << L2 <<"^"<< M+1 <<" = "<< P2 << std::endl;
	
	double product_LHS = P1 * P2;
	
	double product_RHS = 0;
	
	
	
	/*
	//----- V1 -----
	int M1 = M;
	int M2 = M+1;
	for( int ell=abs(L1-L2); ell<=L1+L2; ell++ ){
		double G = ( pow( -1, -M1+M2 ) * ( 2.0*ell+1.0 ) * wigner_3j( L1, L2, ell, 0, 0, 0 ) * wigner_3j( L1, L2, ell, -M1, M2, M1-M2 ) ).real();
		double sqrt_factor = sqrt( tgamma( ell-M2+M1 +1 ) / tgamma( ell+M2-M1 +1 ) );
		double ALF_ell = associated_legendre_function( x, ell, -M1+M2 );
		product_RHS += G * sqrt_factor * ALF_ell;
	}
	product_RHS *= pow( -1, M1 ) * sqrt( tgamma( L1+M1 +1 ) * tgamma( L2+M2 +1 ) / ( tgamma( L1-M1 + 1 ) * tgamma( L2-M2 + 1 ) ) );
	*/
	
	for( int ell=abs(L1-L2); ell<=L1+L2; ell++ ){
		double G = - ( 2.0*ell+1.0 ) * ( wigner_3j( L1, L2, ell, 0, 0, 0 ) * wigner_3j( L1, L2, ell, -M, M+1, -1 ) ).real();
		double sqrt_factor = 1.0 / sqrt( ell*(ell+1.0) );
		double ALF_ell = associated_legendre_function( x, ell, 1 );
		product_RHS += G * sqrt_factor * ALF_ell;
	}
	product_RHS *= pow( -1, M ) * sqrt( (L2-M) * ( L2+M+1.0 ) * tgamma( L1+M +1 ) * tgamma( L2+M +1 ) / ( tgamma( L1-M + 1 ) * tgamma( L2-M + 1 ) ) );
	
	
	
	
	
	
	std::cout << "P1*P2\t" << product_LHS <<"\t"<< product_RHS << std::endl;
	
	std::cout << "\ndone" << std::endl;
	
	
	return 0;
	
}