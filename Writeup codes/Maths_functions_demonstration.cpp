/*
cd OneDrive\PhD\Codes\20240311 Time evolution code\Writeup codes
g++ Maths_functions_demonstration.cpp -o Maths_functions_demonstration
Maths_functions_demonstration

*/

#include <iostream>
#include <iomanip>
//#include "../Numerical_Integration.h"
#include "Maths_functions.h"

//----- Global variables and function declarations -----




//----- Main code -----
int main(){
	
	double x     = 0.321;		// Argument for Legendre polynomials and associated Legendre functions.
	double theta = 0.843*pi;	// Argument for spherical harmonics and vector spherical harmonics (polar).
	double phi   = 1.234*pi;	// Argument for spherical harmonics and vector spherical harmonics (azimuthal).
	int    ell   = 4;
	double h     = 1e-5;     // Finite difference for numerical derivatives.
	
	/*
	std::cout << "Legendre polynomials, x = " << x << std::endl;
	for( int L=0; L<=10; L++ ){
		std::cout << L <<"\t"<< P_ell( x, L ) << std::endl;
	}
	*/
	
	/*
	std::cout << "Factorials aand double factorials up to 20!" << std::endl;
	for( int n=1; n<=20; n++ ){
		std::cout << n <<"\t"<< factorial( n ) <<"\t"<< double_factorial( n ) << std::endl;
	}
	*/
	
	/*
	std::cout << "\nAssociated Legendre functions, x=" << x << ", ell=" << ell <<", varying m" << std::endl;
	for( int m=-ell; m<=ell; m++ ){
		std::cout << m <<"\t"<< P_ell_m( x, ell, m ) << std::endl;
	}
	*/
	
	/*
	ell = 7;
	std::cout << "\nSpherical harmonics, ell=" << ell << ", varying m" << std::endl;
	for( int m=-ell; m<=ell; m++ ){
		std::cout << m <<"\t"<< Y_ell_m( theta, phi, ell, m ) << std::endl;
	}
	*/
	
	/*
	ell = 7;
	std::cout << "\nSpherical harmonics d/dtheta, ell=" << ell << ", varying m" << std::endl;
	for( int m=-ell; m<=ell; m++ ){
		std::complex<double> derivative = ( Y_ell_m( theta + h, phi, ell, m ) - Y_ell_m( theta, phi, ell, m ) ) / h;
		std::cout << m <<"\t"<< Y_ell_m_dtheta( theta, phi, ell, m ) <<"\t"<< derivative << std::endl;
	}
	*/
	
	/*
	std::cout << "\nSpherical harmonics 1/sin(theta) d/dphi, ell=" << ell << ", varying m" << std::endl;
	for( int m=-ell; m<=ell; m++ ){
		std::complex<double> derivative = ( Y_ell_m( theta, phi + h, ell, m ) - Y_ell_m( theta, phi, ell, m ) ) / ( h * sin( theta ) );
		std::cout << m <<"\t"<< Y_ell_m_dphi_over_sintheta( theta, phi, ell, m ) <<"\t"<< derivative << std::endl;
	}
	*/
	
	/*
	ell = 7;
	std::cout << "\nVector spherical harmonics Y, ell=" << ell << ", varying m" << std::endl;
	for( int m=-ell; m<=ell; m++ ){
		std::cout << m <<"\t";
		print_vector( VSH_Y( theta, phi, ell, m ) );
	}
	*/
	
	/*
	ell = 3;
	std::cout << "\nVector spherical harmonics Psi, ell=" << ell << ", varying m" << std::endl;
	for( int m=-ell; m<=ell; m++ ){
		std::cout << m <<"\t";
		print_vector( VSH_Psi( theta, phi, ell, m ) );
	}
	*/
	
	ell = 3;
	std::cout << "\nVector spherical harmonics Phi, ell=" << ell << ", varying m" << std::endl;
	for( int m=-ell; m<=ell; m++ ){
		std::cout << m <<"\t";
		print_vector( VSH_Phi( theta, phi, ell, m ) );
	}
	
	
	
	
	return 0;
	
}