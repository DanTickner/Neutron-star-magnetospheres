// cd OneDrive\PhD\Codes\Test mathematical expressions
// g++ Y_ellminus2_ellminus2.cpp -o Y_ellminus2_ellminus2

// Y_{ell-2}^{ell-2}(theta,phi) = e^{-2phi} 2 sqrt(ell(ell-1)/((2ell+1)(2ell-1))) Y_ell^ell(theta,phi) / sin^2(theta,phi)

#include <iostream>
#include <math.h>
#include "../Generate_Basis_Functions.h"


//----- Global variables and function declarations -----
std::complex<double> LHS( double theta, double phi, int ell );	// Index m not needed.
std::complex<double> RHS( double theta, double phi, int ell );


int main(){
	
	//----- Define variables -----
	int    ell_min = 2;
	int    ell_max = 10;
	double theta = 0.943;
	double phi   = 0.716;
	
	
	//----- Prepare spherical harmonic functions -----
	generate_spherical_harmonic_coefficients( ell_max + 1 );	// First derivative needs one ell higher.
	
	
	
	//----- Output results for only the half-integers -----
	for( int ell=ell_min; ell<=ell_max; ell++ ){
		std::cout << ell << " -> " << ell-2 << ":\t" << LHS( theta, phi, ell ) << "\t" << RHS( theta, phi, ell ) << std::endl;
	}
	
	std::cout << "Finished" << std::endl;
	
	
	return 0;
}




std::complex<double> LHS( double theta, double phi, int ell ){
	return ylm( theta, phi, ell-2, ell-2 );
}

std::complex<double> RHS( double theta, double phi, int ell ){
	//return std::complex<double>{cos(2*phi),-sin(2*phi)} * 2.0 * sqrt(ell*(ell-1)/((double)(2*ell+1)*(2*ell-1))) * ylm(theta,phi,ell,ell) / pow( sin(theta), 2 );
	return std::complex<double>{cos(2*phi),-sin(2*phi)} * 2.0 * sqrt(ell*(ell-1)/((double)(2*ell+1)*(2*ell-1))) * exp(std::complex<double>{0,ell*phi}) * sqrt(tgamma(2*ell+1+1)/(4*pi)) * pow(-1,ell) * pow( pow(2,ell) * tgamma(ell+1), -1 ) * pow(sin(theta),ell-2 );
}