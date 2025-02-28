// cd OneDrive\PhD\Codes\20221108 Vector spherical harmonics\Complex conjugates
// g++ Check_spherical_harmonic_conjugate.cpp -o Check_spherical_harmonic_conjugate

/*
Check that ( Y_ell^m )^* = (-1)^m Y_ell^m, where we consider the vector, not the scalar spherical harmonic.
*/

#include <iostream>
#include "../../Generate_Spherical_Harmonics.h"

int main(){
	
	//----- Define variables -----
	int    ell   = 3;
	int    m     = 2;
	double theta = 0.3 * pi;
	double phi   = 1.1 * pi;
	
	
	//----- Generate spherical harmonic coefficients -----
	generate_coeffs_spherical_harmonic( ell );
	
	
	//----- Calculate correct and guessed values of the VSH -----	
	std::complex<double> Y_ell_m_star_correct = std::conj( ylm( theta, phi, ell, m ) );
	std::complex<double> Y_ell_m_star_guess = pow( -1, m ) * ylm( theta, phi, ell, -m );
	
	//----- Output results -----
	std::cout << Y_ell_m_star_correct <<"\t"<< Y_ell_m_star_guess << std::endl;
	
	return 0;
}