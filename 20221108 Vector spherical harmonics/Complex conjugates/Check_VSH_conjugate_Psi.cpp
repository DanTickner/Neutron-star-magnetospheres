// cd OneDrive\PhD\Codes\20221108 Vector spherical harmonics\Complex conjugates
// g++ Check_VSH_conjugate_Psi.cpp -o Check_VSH_conjugate_Psi

/*
Check that ( Psi_ell^m )^* = (=1)^m Psi_ell^m.
*/

#include <iostream>
#include "../../Generate_Spherical_Harmonics.h"
#include "../Headers_VSH/VSH_Psi_02.h"
#include "../../Vector_Operations.h"						// print_vector()

int main(){
	
	//----- Define variables -----
	int    ell   = 3;
	int    m     = 2;
	double theta = 0.3 * pi;
	double phi   = 1.1 * pi;
	
	
	//----- Generate spherical harmonic coefficients -----
	generate_coeffs_spherical_harmonic( ell );
	
	
	//----- Calculate correct and guessed values of the VSH -----
	std::vector< std::complex<double> > Psi_ell_m = VSH_Psi( theta, phi, ell, m );
	std::vector< std::complex<double> > Psi_ell_m_star_correct( 3 );
	for( int i=0; i<3; i++ ){
		Psi_ell_m_star_correct[i] = std::conj( Psi_ell_m[i] );
	}
	
	std::vector< std::complex<double> > Psi_ell_minusm = VSH_Psi( theta, phi, ell, -m );
	std::vector< std::complex<double> > Psi_ell_m_star_guess( 3 );
	for( int i=0; i<3; i++ ){
		Psi_ell_m_star_guess  [i] = pow( -1, m ) * Psi_ell_minusm[i];
	}
	
	//----- Output results -----
	print_vector( Psi_ell_m_star_correct );
	print_vector( Psi_ell_m_star_guess   );
	
	return 0;
}