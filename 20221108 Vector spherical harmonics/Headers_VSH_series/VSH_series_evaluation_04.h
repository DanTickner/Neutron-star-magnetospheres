// VSH_series_evaluation_04.h
// Evaluate the VSH series of a vector A with no r-dependence at a point (theta,phi), given its series coefficients.
// Try the naive method of just adding the components together.

//#include "../../Vector_Operations.h"
#include "../Headers_VSH/VSH_Y_01.h"
#include "../Headers_VSH/VSH_Psi_02.h"
#include "../Headers_VSH/VSH_Phi_01.h"

std::vector< std::complex<double> > evaluate_vsh_series( double t, double p, std::vector< std::vector< std::complex<double> > > coeffs_VSH_series_r, std::vector< std::vector< std::complex<double> > > coeffs_VSH_series_1, std::vector< std::vector< std::complex<double> > > coeffs_VSH_series_2 ){
	// Generate the C^{r,ell}_m coefficients in the vector spherical harmonic series of a vector function A(theta,phi) with no r-dependence.
	
	int ell_max = coeffs_VSH_series_r.size() - 1;
	std::vector< std::complex<double> > ret { 0, 0, 0 };
	
	for( int ell=0; ell<=ell_max; ell++ ){
		for( int m=-ell; m<=ell; m++ ){
			std::vector< std::complex<double> > Y   = VSH_Y  ( t, p, ell, m );
			std::vector< std::complex<double> > Psi = VSH_Psi( t, p, ell, m );
			std::vector< std::complex<double> > Phi = VSH_Phi( t, p, ell, m );
			for( int i=0; i<3; i++ ){
				ret[i] += coeffs_VSH_series_r[ell][ell+m] * Y[i] + coeffs_VSH_series_1[ell][ell+m] * Psi[i] + coeffs_VSH_series_2[ell][ell+m] * Phi[i];
			}
		}
	}
	
	return ret;
}