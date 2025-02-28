// VSH_series_divergence_01.h
// Evaluate divergence of a vector A with no r-dependence at a point (theta,phi), given its series coefficients.
// Initially, the coefficients are fed in as separate vectors.
// Initially, the vector addition is handled by external functions.
// The commented-out #include "" lines represent the modules needed for this code to work, but were already included in VSH_series_calculation_xx.h.

//#include "../../Vector_Operations.h"
//#include "../Headers_VSH/VSH_Y_01.h"
//#include "../Headers_VSH/VSH_Psi_02.h"
//#include "../Headers_VSH/VSH_Phi_01.h"

std::vector< std::complex<double> > divergence( double t, double p, std::vector< std::vector< std::complex<double> > > coeffs_VSH_series_1 ){
	// Generate the divergence of a vector function A(theta,phi) with no r-dependence, given its VSH series coefficients.
	
	double ret     = 0;
	int    ell_max = coeffs_VSH_series_r.size() - 1;
	
	for( int ell=0; ell<=ell_max; ell++ ){
		std::complex<double> ret_for_this_ell = 0;
		for( int m=-ell; m<=ell; m++ ){
			ret_for_this_ell += coeffs_VSH_series_1[ell][ell+m] * ylm( t, p, ell, m );
		}
		ret += (std::complex<double>) ell*(ell+1) * ret_for_this_ell;
	}
	
	return ret;
}