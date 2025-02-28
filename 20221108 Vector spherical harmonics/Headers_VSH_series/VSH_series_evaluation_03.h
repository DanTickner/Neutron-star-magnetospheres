// VSH_series_evaluation_03.h
// Evaluate the VSH series of a vector A with no r-dependence at a point (theta,phi), given its series coefficients.
// Initially, the coefficients are fed in as separate vectors.
// Initially, the vector addition is handled by external functions.
// The commented-out #include "" lines represent the modules needed for this code to work, but were already included in VSH_series_calculation_xx.h.

//#include "../../Vector_Operations.h"
//#include "../Headers_VSH/VSH_Y_01.h"
//#include "../Headers_VSH/VSH_Psi_02.h"
//#include "../Headers_VSH/VSH_Phi_01.h"

std::vector< std::complex<double> > evaluate_vsh_series( double t, double p, std::vector< std::vector< std::complex<double> > > coeffs_VSH_series_r, std::vector< std::vector< std::complex<double> > > coeffs_VSH_series_1, std::vector< std::vector< std::complex<double> > > coeffs_VSH_series_2 ){
	// Generate the C^{r,ell}_m coefficients in the vector spherical harmonic series of a vector function A(theta,phi) with no r-dependence.
	
	int ell_max = coeffs_VSH_series_r.size() - 1;
	std::vector< std::complex<double> > ret { 0, 0, 0 };
	
	for( int ell=0; ell<=ell_max; ell++ ){
		for( int m=-ell; m<=ell; m++ ){
			//print_vector( ret );
			
			std::vector< std::complex<double> > ret_plus_r = scalar_times_vector_spher( coeffs_VSH_series_r[ell][ell+m], VSH_Y  ( t, p, ell, m ) );
			std::vector< std::complex<double> > ret_plus_1 = scalar_times_vector_spher( coeffs_VSH_series_1[ell][ell+m], VSH_Psi( t, p, ell, m ) );
			std::vector< std::complex<double> > ret_plus_2 = scalar_times_vector_spher( coeffs_VSH_series_2[ell][ell+m], VSH_Phi( t, p, ell, m ) );
			
			ret = vector_sum_spher( ret, ret_plus_r );
			ret = vector_sum_spher( ret, ret_plus_1 );
			ret = vector_sum_spher( ret, ret_plus_2 );
		}
	}
	
	return ret;
}