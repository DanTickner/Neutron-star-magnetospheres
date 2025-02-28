// VSH_series_evaluation_02.h
// Evaluate the VSH series of a vector A with no r-dependence at a point (theta,phi), given its series coefficients.
// The radial component is simply the sum over VSH Y.

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
			ret[0] += coeffs_VSH_series_r[ell][ell+m] * ylm( t, p, ell, m );
			std::vector< std::complex<double> > ret_plus_1 = scalar_times_vector_spher( coeffs_VSH_series_1[ell][ell+m], VSH_Psi( t, p, ell, m ) );
			std::vector< std::complex<double> > ret_plus_2 = scalar_times_vector_spher( coeffs_VSH_series_2[ell][ell+m], VSH_Phi( t, p, ell, m ) );
			
			// Development on hiatus until vector_sum_spher() has been defined for complex vectors. 20221118	
			// Trial simple addition rule for the vectors.
			//ret[1] += ret_plus_1[1] + ret_plus_2[1];
			//ret[2] += ret_plus_1[2] + ret_plus_2[2];
		}
	}
	
	return ret;
}