// VSH_series_calculation_101.h
// Calculate the VSH series coefficients for a vector A with no r-dependence.
// Start again.

#include "../../Numerical_Integration.h"
#include "../../Vector_Operations.h"
#include "../Headers_VSH/VSH_Y_01.h"
#include "../Headers_VSH/VSH_Psi_02.h"
#include "../Headers_VSH/VSH_Phi_01.h"
	
std::vector< std::vector< std::vector< std::complex<double> > > > generate_vsh_series( std::vector< std::vector< std::vector< std::complex<double> > > > A,
																					   std::vector<double> t_list,
																					   std::vector<double> p_list,
																					   int ell_max ){
	// Generate the C^{r,ell}_m coefficients in the vector spherical harmonic series of a vector function A(theta,phi) with no r-dependence.
	
	int N_t = t_list.size();
	int N_p = p_list.size();
	
	std::vector< std::vector< std::vector< std::complex<double> > > > coeffs ( 3, ( std::vector< std::vector< std::complex<double> > > ( ell_max+1, ( std::vector< std::complex<double> > ( 2*ell_max+1 ) ) ) ) );
	
	
	for( int ell=0; ell<=ell_max; ell++ ){
		for( int m=-ell; m<=ell; m++ ){
			
			std::vector< std::vector< std::complex<double> > > integrand_r( N_t, std::vector< std::complex<double> > ( N_p ) );
			std::vector< std::vector< std::complex<double> > > integrand_1( N_t, std::vector< std::complex<double> > ( N_p ) );
			std::vector< std::vector< std::complex<double> > > integrand_2( N_t, std::vector< std::complex<double> > ( N_p ) );
			
			for( int i=0; i<N_t; i++ ){
				for( int j=0; j<N_p; j++ ){
					std::vector< std::complex<double> > A_at_this_point { A[0][i][j], A[1][i][j], A[2][i][j] };
					integrand_r[i][j] = dot_product_cart( A_at_this_point, vector_conj( VSH_Y  (t_list[i],p_list[j],ell,m) ) );
					integrand_1[i][j] = dot_product_cart( A_at_this_point, vector_conj( VSH_Psi(t_list[i],p_list[j],ell,m) ) );
					integrand_2[i][j] = dot_product_cart( A_at_this_point, vector_conj( VSH_Phi(t_list[i],p_list[j],ell,m) ) );
				}
			}
			
			coeffs[0][ell][ell+m] = integral_spherical_surface( integrand_r, t_list, p_list );
			coeffs[1][ell][ell+m] = integral_spherical_surface( integrand_1, t_list, p_list );
			coeffs[2][ell][ell+m] = integral_spherical_surface( integrand_2, t_list, p_list );
			
			if( ell != 0 ){
				coeffs[1][ell][ell+m] /= (double) ell*(ell+1);
				coeffs[2][ell][ell+m] /= (double) ell*(ell+1);
			}
		}
	}
	
	return coeffs;
}