// VSH_series_calculation_18.h
// Calculate the VSH series coefficients for a vector A with no r-dependence.
// Get rid of the separate steps for positive and negative m. Still need separate step for zero m, because otherwise get 1/0 error.

#include "../../Numerical_Integration.h"
	
std::vector< std::vector< std::vector< std::complex<double> > > > generate_vsh_series(
                                                                                       std::vector< std::vector< std::vector< std::complex<double> > > >& A,
																					   std::vector<double>& t_list,
																					   std::vector<double>& p_list,
																					   std::vector< std::vector< std::vector< std::vector< std::complex<double> > > > >& ylm_vals,
																					   std::vector< std::vector< std::vector< std::vector< std::complex<double> > > > >& ylm_dtheta_vals,
																					   std::vector< std::vector< std::vector< std::vector< std::complex<double> > > > >& ylm_over_sintheta_vals ){
	// Generate the C^{r,ell}_m coefficients in the vector spherical harmonic series of a vector function A(theta,phi) with no r-dependence.
	
	
	int ell_max = ylm_vals.size() - 1;
	int N_t = t_list.size();
	int N_p = p_list.size();
	
	std::vector< std::vector< std::vector< std::complex<double> > > > coeffs ( 3, ( std::vector< std::vector< std::complex<double> > > ( ell_max+1, ( std::vector< std::complex<double> > ( 2*ell_max+1 ) ) ) ) );
	
	//----- ell=0 (only an r-coefficient) -----
	coeffs[0][0][0] = integral_spherical_surface( A[0], t_list, p_list ) / sqrt(4.0*pi);
	coeffs[1][0][0] = 0;
	coeffs[2][0][0] = 0;
	
	//----- Loop through ell -----
	for( int ell=1; ell<=ell_max; ell++ ){
		
		//----- m=0 coefficients -----
		std::vector< std::vector< std::complex<double> > > integrand_r_zero( N_t, std::vector< std::complex<double> > ( N_p ) );
		std::vector< std::vector< std::complex<double> > > integrand_1_zero( N_t, std::vector< std::complex<double> > ( N_p ) );
		std::vector< std::vector< std::complex<double> > > integrand_2_zero( N_t, std::vector< std::complex<double> > ( N_p ) );
		for( int i=0; i<N_t; i++ ){
			for( int j=0; j<N_p; j++ ){
				integrand_r_zero[i][j] = A[0][i][j] * ylm_vals       [ell][ell][i][j];
				integrand_1_zero[i][j] = A[1][i][j] * ylm_dtheta_vals[ell][ell][i][j];
				integrand_2_zero[i][j] = A[2][i][j] * ylm_dtheta_vals[ell][ell][i][j];
			}
		}
		coeffs[0][ell][ell] = integral_spherical_surface( integrand_r_zero, t_list, p_list );
		coeffs[1][ell][ell] = integral_spherical_surface( integrand_1_zero, t_list, p_list ) / ( (double) ell*(ell+1) );
		coeffs[2][ell][ell] = integral_spherical_surface( integrand_2_zero, t_list, p_list ) / ( (double) ell*(ell+1) );
		
		//----- m!=0 coefficients -----
		for( int m=-ell; m<=ell; m++ ){
			if( m==0 ){ continue; }
			std::vector< std::vector< std::complex<double> > > integrand_r( N_t, std::vector< std::complex<double> > ( N_p ) );
			std::vector< std::vector< std::complex<double> > > integrand_1( N_t, std::vector< std::complex<double> > ( N_p ) );
			std::vector< std::vector< std::complex<double> > > integrand_2( N_t, std::vector< std::complex<double> > ( N_p ) );
			for( int i=0; i<N_t; i++ ){
				for( int j=0; j<N_p; j++ ){
					integrand_r[i][j] = A[0][i][j] * std::conj( ylm_vals[ell][ell+m][i][j] );//ylm(t_list[i],p_list[j],ell,m));//
					integrand_1[i][j] = A[1][i][j] * std::conj( ylm_dtheta_vals[ell][ell+m][i][j] ) + A[2][i][j] * std::complex<double>{0,(double)-m} * std::conj( ylm_over_sintheta_vals[ell][ell+m][i][j] );
					integrand_2[i][j] = A[2][i][j] * std::conj( ylm_dtheta_vals[ell][ell+m][i][j] ) + A[1][i][j] * std::complex<double>{0,(double) m} * std::conj( ylm_over_sintheta_vals[ell][ell+m][i][j] );
				}
			}
			
			coeffs[0][ell][ell+m] = integral_spherical_surface( integrand_r, t_list, p_list );
			coeffs[1][ell][ell+m] = integral_spherical_surface( integrand_1, t_list, p_list ) / ( (double) ell*(ell+1) );
			coeffs[2][ell][ell+m] = integral_spherical_surface( integrand_2, t_list, p_list ) / ( (double) ell*(ell+1) );
		}
		
	}
	
	return coeffs;
}