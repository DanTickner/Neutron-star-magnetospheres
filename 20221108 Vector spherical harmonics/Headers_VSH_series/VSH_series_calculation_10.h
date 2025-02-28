// VSH_series_calculation_10.h
// Calculate the VSH series coefficients for a vector A with no r-dependence.
// Put C_r coefficient in m=0, m!=0 form, and separate the ell=0 term, ready for future combination of all the functions.
// Header files haven't been needed for a while now.

#include "../../Numerical_Integration.h"

std::vector< std::vector< std::complex<double> > > generate_vsh_series_Cr( std::vector< std::vector< std::vector< std::complex<double> > > > A, std::vector<double> t_list, std::vector<double> p_list, int ell_max ){
	// Generate the C^{r,ell}_m coefficients in the vector spherical harmonic series of a vector function A(theta,phi) with no r-dependence.
	
	std::vector< std::vector< std::complex<double> > > coeffs;
	
	int N_t = t_list.size();
	int N_p = p_list.size();
	
	coeffs.push_back( std::vector< std::complex<double> > { integral_spherical_surface( A[0], t_list, p_list ) / sqrt(4.0*pi) } );
	
	for( int ell=1; ell<=ell_max; ell++ ){
		// The range of the index m depends on ell, so the second index of our 2D vectors has varying size and we cannot define the entire vector skeleton in one go.
		std::vector< std::complex<double> > this_ell_coeffs( 2*ell+1 );
		
		//----- m=0 coefficient -----
		std::vector< std::vector< std::complex<double> > > integrand( N_t, std::vector< std::complex<double> > ( N_p ) );
		for( int i=0; i<N_t; i++ ){
			for( int j=0; j<N_p; j++ ){
				integrand[i][j] = A[0][i][j] * ylm( t_list[i], p_list[j], ell, 0 );
			}
		}
		this_ell_coeffs[ell] = integral_spherical_surface( integrand, t_list, p_list );
		
		//----- m!=0 coefficients -----
		for( int m=1; m<=ell; m++ ){
			std::vector< std::vector< std::complex<double> > > integrand_plus ( N_t, std::vector< std::complex<double> > ( N_p ) );
			std::vector< std::vector< std::complex<double> > > integrand_minus( N_t, std::vector< std::complex<double> > ( N_p ) );
			for( int i=0; i<N_t; i++ ){
				for( int j=0; j<N_p; j++ ){
					std::complex<double> SH = ylm( t_list[i], p_list[j], ell, m );
					integrand_plus [i][j] = A[0][i][j] * std::conj( SH );
					integrand_minus[i][j] = A[0][i][j] *            SH;
				}
			}
			this_ell_coeffs[ell+m] = integral_spherical_surface( integrand_plus , t_list, p_list );
			this_ell_coeffs[ell-m] = integral_spherical_surface( integrand_minus, t_list, p_list );
		}
		
		coeffs.push_back( this_ell_coeffs );
	}
	
	return coeffs;
}




std::vector< std::vector< std::complex<double> > > generate_vsh_series_C1( std::vector< std::vector< std::vector< std::complex<double> > > > A, std::vector<double> t_list, std::vector<double> p_list, int ell_max ){
	// Generate the C^{1,ell}_m coefficients in the vector spherical harmonic series of a vector function A(theta,phi) with no r-dependence.
	
	std::vector< std::vector< std::complex<double> > > coeffs { std::vector< std::complex<double> > { 0 } };	// (ell,m)=(0,0) coefficient is zero.
	
	int N_t = t_list.size();
	int N_p = p_list.size();
	
	for( int ell=1; ell<=ell_max; ell++ ){
		// The range of the index m depends on ell, so the second index of our 2D vectors has varying size and we cannot define the entire vector skeleton in one go.
		std::vector< std::complex<double> > this_ell_coeffs( 2*ell+1 );
		
		//----- m=0 coefficient -----
		std::vector< std::vector< std::complex<double> > > integrand( N_t, std::vector< std::complex<double> > ( N_p ) );
		for( int i=0; i<N_t; i++ ){
			for( int j=0; j<N_p; j++ ){
				integrand[i][j] = A[1][i][j] * ylm_dtheta( t_list[i], p_list[j], ell, 0 );
			}
		}
		this_ell_coeffs[ell] = integral_spherical_surface( integrand, t_list, p_list ) / ( (double) ell*(ell+1) );
		
		//----- m!=0 coefficients -----
		for( int m=1; m<=ell; m++ ){
			std::vector< std::vector< std::complex<double> > > integrand_plus ( N_t, std::vector< std::complex<double> > ( N_p ) );
			std::vector< std::vector< std::complex<double> > > integrand_minus( N_t, std::vector< std::complex<double> > ( N_p ) );
			for( int i=0; i<N_t; i++ ){
				for( int j=0; j<N_p; j++ ){
					std::complex<double> SH_dtheta = ylm_dtheta( t_list[i], p_list[j], ell, m );
					std::complex<double> SH_over_sintheta = ylm_over_sintheta( t_list[i], p_list[j], ell, m );
					integrand_plus [i][j] = A[1][i][j] * std::conj( SH_dtheta ) + A[2][i][j] * std::complex<double>{0,(double)-m} * std::conj( SH_over_sintheta );
					integrand_minus[i][j] = A[1][i][j] *            SH_dtheta   + A[2][i][j] * std::complex<double>{0,(double) m} *            SH_over_sintheta;
				}
			}
			this_ell_coeffs[ell+m] =             integral_spherical_surface( integrand_plus , t_list, p_list ) / ( (double) ell*(ell+1) );
			this_ell_coeffs[ell-m] = pow(-1,m) * integral_spherical_surface( integrand_minus, t_list, p_list ) / ( (double) ell*(ell+1) );
		}
		
		coeffs.push_back( this_ell_coeffs );
	}
	
	return coeffs;
}




std::vector< std::vector< std::complex<double> > > generate_vsh_series_C2( std::vector< std::vector< std::vector< std::complex<double> > > > A, std::vector<double> t_list, std::vector<double> p_list, int ell_max ){
	// Generate the C^{2,ell}_m coefficients in the vector spherical harmonic series of a vector function A(theta,phi) with no r-dependence.
	
	std::vector< std::vector< std::complex<double> > > coeffs { std::vector< std::complex<double> > { 0 } };	// (ell,m)=(0,0) coefficient is zero.
	
	int N_t = t_list.size();
	int N_p = p_list.size();
	
	for( int ell=1; ell<=ell_max; ell++ ){
		// The range of the index m depends on ell, so the second index of our 2D vectors has varying size and we cannot define the entire vector skeleton in one go.
		std::vector< std::complex<double> > this_ell_coeffs( 2*ell+1 );
		
		//----- m=0 coefficient -----
		std::vector< std::vector< std::complex<double> > > integrand( N_t, std::vector< std::complex<double> > ( N_p ) );
		for( int i=0; i<N_t; i++ ){
			for( int j=0; j<N_p; j++ ){
				integrand[i][j] = A[2][i][j] * ylm_dtheta( t_list[i], p_list[j], ell, 0 );
			}
		}
		this_ell_coeffs[ell] = integral_spherical_surface( integrand, t_list, p_list ) / ( (double) ell*(ell+1) );
		
		//----- m!=0 coefficients -----
		for( int m=1; m<=ell; m++ ){
			std::vector< std::vector< std::complex<double> > > integrand_plus ( N_t, std::vector< std::complex<double> > ( N_p ) );
			std::vector< std::vector< std::complex<double> > > integrand_minus( N_t, std::vector< std::complex<double> > ( N_p ) );
			for( int i=0; i<N_t; i++ ){
				for( int j=0; j<N_p; j++ ){
					std::complex<double> SH_dtheta = ylm_dtheta( t_list[i], p_list[j], ell, m );
					std::complex<double> SH_over_sintheta = ylm_over_sintheta( t_list[i], p_list[j], ell, m );
					integrand_plus [i][j] = A[2][i][j] * std::conj( SH_dtheta ) + A[1][i][j] * std::complex<double>{0,(double) m} * std::conj( SH_over_sintheta );
					integrand_minus[i][j] = A[2][i][j] *            SH_dtheta   + A[1][i][j] * std::complex<double>{0,(double)-m} *            SH_over_sintheta;
				}
			}
			this_ell_coeffs[ell+m] =             integral_spherical_surface( integrand_plus , t_list, p_list ) / ( (double) ell*(ell+1) );
			this_ell_coeffs[ell-m] = pow(-1,m) * integral_spherical_surface( integrand_minus, t_list, p_list ) / ( (double) ell*(ell+1) );
		}
		
		coeffs.push_back( this_ell_coeffs );
	}
	
	return coeffs;
}