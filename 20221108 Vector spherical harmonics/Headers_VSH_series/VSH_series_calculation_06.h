// VSH_series_calculation_06.h
// Calculate the VSH series coefficients for a vector A with no r-dependence.
// Separate the m=0 term to avoid the if-statement, then recombine the two additions for integrand[i][j]. Lose the unnecessary coeff variable.
// For now, only change the C_1 coefficient.

#include "../../Vector_Operations.h"
#include "../Headers_VSH/VSH_Y_01.h"
#include "../Headers_VSH/VSH_Psi_02.h"
#include "../Headers_VSH/VSH_Phi_01.h"
#include "../../Numerical_Integration.h"

std::vector< std::vector< std::complex<double> > > generate_vsh_series_Cr( std::vector< std::vector< std::vector< std::complex<double> > > > A, std::vector<double> t_list, std::vector<double> p_list, int ell_max ){
	// Generate the C^{r,ell}_m coefficients in the vector spherical harmonic series of a vector function A(theta,phi) with no r-dependence.
	
	std::vector< std::vector< std::complex<double> > > coeffs;
	
	int N_t = t_list.size();
	int N_p = p_list.size();
	
	for( int ell=0; ell<=ell_max; ell++ ){
		// The range of the index m depends on ell, so the second index of our 2D vectors has varying size and we cannot define the entire vector skeleton in one go.
		std::vector< std::complex<double> > this_ell_coeffs( 2*ell+1 );
		for( int m=-ell; m<=ell; m++ ){
	
			//----- Precalculate the integrand -----
			std::vector< std::vector< std::complex<double> > > integrand( N_t, std::vector< std::complex<double> > ( N_p ) );
			for( int i=0; i<N_t; i++ ){
				for( int j=0; j<N_p; j++ ){
					integrand[i][j] = A[0][i][j] * std::conj( ylm( t_list[i], p_list[j], ell, m ) );
				}
			}
			
			//----- Perform the integral -----
			std::complex<double> coeff = integral_spherical_surface( integrand, t_list, p_list );
			
			//----- Append the coefficient -----
			this_ell_coeffs[ell+m] = coeff;
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
		
		//----- m<0 coefficient -----
		for( int m=-ell; m<=-1; m++ ){
			std::vector< std::vector< std::complex<double> > > integrand( N_t, std::vector< std::complex<double> > ( N_p ) );
			for( int i=0; i<N_t; i++ ){
				for( int j=0; j<N_p; j++ ){
					integrand[i][j] = A[1][i][j] * std::conj( ylm_dtheta( t_list[i], p_list[j], ell, m ) ) + A[2][i][j] * std::complex<double>{0,-(double)m} * std::conj( ylm_over_sintheta( t_list[i], p_list[j], ell, m ) );
				}
			}
			this_ell_coeffs[ell+m] = integral_spherical_surface( integrand, t_list, p_list ) / ( (double) ell*(ell+1) );
		}
		
		//----- m>0 coefficient -----
		for( int m=1; m<=ell; m++ ){
			std::vector< std::vector< std::complex<double> > > integrand( N_t, std::vector< std::complex<double> > ( N_p ) );
			for( int i=0; i<N_t; i++ ){
				for( int j=0; j<N_p; j++ ){
					integrand[i][j] = A[1][i][j] * std::conj( ylm_dtheta( t_list[i], p_list[j], ell, m ) ) + A[2][i][j] * std::complex<double>{0,-(double)m} * std::conj( ylm_over_sintheta( t_list[i], p_list[j], ell, m ) );
				}
			}
			this_ell_coeffs[ell+m] = integral_spherical_surface( integrand, t_list, p_list ) / ( (double) ell*(ell+1) );
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
		for( int m=-ell; m<=ell; m++ ){
	
			//----- Precalculate the integrand -----
			std::vector< std::vector< std::complex<double> > > integrand( N_t, std::vector< std::complex<double> > ( N_p ) );
			for( int i=0; i<N_t; i++ ){
				for( int j=0; j<N_p; j++ ){
					integrand[i][j] = A[2][i][j] * std::conj( ylm_dtheta( t_list[i], p_list[j], ell, m ) );
					if( m != 0 ){
						integrand[i][j] += A[1][i][j] * std::complex<double>{0,(double)m} * std::conj( ylm_over_sintheta( t_list[i], p_list[j], ell, m ) );
					}
				}
			}
			
			//----- Perform the integral -----
			std::complex<double> coeff = integral_spherical_surface( integrand, t_list, p_list );
			
			//----- Append the coefficient -----
			this_ell_coeffs[ell+m] = coeff / ( (double) ell*(ell+1) );
		}
		
		coeffs.push_back( this_ell_coeffs );
	}
	
	return coeffs;
}