// VSH_series_calculation_02.h
// Calculate the VSH series coefficients for a vector A with no r-dependence.
// As code _01, but allowing for A to be complex.

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
					std::vector< std::complex<double> > A_at_this_point { A[0][i][j], A[1][i][j], A[2][i][j] };
					integrand[i][j] = dot_product_cart( A_at_this_point, vector_conj( VSH_Y( t_list[i], p_list[j], ell, m ) ) );
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
		for( int m=-ell; m<=ell; m++ ){
	
			//----- Precalculate the integrand -----
			std::vector< std::vector< std::complex<double> > > integrand( N_t, std::vector< std::complex<double> > ( N_p ) );
			for( int i=0; i<N_t; i++ ){
				for( int j=0; j<N_p; j++ ){
					std::vector< std::complex<double> > A_at_this_point { A[0][i][j], A[1][i][j], A[2][i][j] };
					integrand[i][j] = dot_product_cart( A_at_this_point, vector_conj( VSH_Psi( t_list[i], p_list[j], ell, m ) ) );
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
					std::vector< std::complex<double> > A_at_this_point { A[0][i][j], A[1][i][j], A[2][i][j] };
					integrand[i][j] = dot_product_cart( A_at_this_point, vector_conj( VSH_Phi( t_list[i], p_list[j], ell, m ) ) );
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