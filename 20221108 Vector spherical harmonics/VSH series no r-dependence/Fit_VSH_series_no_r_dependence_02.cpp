// cd OneDrive\PhD\Codes\20221108 Vector spherical harmonics\VSH series no r-dependence
// g++ Fit_VSH_series_no_r_dependence_02.cpp -o Fit_VSH_series_no_r_dependence_02

/*
Fit a vector spherical harmonic series to a vector A(theta,phi).
In this version, we precalculate the spherical harmonics and their derivatives instead of using functions to evaluate the VSHs.
*/

#include <iostream>
#include <iomanip>
#include "../../Generate_Spherical_Harmonics.h"
#include "../Headers_VSH/VSH_Y_01.h"
#include "../Headers_VSH/VSH_Psi_02.h"
#include "../Headers_VSH/VSH_Phi_01.h"
#include "../../Vector_Operations.h"
#include "../../Numerical_Integration.h"
#include "../Headers_vector_function/Vector_function_01.h"

int main(){
	
	//----- Define variables -----
	int ell_max = 1;
	int N_t     = 100;
	int N_p     = 100;
	
	std::vector<double> t_list ( N_t );
	std::vector<double> p_list ( N_p );
	
	std::vector< std::vector< std::vector<double> > > A( 3, std::vector< std::vector<double> > ( N_t, std::vector<double> ( N_p ) ) );		// Precalculated values of the function that we wish to fit.
	std::vector< std::vector< std::vector< std::vector< std::complex< double > > > > > sh;																	// Precalculated spherical harmonics.
	std::vector< std::vector< std::vector< std::vector< std::complex< double > > > > > sh_dt;																	// Precalcualted spherical harmonic theta-derivatives.
	
	
	//----- Build coordinate lists -----
	for( int i=0; i<N_t; i++ ){
		t_list[i] = i * pi / ( (double) N_t - 1 );
	}
	for( int j=0; j<N_p; j++ ){
		p_list[j] = j * 2*pi / ( (double) N_p - 1 );
	}
	
	
	//----- Generate spherical harmonic coefficients -----
	generate_coeffs_spherical_harmonic( ell_max );
	
	
	//----- Precalculate vector A to be fitted -----
	for( int i=0; i<N_t; i++ ){
		for( int j=0; j<N_p; j++ ){
			A[0][i][j] = A_r_function( t_list[i], p_list[j] );
			A[1][i][j] = A_t_function( t_list[i], p_list[j] );
			A[2][i][j] = A_p_function( t_list[i], p_list[j] );
		}
	}
	
	
	//----- Precalculate the spherical harmonics and their theta-derivatives -----
	for( int ell=0; ell<=ell_max; ell++ ){
		// The range of the index m depends on ell, so the second index of our 4D vectors has varying size and we cannot define the entire vector skeleton in one go.
		std::vector< std::vector< std::vector< std::complex<double> > > > this_ell_evaluated_sh    ( 2*ell+1, std::vector< std::vector< std::complex<double> > > ( N_t, std::vector< std::complex<double> > ( N_p ) ) );
		std::vector< std::vector< std::vector< std::complex<double> > > > this_ell_evaluated_sh_dt ( 2*ell+1, std::vector< std::vector< std::complex<double> > > ( N_t, std::vector< std::complex<double> > ( N_p ) ) );
		for( int m=-ell; m<=ell; m++ ){
			for( int i=0; i<N_t; i++ ){
				for( int j=0; j<N_p; j++ ){
					this_ell_evaluated_sh   [ell+m][i][j] = ylm       ( t_list[i], p_list[j], ell, m );
					this_ell_evaluated_sh_dt[ell+m][i][j] = ylm_dtheta( t_list[i], p_list[j], ell, m );
				}
			}
		}
		sh   .push_back( this_ell_evaluated_sh    );
		sh_dt.push_back( this_ell_evaluated_sh_dt );
	}
	
	std::cout << ylm       ( t_list[3], p_list[5], 1, 0 ) <<"\t"<< sh   [1][1+0][3][5] << std::endl;	// Test calculations and indexing are correct.
	std::cout << ylm_dtheta( t_list[3], p_list[5], 1, 0 ) <<"\t"<< sh_dt[1][1+0][3][5] << std::endl;	// Test calculations and indexing are correct.
	
	std::cout << "Finished" << std::endl;
}