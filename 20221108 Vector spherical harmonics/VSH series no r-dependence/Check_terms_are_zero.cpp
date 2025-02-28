// cd OneDrive\PhD\Codes\20221108 Vector spherical harmonics\VSH series no r-dependence
// g++ Check_terms_are_zero.cpp -o Check_terms_are_zero

/*
Fit a vector spherical harmonic series to a vector A(theta,phi).
New version in which the VSHs are precalculated and all fitting is done within the code, i.e. not using external functions, to avoid the read-time issues.
Determine whether
sum_ell ( sum_m ( A^(1) Phi - A(2) Psi ) ) = 0
Note that sum_ell( sum_m( A^1 Psi ) ) is the (1) part of A, and likewise for (2), but here we have cross-terms.
*/

#include <iostream>
#include <iomanip>
#include <chrono>
#include "../../Generate_Spherical_Harmonics.h"
#include "../Headers_vector_function/Vector_function_05.h"
#include "../../Numerical_Integration.h"
#include "../Headers_VSH/VSH.h"

int main(){
	
	auto time_start = std::chrono::high_resolution_clock::now();
	
	//----- Define variables -----
	int ell_max = 10;
	int N_t     = 300;
	int N_p     = 300;
	
	int i_test = 18;
	int j_test = 20;
	
	std::vector<double> t_list ( N_t );
	std::vector<double> p_list ( N_p );
	
	std::vector< std::vector< std::vector< std::complex<double> > > > A( 3, std::vector< std::vector< std::complex<double> > > ( N_t, std::vector< std::complex<double> > ( N_p ) ) );		// Precalculated values of the function that we wish to fit.
	std::vector< std::vector< std::vector< std::vector< std::complex<double> > > > > SH;				// Precalculated spherical harmonics
	std::vector< std::vector< std::vector< std::vector< std::complex<double> > > > > SH_dtheta;			// Precalculated spherical harmonic theta-derivatives
	std::vector< std::vector< std::vector< std::vector< std::complex<double> > > > > SH_over_sintheta;	// Precalculated spherical harmonics over sin(theta)
	std::vector< std::vector< std::vector< std::complex<double> > > > coeffs_VSH_series ( 3, ( std::vector< std::vector< std::complex<double> > > ( ell_max+1, ( std::vector< std::complex<double> > ( 2*ell_max+1 ) ) ) ) );	// Coefficients of VSH series to be calculated. Many of these are zero and redundant, but the only way to avoid is by building three separate lists. Since these are done in the mainspace, it will waste memory allocation even more.
	
	// Totals of absolute values of fitted coefficients in each direction. No physical interpretation but useful as a quick verification tool.
	std::complex<double> coeffs_total_r = 0;
	std::complex<double> coeffs_total_1 = 0;
	std::complex<double> coeffs_total_2 = 0;
	
	
	//----- Fail if the test indices are greater than the maximum number of points -----
	if( i_test > N_t-1 ){
		std::cout << "Test theta index (" << i_test << ") exceeds number of datapoints (" << N_t << "). Stopping the code to prevent undefined behaviour." << std::endl;
		return 1;
	}
	if( j_test > N_p-1 ){
		std::cout << "Test phi index (" << j_test << ") exceeds number of datapoints (" << N_p << "). Stopping the code to prevent undefined behaviour." << std::endl;
		return 1;
	}
	
	
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
	
	
	//----- Precalculate spherical harmonics -----
	for( int ell=0; ell<=ell_max; ell++ ){
		SH              .push_back( std::vector< std::vector< std::vector< std::complex<double> > > > ( 2*ell+1, std::vector< std::vector< std::complex<double> > > ( N_t, std::vector< std::complex<double> > ( N_p ) ) ) );
		SH_dtheta       .push_back( std::vector< std::vector< std::vector< std::complex<double> > > > ( 2*ell+1, std::vector< std::vector< std::complex<double> > > ( N_t, std::vector< std::complex<double> > ( N_p ) ) ) );
		SH_over_sintheta.push_back( std::vector< std::vector< std::vector< std::complex<double> > > > ( 2*ell+1, std::vector< std::vector< std::complex<double> > > ( N_t, std::vector< std::complex<double> > ( N_p ) ) ) );
		for( int m=-ell; m<=ell; m++ ){
			for( int i=0; i<N_t; i++ ){
				for( int j=0; j<N_p; j++ ){
					SH              [ell][ell+m][i][j] = ylm              ( t_list[i], p_list[j], ell, m );
					SH_dtheta       [ell][ell+m][i][j] = ylm_dtheta       ( t_list[i], p_list[j], ell, m );
					SH_over_sintheta[ell][ell+m][i][j] = ylm_over_sintheta( t_list[i], p_list[j], ell, m );
				}
			}
		}
	}
	//std::cout << SH[2][2+1][3][7] <<"\t"<< ylm(t_list[3],p_list[7],2,1) << std::endl;
	
	
	//----- Calculate VSH series -----
	for( int ell=0; ell<=ell_max; ell++ ){
		for( int m=-ell; m<=ell; m++ ){
			
			//--- Define arrays to store integrand values ---
			std::vector< std::vector< std::complex<double> > > integrand_r( N_t, std::vector< std::complex<double> > ( N_p ) );
			std::vector< std::vector< std::complex<double> > > integrand_1( N_t, std::vector< std::complex<double> > ( N_p ) );
			std::vector< std::vector< std::complex<double> > > integrand_2( N_t, std::vector< std::complex<double> > ( N_p ) );
			
			//--- Calculate the VSH series integrands by their definitions ---
			for( int i=0; i<N_t; i++ ){
				for( int j=0; j<N_p; j++ ){
					integrand_r[i][j] += A[0][i][j] * std::conj( SH[ell][ell+m][i][j] );
					integrand_1[i][j] += A[1][i][j] * std::conj( SH_dtheta[ell][ell+m][i][j] );
					integrand_2[i][j] += A[2][i][j] * std::conj( SH_dtheta[ell][ell+m][i][j] );
					if( m != 0 ){
						integrand_1[i][j] += A[2][i][j] * std::conj( std::complex<double>{0,(double) m} * SH_over_sintheta[ell][ell+m][i][j] );
						integrand_2[i][j] += A[1][i][j] * std::conj( std::complex<double>{0,(double)-m} * SH_over_sintheta[ell][ell+m][i][j] );
					}
				}
			}
			
			//--- Calculate the integral. Must be done in separate loop because we look ahead to (i+1,j+). ---
			for( int i=0; i<N_t-1; i++ ){
				for( int j=0; j<N_p-1; j++ ){
					double dA = 0.25 * ( cos(t_list[i]) - cos(t_list[i+1]) ) * ( p_list[j+1]-p_list[j] );
					coeffs_VSH_series[0][ell][ell+m] += ( integrand_r[i][j] + integrand_r[i][j+1] + integrand_r[i+1][j] + integrand_r[i+1][j+1] ) * dA;
					coeffs_VSH_series[1][ell][ell+m] += ( integrand_1[i][j] + integrand_1[i][j+1] + integrand_1[i+1][j] + integrand_1[i+1][j+1] ) * dA;
					coeffs_VSH_series[2][ell][ell+m] += ( integrand_2[i][j] + integrand_2[i][j+1] + integrand_2[i+1][j] + integrand_2[i+1][j+1] ) * dA;
				}
			}
			
			if( ell != 0 ){
				coeffs_VSH_series[1][ell][ell+m] /= (double) ell*(ell+1);
				coeffs_VSH_series[2][ell][ell+m] /= (double) ell*(ell+1);
			}
		}
	}
	
	
	//----- Output VSH series -----
	std::cout <<std::left<<std::setw(4)<< "ell" <<std::left<<std::setw(4)<< "m"   <<":\t"
	<<std::left<<std::setw(15)<< "Re(C^r)" <<std::left<<std::setw(15)<< "Im(C^r)" <<"\t/\t"
	<<std::left<<std::setw(15)<< "Re(C^1)" <<std::left<<std::setw(15)<< "Im(C^1)" <<"\t/\t"
	<<std::left<<std::setw(15)<< "Re(C^2)" <<std::left<<std::setw(15)<< "Im(C^2)" << std::endl;
			
	for( int ell=0; ell<=ell_max; ell++ ){
		for( int m=-ell; m<=ell; m++ ){
			std::cout <<std::left<<std::setw(4)<< ell <<std::left<<std::setw(4)<< m <<":\t"
			<<std::left<<std::setw(15)<< coeffs_VSH_series[0][ell][ell+m].real() <<std::left<<std::setw(15)<< coeffs_VSH_series[0][ell][ell+m].imag() <<"\t/\t"
			<<std::left<<std::setw(15)<< coeffs_VSH_series[1][ell][ell+m].real() <<std::left<<std::setw(15)<< coeffs_VSH_series[1][ell][ell+m].imag() <<"\t/\t"
			<<std::left<<std::setw(15)<< coeffs_VSH_series[2][ell][ell+m].real() <<std::left<<std::setw(15)<< coeffs_VSH_series[2][ell][ell+m].imag()	<< std::endl;

			coeffs_total_r += std::complex<double> { abs( coeffs_VSH_series[0][ell][ell+m].real() ), abs( coeffs_VSH_series[0][ell][ell+m].imag() ) };
			coeffs_total_1 += std::complex<double> { abs( coeffs_VSH_series[1][ell][ell+m].real() ), abs( coeffs_VSH_series[1][ell][ell+m].imag() ) };
			coeffs_total_2 += std::complex<double> { abs( coeffs_VSH_series[2][ell][ell+m].real() ), abs( coeffs_VSH_series[2][ell][ell+m].imag() ) };
		}
		std::cout << std::endl;
	}
	
	std::cout << "Coeffs total:\t" << coeffs_total_r <<"\t"<< coeffs_total_1 <<"\t"<< coeffs_total_2 <<"\t\t"<< coeffs_total_r+coeffs_total_1+coeffs_total_2 << std::endl;
			
	
	//----- Evaluate VSH series at the test point -----
	std::cout << "\nTest fit for (theta,phi) = (" << t_list[i_test] <<","<< p_list[j_test] << "):" << std::endl;
	
	std::vector< std::complex<double> > A_test_calc (3);
	
	for( int ell=0; ell<=ell_max; ell++ ){
		for( int m=-ell; m<=ell; m++ ){
			A_test_calc[0] += coeffs_VSH_series[0][ell][ell+m] * SH[ell][ell+m][i_test][j_test];
			A_test_calc[1] += coeffs_VSH_series[1][ell][ell+m] * SH_dtheta[ell][ell+m][i_test][j_test];
			A_test_calc[2] += coeffs_VSH_series[2][ell][ell+m] * SH_dtheta[ell][ell+m][i_test][j_test];
			if( m!=0 ){
				A_test_calc[1] += coeffs_VSH_series[2][ell][ell+m] * std::complex<double>{0,(double)-m} * SH_over_sintheta[ell][ell+m][i_test][j_test];
				A_test_calc[2] += coeffs_VSH_series[1][ell][ell+m] * std::complex<double>{0,(double) m} * SH_over_sintheta[ell][ell+m][i_test][j_test];
			}
		}
	}
	
	std::vector< std::complex<double> > A_test_exact { A[0][i_test][j_test], A[1][i_test][j_test], A[2][i_test][j_test] };
	
	std::cout <<std::left<<std::setw(7)<< "Cmpnt." <<std::left<<std::setw(28)<< "Calc."        <<std::left<<std::setw(28)<< "Exact"         << std::endl;
	std::cout <<std::left<<std::setw(7)<< "r"      <<std::left<<std::setw(28)<< A_test_calc[0] <<std::left<<std::setw(28)<< A_test_exact[0] << std::endl;
	std::cout <<std::left<<std::setw(7)<< "theta"  <<std::left<<std::setw(28)<< A_test_calc[1] <<std::left<<std::setw(28)<< A_test_exact[1] << std::endl;
	std::cout <<std::left<<std::setw(7)<< "phi"    <<std::left<<std::setw(28)<< A_test_calc[2] <<std::left<<std::setw(28)<< A_test_exact[2] << std::endl;
	
	
	//----- Output execution time and finish code -----
	auto time_stop = std::chrono::high_resolution_clock::now();
	double exec_time = std::chrono::duration_cast<std::chrono::nanoseconds>( time_stop - time_start ).count() * 1e-9;
	std::cout << "\nExecution time:\t" << floor(exec_time/60) << " m " << exec_time-60*floor(exec_time/60) << " s." << std::endl;
	std::cout << "\nFinished" << std::endl;
	
	
	//----- Calculate the sum -----
	std::vector< std::complex<double> > sum_to_check (3);
	
	double t = t_list[i_test];
	double p = p_list[j_test];
	
	for( int ell=0; ell<=ell_max; ell++ ){
		for( int m=-ell; m<=ell; m++ ){
			for( int i=1; i<=2; i++ ){
				sum_to_check[i] += coeffs_VSH_series[1][ell][ell+m]*VSH_Phi(t,p,ell,m)[i] - coeffs_VSH_series[2][ell][ell+m]*VSH_Psi(t,p,ell,m)[i];
			}
		}
	}
	std::cout << "Sum to check:\t" << sum_to_check[0] <<"\t"<< sum_to_check[1] <<"\t"<< sum_to_check[2] << std::endl;
}