// cd OneDrive\PhD\Codes\20221108 Vector spherical harmonics\VSH series axisymmetric
// g++ Fit_VSH_series_axisymmetric_dot_and_cross_products_01.cpp -o Fit_VSH_series_axisymmetric_dot_and_cross_products_01

/*
Fit a vector spherical harmonic series to an axisymmetric (no phi-dependence) vector A(theta).
*/

#include <iostream>
#include <iomanip>
#include <chrono>
#include "../../Generate_Associated_Legendre_Functions.h"
#include "../Headers_vector_function/Vector_function_11_axisymmetric_r_dependence.h"
#include "../Headers_vector_function/Vector_function_12_axisymmetric_r_dependence_b.h"
#include "../../Numerical_Integration.h"
#include "../../Vector_Operations.h"
#include "../../Numerical_Differentiation.h"

int main(){
	
	auto time_start = std::chrono::high_resolution_clock::now();
	
	//----- Define variables -----
	int n_max   =  2;	// Maximum order of Chebyshev polynomials.
	int ell_max = 5;	// Maximum order of associated Legendre functions.
	int N_r     = 20;	// Number of coordinates in each dimension.
	int N_t     = 10;
	int N_p     = 10;
	int i_test  = 6;	// Coordinate indices of a test point.
	int j_test  = 3;
	int k_test  = 8;
	
	std::vector<double> r_list ( N_r );
	std::vector<double> t_list ( N_t );
	std::vector<double> p_list ( N_p );
	std::vector<double> acos_r_list( N_r );
	std::vector<double> cos_t_list ( N_t );
	
	std::vector< std::vector< std::vector< std::vector<double> > > > A( 3, std::vector< std::vector< std::vector<double> > > ( N_r, std::vector< std::vector<double> > ( N_t, std::vector<double> ( N_p ) ) ) );		// Precalculated values of the function that we wish to fit.
	std::vector< std::vector< std::vector< std::vector<double> > > > B( 3, std::vector< std::vector< std::vector<double> > > ( N_r, std::vector< std::vector<double> > ( N_t, std::vector<double> ( N_p ) ) ) );
	std::vector< std::vector< std::vector<double> > > coeffs_VSH_series_A ( 3, std::vector< std::vector<double> > ( n_max+1, ( std::vector<double> ( ell_max+1 ) ) ) );	// Coefficients of VSH series to be calculated.
	std::vector< std::vector< std::vector<double> > > coeffs_VSH_series_B ( 3, std::vector< std::vector<double> > ( n_max+1, ( std::vector<double> ( ell_max+1 ) ) ) );
	std::vector< std::vector<double> > Tn( n_max  +1, ( std::vector<double> ( N_r ) ) );				// Precalculated Chebyshev polynomials of the first  kind.
	std::vector< std::vector<double> > Un( n_max  +1, ( std::vector<double> ( N_r ) ) );				// Precalculated Chebyshev polynomials of the second kind.
	std::vector< std::vector<double> > P0( ell_max+1, ( std::vector<double> ( N_t ) ) );				// Precalculated associated Legendre functions m=0 (i.e. Legendre polynomials).
	std::vector< std::vector<double> > P1( ell_max+1, ( std::vector<double> ( N_t ) ) );				// Precalculated associated Legendre functions m=1.
	
	// Totals of absolute values of fitted coefficients in each direction. No physical interpretation but useful as a quick verification tool.
	double coeffs_total_A_r = 0;
	double coeffs_total_A_1 = 0;
	double coeffs_total_A_2 = 0;
	double coeffs_total_B_r = 0;
	double coeffs_total_B_1 = 0;
	double coeffs_total_B_2 = 0;
	
	
	//----- Fail if the test indices are greater than the maximum number of points -----
	if( i_test > N_r-1 ){
		std::cout << "Test r index     (" << i_test << ") exceeds number of datapoints (" << N_r << "). Stopping the code to prevent undefined behaviour." << std::endl;
		return 1;
	}
	if( j_test > N_t-1 ){
		std::cout << "Test theta index (" << j_test << ") exceeds number of datapoints (" << N_t << "). Stopping the code to prevent undefined behaviour." << std::endl;
		return 1;
	}
	if( k_test > N_p-1 ){
		std::cout << "Test phi index   (" << k_test << ") exceeds number of datapoints (" << N_p << "). Stopping the code to prevent undefined behaviour." << std::endl;
		return 1;
	}
	
	
	//----- Build coordinate lists -----
	for( int i=0; i<N_r; i++ ){
		r_list[i] = -1 + i * 2.0 / ( (double) N_r - 1 );
		acos_r_list[i] = acos( r_list[i] );
	}
	for( int j=0; j<N_t; j++ ){
		t_list[j] = j * pi / ( (double) N_t - 1 );
		cos_t_list[j] = cos( t_list[j] );
	}
	for( int k=0; k<N_p; k++ ){
		p_list[k] = k * 2*pi / ( (double) N_p - 1 );
	}
	
	
	//----- Generate associated Legendre function coefficients -----
	generate_coeffs_associated_legendre( ell_max );
	
	
	//----- Precalculate vectors A and B to be fitted -----
	for( int i=0; i<N_r; i++ ){
		for( int j=0; j<N_t; j++ ){
			for( int k=0; k<N_p; k++ ){
				A[0][i][j][k] = A_r_function( r_list[i], t_list[j], p_list[k] );
				A[1][i][j][k] = A_t_function( r_list[i], t_list[j], p_list[k] );
				A[2][i][j][k] = A_p_function( r_list[i], t_list[j], p_list[k] );
				B[0][i][j][k] = B_r_function( r_list[i], t_list[j], p_list[k] );
				B[1][i][j][k] = B_t_function( r_list[i], t_list[j], p_list[k] );
				B[2][i][j][k] = B_p_function( r_list[i], t_list[j], p_list[k] );
			}
		}
	}
	
	
	//----- Precalculate Chebyshev polynomials and associated Legendre functions -----
	for( int n=0; n<=n_max; n++ ){
		for( int i=0; i<N_r; i++ ){
			Tn[n][i] = cos( (double) n * acos_r_list[i] );
			Un[n][i] = sin( (double) (n+1) * acos_r_list[i] ) / sin( acos_r_list[i] );
		}
	}
	//std::cout << "Chebyshev polynomial test:\t"<< r_list[12] <<"\t"<< Tn[2][12] <<"\t"<< Un[2][12] << std::endl;
	
	for( int ell=0; ell<=ell_max; ell++ ){
		for( int j=0; j<N_t; j++ ){
			P0[ell][j] = associated_legendre_function( cos_t_list[j], ell, 0 );
			P1[ell][j] = associated_legendre_function( cos_t_list[j], ell, 1 );
		}
	}
	
	
	//----- Calculate VSH series of A and output as we go -----
	std::cout <<std::left<<std::setw(4)<< "n"<< std::left<<std::setw(4)<< "ell" <<":\t"<< std::left<<std::setw(15)<< "A^r" <<"\t"<< std::left<<std::setw(15)<< "A^1" <<"\t"<< std::left<<std::setw(15)<< "A^2" << std::endl;
	
	for( int n=0; n<=n_max; n++ ){
		for( int ell=0; ell<=ell_max; ell++ ){
				
			//--- Define arrays to store integrand values ---
			std::vector< std::vector< std::vector<double> > > integrand_r( N_r, ( std::vector< std::vector<double> > ( N_t, std::vector<double> ( N_p ) ) ) );
			std::vector< std::vector< std::vector<double> > > integrand_1( N_r, ( std::vector< std::vector<double> > ( N_t, std::vector<double> ( N_p ) ) ) );
			std::vector< std::vector< std::vector<double> > > integrand_2( N_r, ( std::vector< std::vector<double> > ( N_t, std::vector<double> ( N_p ) ) ) );
			
			//--- Calculate the VSH series integrands by their definitions ---
			for( int i=0; i<N_r; i++ ){
				for( int j=0; j<N_t; j++ ){
					for( int k=0; k<N_p; k++ ){
						integrand_r[i][j][k] += A[0][i][j][k] * Tn[n][i] * P0[ell][j];
						integrand_1[i][j][k] += A[1][i][j][k] * Tn[n][i] * P1[ell][j];
						integrand_2[i][j][k] += A[2][i][j][k] * Tn[n][i] * P1[ell][j];
					}
				}
			}
		
			//--- Calculate the integral. Must be done in separate loop because we look ahead to (i+1,j+1,k+1). ---
			for( int i=0; i<N_r-1; i++ ){
				for( int j=0; j<N_t-1; j++ ){
					for( int k=0; k<N_p-1; k++ ){
						double dV = 2.0 / pi * sqrt( ( 2.0*ell+1.0 ) / ( 4.0*pi ) ) * 0.125 * ( acos_r_list[i] - acos_r_list[i+1] ) * ( cos_t_list[j] - cos_t_list[j+1] ) * ( p_list[k+1]-p_list[k] );
						coeffs_VSH_series_A[0][n][ell] += ( integrand_r[i][j][k] + integrand_r[i][j][k+1] + integrand_r[i][j+1][k] + integrand_r[i][j+1][k+1] + integrand_r[i+1][j][k] + integrand_r[i+1][j][k+1] + integrand_r[i+1][j+1][k] + integrand_r[i+1][j+1][k+1] ) * dV;
						coeffs_VSH_series_A[1][n][ell] += ( integrand_1[i][j][k] + integrand_1[i][j][k+1] + integrand_1[i][j+1][k] + integrand_1[i][j+1][k+1] + integrand_1[i+1][j][k] + integrand_1[i+1][j][k+1] + integrand_1[i+1][j+1][k] + integrand_1[i+1][j+1][k+1] ) * dV;
						coeffs_VSH_series_A[2][n][ell] += ( integrand_2[i][j][k] + integrand_2[i][j][k+1] + integrand_2[i][j+1][k] + integrand_2[i][j+1][k+1] + integrand_2[i+1][j][k] + integrand_2[i+1][j][k+1] + integrand_2[i+1][j+1][k] + integrand_2[i+1][j+1][k+1] ) * dV;
					}
				}
			}
			
			if( ell != 0 ){
				coeffs_VSH_series_A[1][n][ell] /= (double) ell*(ell+1);
				coeffs_VSH_series_A[2][n][ell] /= (double) ell*(ell+1);
			}
			
			//--- Output and increment coefficient totals (for comparison between code runs) ---
			std::cout <<std::left<<std::setw(4)<< n << std::left<<std::setw(4)<< ell <<":\t"<< std::left<<std::setw(15)<< coeffs_VSH_series_A[0][n][ell] <<"\t"<< std::left<<std::setw(15)<< coeffs_VSH_series_A[1][n][ell] <<"\t"<< std::left<<std::setw(15)<< coeffs_VSH_series_A[2][n][ell] << std::endl;

			coeffs_total_A_r += coeffs_VSH_series_A[0][n][ell];
			coeffs_total_A_1 += coeffs_VSH_series_A[1][n][ell];
			coeffs_total_A_2 += coeffs_VSH_series_A[2][n][ell];
		}
	}
	
	std::cout << "\nCoeffs total A:\t" << coeffs_total_A_r <<"\t"<< coeffs_total_A_1 <<"\t"<< coeffs_total_A_2 <<"\t\t"<< coeffs_total_A_r+coeffs_total_A_1+coeffs_total_A_2 << std::endl;
	
	
	//----- Calculate VSH series of B and output as we go -----
	std::cout <<std::left<<std::setw(4)<< "\nn"<< std::left<<std::setw(4)<< "ell" <<":\t"<< std::left<<std::setw(15)<< "B^r" <<"\t"<< std::left<<std::setw(15)<< "B^1" <<"\t"<< std::left<<std::setw(15)<< "B^2" << std::endl;
	
	for( int n=0; n<=n_max; n++ ){
		for( int ell=0; ell<=ell_max; ell++ ){
				
			//--- Define arrays to store integrand values ---
			std::vector< std::vector< std::vector<double> > > integrand_r( N_r, ( std::vector< std::vector<double> > ( N_t, std::vector<double> ( N_p ) ) ) );
			std::vector< std::vector< std::vector<double> > > integrand_1( N_r, ( std::vector< std::vector<double> > ( N_t, std::vector<double> ( N_p ) ) ) );
			std::vector< std::vector< std::vector<double> > > integrand_2( N_r, ( std::vector< std::vector<double> > ( N_t, std::vector<double> ( N_p ) ) ) );
			
			//--- Calculate the VSH series integrands by their definitions ---
			for( int i=0; i<N_r; i++ ){
				for( int j=0; j<N_t; j++ ){
					for( int k=0; k<N_p; k++ ){
						integrand_r[i][j][k] += B[0][i][j][k] * Tn[n][i] * P0[ell][j];
						integrand_1[i][j][k] += B[1][i][j][k] * Tn[n][i] * P1[ell][j];
						integrand_2[i][j][k] += B[2][i][j][k] * Tn[n][i] * P1[ell][j];
					}
				}
			}
		
			//--- Calculate the integral. Must be done in separate loop because we look ahead to (i+1,j+1,k+1). ---
			for( int i=0; i<N_r-1; i++ ){
				for( int j=0; j<N_t-1; j++ ){
					for( int k=0; k<N_p-1; k++ ){
						double dV = 2.0 / pi * sqrt( ( 2.0*ell+1.0 ) / ( 4.0*pi ) ) * 0.125 * ( acos_r_list[i] - acos_r_list[i+1] ) * ( cos_t_list[j] - cos_t_list[j+1] ) * ( p_list[k+1]-p_list[k] );
						coeffs_VSH_series_B[0][n][ell] += ( integrand_r[i][j][k] + integrand_r[i][j][k+1] + integrand_r[i][j+1][k] + integrand_r[i][j+1][k+1] + integrand_r[i+1][j][k] + integrand_r[i+1][j][k+1] + integrand_r[i+1][j+1][k] + integrand_r[i+1][j+1][k+1] ) * dV;
						coeffs_VSH_series_B[1][n][ell] += ( integrand_1[i][j][k] + integrand_1[i][j][k+1] + integrand_1[i][j+1][k] + integrand_1[i][j+1][k+1] + integrand_1[i+1][j][k] + integrand_1[i+1][j][k+1] + integrand_1[i+1][j+1][k] + integrand_1[i+1][j+1][k+1] ) * dV;
						coeffs_VSH_series_B[2][n][ell] += ( integrand_2[i][j][k] + integrand_2[i][j][k+1] + integrand_2[i][j+1][k] + integrand_2[i][j+1][k+1] + integrand_2[i+1][j][k] + integrand_2[i+1][j][k+1] + integrand_2[i+1][j+1][k] + integrand_2[i+1][j+1][k+1] ) * dV;
					}
				}
			}
			
			if( ell != 0 ){
				coeffs_VSH_series_B[1][n][ell] /= (double) ell*(ell+1);
				coeffs_VSH_series_B[2][n][ell] /= (double) ell*(ell+1);
			}
			
			//--- Output and increment coefficient totals (for comparison between code runs) ---
			std::cout <<std::left<<std::setw(4)<< n << std::left<<std::setw(4)<< ell <<":\t"<< std::left<<std::setw(15)<< coeffs_VSH_series_B[0][n][ell] <<"\t"<< std::left<<std::setw(15)<< coeffs_VSH_series_B[1][n][ell] <<"\t"<< std::left<<std::setw(15)<< coeffs_VSH_series_B[2][n][ell] << std::endl;

			coeffs_total_B_r += coeffs_VSH_series_B[0][n][ell];
			coeffs_total_B_1 += coeffs_VSH_series_B[1][n][ell];
			coeffs_total_B_2 += coeffs_VSH_series_B[2][n][ell];
		}
	}
	
	std::cout << "\nCoeffs total B:\t" << coeffs_total_B_r <<"\t"<< coeffs_total_B_1 <<"\t"<< coeffs_total_B_2 <<"\t\t"<< coeffs_total_B_r+coeffs_total_B_1+coeffs_total_B_2 << std::endl;
	
	
	//----- Evaluate VSH series at the test point -----
	std::vector<double> A_test_exact { A[0][i_test][j_test][k_test], A[1][i_test][j_test][k_test], A[2][i_test][j_test][k_test] };
	std::vector<double> B_test_exact { B[0][i_test][j_test][k_test], B[1][i_test][j_test][k_test], B[2][i_test][j_test][k_test] };
	std::vector<double> A_test_calc (3);
	std::vector<double> B_test_calc (3);
	
	//--- n=0 terms ---
	for( int ell=0; ell<=ell_max; ell++ ){
		double factor = sqrt( ( 2.0*ell+1.0 ) / ( 4.0*pi ) );
		A_test_calc[0] += coeffs_VSH_series_A[0][0][ell] * 0.5 * factor * P0[ell][j_test];
		A_test_calc[1] += coeffs_VSH_series_A[1][0][ell] * 0.5 * factor * P1[ell][j_test];
		A_test_calc[2] += coeffs_VSH_series_A[2][0][ell] * 0.5 * factor * P1[ell][j_test];
		B_test_calc[0] += coeffs_VSH_series_B[0][0][ell] * 0.5 * factor * P0[ell][j_test];
		B_test_calc[1] += coeffs_VSH_series_B[1][0][ell] * 0.5 * factor * P1[ell][j_test];
		B_test_calc[2] += coeffs_VSH_series_B[2][0][ell] * 0.5 * factor * P1[ell][j_test];
	}
	//--- n>0 terms ---
	for( int n=1; n<=n_max; n++ ){
		for( int ell=0; ell<=ell_max; ell++ ){
			double factor = sqrt( ( 2.0*ell+1.0 ) / ( 4.0*pi ) );
			A_test_calc[0] += coeffs_VSH_series_A[0][n][ell] * factor * Tn[n][i_test] * P0[ell][j_test];
			A_test_calc[1] += coeffs_VSH_series_A[1][n][ell] * factor * Tn[n][i_test] * P1[ell][j_test];
			A_test_calc[2] += coeffs_VSH_series_A[2][n][ell] * factor * Tn[n][i_test] * P1[ell][j_test];
			B_test_calc[0] += coeffs_VSH_series_B[0][n][ell] * factor * Tn[n][i_test] * P0[ell][j_test];
			B_test_calc[1] += coeffs_VSH_series_B[1][n][ell] * factor * Tn[n][i_test] * P1[ell][j_test];
			B_test_calc[2] += coeffs_VSH_series_B[2][n][ell] * factor * Tn[n][i_test] * P1[ell][j_test];
		}
	}
	
	std::cout << "\nTest fit at (r,theta,phi) = (" << r_list[i_test] <<","<< t_list[j_test] <<","<< p_list[k_test] << "):\n" << std::endl;
	std::cout <<std::left<<std::setw(7)<< "Cmpnt." <<std::left<<std::setw(28)<< "A Calc."      <<std::left<<std::setw(28)<< "A Exact"       <<std::left<<std::setw(28)<< "B Calc."      <<std::left<<std::setw(28)<< "B Exact"       << std::endl;
	std::cout <<std::left<<std::setw(7)<< "r"      <<std::left<<std::setw(28)<< A_test_calc[0] <<std::left<<std::setw(28)<< A_test_exact[0] <<std::left<<std::setw(28)<< B_test_calc[0] <<std::left<<std::setw(28)<< B_test_exact[0] << std::endl;
	std::cout <<std::left<<std::setw(7)<< "theta"  <<std::left<<std::setw(28)<< A_test_calc[1] <<std::left<<std::setw(28)<< A_test_exact[1] <<std::left<<std::setw(28)<< B_test_calc[1] <<std::left<<std::setw(28)<< B_test_exact[1] << std::endl;
	std::cout <<std::left<<std::setw(7)<< "phi"    <<std::left<<std::setw(28)<< A_test_calc[2] <<std::left<<std::setw(28)<< A_test_exact[2] <<std::left<<std::setw(28)<< B_test_calc[2] <<std::left<<std::setw(28)<< B_test_exact[2] << std::endl;
	
	
	//----- Calculate A dot B -----
	double r = r_list[i_test];
	double t = t_list[j_test];
	double p = p_list[k_test];
	double A_dot_B_exact = A_test_exact[0]*B_test_exact[0] + A_test_exact[1]*B_test_exact[1] + A_test_exact[2]*B_test_exact[2];
	//double A_dot_B_exact_2 = A_r_function(r,t,p)*B_r_function(r,t,p) + A_t_function(r,t,p)*B_t_function(r,t,p) + A_p_function(r,t,p)*B_p_function(r,t,p);	// It's correct. Keeping in case future checks needed.
	double A_dot_B_test  = A_test_calc[0]*B_test_calc[0] + A_test_calc[1]*B_test_calc[1] + A_test_calc[2]*B_test_calc[2];
	std::cout << "Code gets here" << std::endl;
	double A_dot_B_calc = 0;
	for( int ell_1=0; ell_1<=ell_max; ell_1++ ){
		
		//--- V1 ---
		double A_r = 0.5 * coeffs_VSH_series_A[0][0][ell_1];
		double A_1 = 0.5 * coeffs_VSH_series_A[1][0][ell_1];
		double A_2 = 0.5 * coeffs_VSH_series_A[2][0][ell_1];
		for( int n=1; n<=n_max; n++ ){
			A_r += coeffs_VSH_series_A[0][n][ell_1] * Tn[n][i_test];
			A_1 += coeffs_VSH_series_A[1][n][ell_1] * Tn[n][i_test];
			A_2 += coeffs_VSH_series_A[2][n][ell_1] * Tn[n][i_test];
		}
		
		for( int ell_2=0; ell_2<=ell_max; ell_2++ ){	
			double B_r = 0.5 * coeffs_VSH_series_B[0][0][ell_2];
			double B_1 = 0.5 * coeffs_VSH_series_B[1][0][ell_2];
			double B_2 = 0.5 * coeffs_VSH_series_B[2][0][ell_2];
			for( int n=1; n<=n_max; n++ ){				
				B_r += coeffs_VSH_series_B[0][n][ell_2] * Tn[n][i_test];
				B_1 += coeffs_VSH_series_B[1][n][ell_2] * Tn[n][i_test];
				B_2 += coeffs_VSH_series_B[2][n][ell_2] * Tn[n][i_test];
			}
			
			double term = A_r * B_r * P0[ell_1][j_test] * P0[ell_2][j_test];
			term += ( A_1*B_1 + A_2*B_2 ) * P1[ell_1][j_test] * P1[ell_2][j_test];
			A_dot_B_calc += term * sqrt( (double)(2.0*ell_1+1.0) * (double)(2.0*ell_2+1.0) ) / ( 4.0 * pi );
		}
	}
	
	std::cout << "\nA dot B:\t" << A_dot_B_exact <<"\t"<< A_dot_B_test <<"\t"<< A_dot_B_calc << std::endl;
	
	
	//----- Calculate A cross B -----
	std::vector<double> A_cross_B_exact       = cross_product( A_test_exact, B_test_exact );
	std::vector<double> A_cross_B_calc_exact  = cross_product( A_test_calc , B_test_calc  );
	std::cout << "\nA cross B:" << std::endl;
	print_vector( A_cross_B_exact );
	print_vector( A_cross_B_calc_exact );
	
	std::vector<double> A_cross_B_calc (3);
	
	//--- V1 ---
	/*for( int ell_1=0; ell_1<=ell_max; ell_1++ ){
		
		double A_r = 0.5 * coeffs_VSH_series_A[0][0][ell_1];
		double A_1 = 0.5 * coeffs_VSH_series_A[1][0][ell_1];
		double A_2 = 0.5 * coeffs_VSH_series_A[2][0][ell_1];
		for( int n=1; n<=n_max; n++ ){
			A_r += coeffs_VSH_series_A[0][n][ell_1] * Tn[n][i_test];
			A_1 += coeffs_VSH_series_A[1][n][ell_1] * Tn[n][i_test];
			A_2 += coeffs_VSH_series_A[2][n][ell_1] * Tn[n][i_test];
		}
		
		for( int ell_2=0; ell_2<=ell_max; ell_2++ ){	
			double B_r = 0.5 * coeffs_VSH_series_B[0][0][ell_2];
			double B_1 = 0.5 * coeffs_VSH_series_B[1][0][ell_2];
			double B_2 = 0.5 * coeffs_VSH_series_B[2][0][ell_2];
			for( int n=1; n<=n_max; n++ ){				
				B_r += coeffs_VSH_series_B[0][n][ell_2] * Tn[n][i_test];
				B_1 += coeffs_VSH_series_B[1][n][ell_2] * Tn[n][i_test];
				B_2 += coeffs_VSH_series_B[2][n][ell_2] * Tn[n][i_test];
			}
			
			double factor = pow( 4.0*pi, -1 ) * sqrt( (double)(2.0*ell_1+1.0) * (double)(2.0*ell_2+1.0) );
			
			A_cross_B_calc[0] += factor * ( A_1*B_2 - A_2*B_1 ) * P1[ell_1][j_test] * P1[ell_2][j_test];
			A_cross_B_calc[1] += factor * ( A_2*B_r*P1[ell_1][j_test]*P0[ell_2][j_test] - A_r*B_2*P0[ell_1][j_test]*P1[ell_2][j_test] );
			A_cross_B_calc[2] += factor * ( A_r*B_1*P0[ell_1][j_test]*P1[ell_2][j_test] - A_1*B_r*P1[ell_1][j_test]*P0[ell_2][j_test] );
		}
	}*/
	
	//--- V2 (Separate the ell_2=0 term) ---
	/*for( int ell_1=0; ell_1<=ell_max; ell_1++ ){
		
		double A_r = 0.5 * coeffs_VSH_series_A[0][0][ell_1];
		double A_1 = 0.5 * coeffs_VSH_series_A[1][0][ell_1];
		double A_2 = 0.5 * coeffs_VSH_series_A[2][0][ell_1];
		for( int n=1; n<=n_max; n++ ){
			A_r += coeffs_VSH_series_A[0][n][ell_1] * Tn[n][i_test];
			A_1 += coeffs_VSH_series_A[1][n][ell_1] * Tn[n][i_test];
			A_2 += coeffs_VSH_series_A[2][n][ell_1] * Tn[n][i_test];
		}
		
		//--- ell_2 = 0 ---
		double B_r0 = 0.5 * coeffs_VSH_series_B[0][0][0];
		for( int n=1; n<=n_max; n++ ){				
			B_r0 += coeffs_VSH_series_B[0][n][0] * Tn[n][i_test];
		}
		double factor_0 = sqrt( ( 2.0*ell_1 + 1 ) / ( 4.0*pi ) ) * pow(4.0*pi,-0.5);
		A_cross_B_calc[1] += factor_0 * P1[ell_1][j_test] * B_r0 *  A_2;
		A_cross_B_calc[2] += factor_0 * P1[ell_1][j_test] * B_r0 * -A_1;
		
		
		//--- ell_2 > 0 ---
		for( int ell_2=1; ell_2<=ell_max; ell_2++ ){	
			double B_r = 0.5 * coeffs_VSH_series_B[0][0][ell_2];
			double B_1 = 0.5 * coeffs_VSH_series_B[1][0][ell_2];
			double B_2 = 0.5 * coeffs_VSH_series_B[2][0][ell_2];
			for( int n=1; n<=n_max; n++ ){				
				B_r += coeffs_VSH_series_B[0][n][ell_2] * Tn[n][i_test];
				B_1 += coeffs_VSH_series_B[1][n][ell_2] * Tn[n][i_test];
				B_2 += coeffs_VSH_series_B[2][n][ell_2] * Tn[n][i_test];
			}
			
			double factor = pow( 4.0*pi, -1 ) * sqrt( (double)(2.0*ell_1+1.0) * (double)(2.0*ell_2+1.0) );
			
			A_cross_B_calc[0] += factor * ( A_1*B_2 - A_2*B_1 ) * P1[ell_1][j_test] * P1[ell_2][j_test];
			A_cross_B_calc[1] += factor * ( A_2*B_r*P1[ell_1][j_test]*P0[ell_2][j_test] - A_r*B_2*P0[ell_1][j_test]*P1[ell_2][j_test] );
			A_cross_B_calc[2] += factor * ( A_r*B_1*P0[ell_1][j_test]*P1[ell_2][j_test] - A_1*B_r*P1[ell_1][j_test]*P0[ell_2][j_test] );
		}
	}*/
	
	//--- V3 (Separate the ell_1=0 and ell_2=0 terms) ---
	//--- ell_1 = 0 ---
	double A_r0 = 0.5 * coeffs_VSH_series_A[0][0][0];
	for( int n=1; n<=n_max; n++ ){
		A_r0 += coeffs_VSH_series_A[0][n][0] * Tn[n][i_test];
	}
	
	for( int ell_2=1; ell_2<=ell_max; ell_2++ ){	
		double B_1 = 0.5 * coeffs_VSH_series_B[1][0][ell_2];
		double B_2 = 0.5 * coeffs_VSH_series_B[2][0][ell_2];
		for( int n=1; n<=n_max; n++ ){				
			B_1 += coeffs_VSH_series_B[1][n][ell_2] * Tn[n][i_test];
			B_2 += coeffs_VSH_series_B[2][n][ell_2] * Tn[n][i_test];
		}
		
		double factor = sqrt(2.0*ell_2+1) * pow(4.0*pi,-1);
		A_cross_B_calc[1] += factor * P1[ell_2][j_test] * A_r0 * -B_2;
		A_cross_B_calc[2] += factor * P1[ell_2][j_test] * A_r0 *  B_1;
	}
	
	//--- ell_1 > 0 ---
	for( int ell_1=1; ell_1<=ell_max; ell_1++ ){
		
		double A_r = 0.5 * coeffs_VSH_series_A[0][0][ell_1];
		double A_1 = 0.5 * coeffs_VSH_series_A[1][0][ell_1];
		double A_2 = 0.5 * coeffs_VSH_series_A[2][0][ell_1];
		for( int n=1; n<=n_max; n++ ){
			A_r += coeffs_VSH_series_A[0][n][ell_1] * Tn[n][i_test];
			A_1 += coeffs_VSH_series_A[1][n][ell_1] * Tn[n][i_test];
			A_2 += coeffs_VSH_series_A[2][n][ell_1] * Tn[n][i_test];
		}
		
		//--- ell_2 = 0 ---
		double B_r0 = 0.5 * coeffs_VSH_series_B[0][0][0];
		for( int n=1; n<=n_max; n++ ){				
			B_r0 += coeffs_VSH_series_B[0][n][0] * Tn[n][i_test];
		}
		double factor_0 = sqrt( ( 2.0*ell_1 + 1 ) / ( 4.0*pi ) ) * pow(4.0*pi,-0.5);
		A_cross_B_calc[1] += factor_0 * P1[ell_1][j_test] * B_r0 *  A_2;
		A_cross_B_calc[2] += factor_0 * P1[ell_1][j_test] * B_r0 * -A_1;
		
		
		//--- ell_2 > 0 ---
		for( int ell_2=1; ell_2<=ell_max; ell_2++ ){	
			double B_r = 0.5 * coeffs_VSH_series_B[0][0][ell_2];
			double B_1 = 0.5 * coeffs_VSH_series_B[1][0][ell_2];
			double B_2 = 0.5 * coeffs_VSH_series_B[2][0][ell_2];
			for( int n=1; n<=n_max; n++ ){				
				B_r += coeffs_VSH_series_B[0][n][ell_2] * Tn[n][i_test];
				B_1 += coeffs_VSH_series_B[1][n][ell_2] * Tn[n][i_test];
				B_2 += coeffs_VSH_series_B[2][n][ell_2] * Tn[n][i_test];
			}
			
			double factor = pow( 4.0*pi, -1 ) * sqrt( (double)(2.0*ell_1+1.0) * (double)(2.0*ell_2+1.0) );
			
			A_cross_B_calc[0] += factor * ( A_1*B_2 - A_2*B_1 ) * P1[ell_1][j_test] * P1[ell_2][j_test];
			A_cross_B_calc[1] += factor * ( A_2*B_r*P1[ell_1][j_test]*P0[ell_2][j_test] - A_r*B_2*P0[ell_1][j_test]*P1[ell_2][j_test] );
			A_cross_B_calc[2] += factor * ( A_r*B_1*P0[ell_1][j_test]*P1[ell_2][j_test] - A_1*B_r*P1[ell_1][j_test]*P0[ell_2][j_test] );
		}
	}
	
	print_vector( A_cross_B_calc );
	
	
	//----- Output execution time and finish code -----
	auto time_stop = std::chrono::high_resolution_clock::now();
	double exec_time = std::chrono::duration_cast<std::chrono::nanoseconds>( time_stop - time_start ).count() * 1e-9;
	std::cout << "\nExecution time:\t" << floor(exec_time/60) << " m " << exec_time-60*floor(exec_time/60) << " s." << std::endl;
	return 0;
}