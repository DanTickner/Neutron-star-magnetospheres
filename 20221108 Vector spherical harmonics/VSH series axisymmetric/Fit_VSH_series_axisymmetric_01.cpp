// cd OneDrive\PhD\Codes\20221108 Vector spherical harmonics\VSH series axisymmetric
// g++ Fit_VSH_series_axisymmetric_01.cpp -o Fit_VSH_series_axisymmetric_01

/*
Fit a vector spherical harmonic series to an axisymmetric (no phi-dependence) vector A(theta).
Calculate also the divergence and curl.
*/

#include <iostream>
#include <iomanip>
#include <chrono>
#include "../../Generate_Associated_Legendre_Functions.h"
#include "../Headers_vector_function/Vector_function_08_axisymmetric_r_dependence.h"
#include "../../Numerical_Integration.h"
#include "../../Vector_Operations.h"
#include "../../Numerical_Differentiation.h"

int main(){
	
	auto time_start = std::chrono::high_resolution_clock::now();
	
	//----- Define variables -----
	int n_max   =  2;	// Maximum order of Chebyshev polynomials.
	int ell_max = 10;	// Maximum order of associated Legendre functions.
	int N_r     = 200;	// Number of coordinates in each dimension.
	int N_t     = 100;
	int N_p     = 100;
	int i_test  = 64;	// Coordinate indices of a test point.
	int j_test  = 38;
	int k_test  = 80;
	
	std::vector<double> r_list ( N_r );
	std::vector<double> t_list ( N_t );
	std::vector<double> p_list ( N_p );
	std::vector<double> acos_r_list( N_r );
	std::vector<double> cos_t_list ( N_t );
	
	std::vector< std::vector< std::vector< std::vector<double> > > > A( 3, std::vector< std::vector< std::vector<double> > > ( N_r, std::vector< std::vector<double> > ( N_t, std::vector<double> ( N_p ) ) ) );		// Precalculated values of the function that we wish to fit.
	std::vector< std::vector< std::vector<double> > > coeffs_VSH_series ( 3, std::vector< std::vector<double> > ( n_max+1, ( std::vector<double> ( ell_max+1 ) ) ) );	// Coefficients of VSH series to be calculated.
	std::vector< std::vector<double> > Tn( n_max  +1, ( std::vector<double> ( N_r ) ) );				// Precalculated Chebyshev polynomials of the first  kind.
	std::vector< std::vector<double> > Un( n_max  +1, ( std::vector<double> ( N_r ) ) );				// Precalculated Chebyshev polynomials of the second kind.
	std::vector< std::vector<double> > P0( ell_max+1, ( std::vector<double> ( N_t ) ) );				// Precalculated associated Legendre functions m=0 (i.e. Legendre polynomials).
	std::vector< std::vector<double> > P1( ell_max+1, ( std::vector<double> ( N_t ) ) );				// Precalculated associated Legendre functions m=1.
	
	// Totals of absolute values of fitted coefficients in each direction. No physical interpretation but useful as a quick verification tool.
	double coeffs_total_r = 0;
	double coeffs_total_1 = 0;
	double coeffs_total_2 = 0;
	
	
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
	
	
	//----- Precalculate vector A to be fitted -----
	for( int i=0; i<N_r; i++ ){
		for( int j=0; j<N_t; j++ ){
			for( int k=0; k<N_p; k++ ){
				A[0][i][j][k] = A_r_function( r_list[i], t_list[j], p_list[k] );
				A[1][i][j][k] = A_t_function( r_list[i], t_list[j], p_list[k] );
				A[2][i][j][k] = A_p_function( r_list[i], t_list[j], p_list[k] );
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
	
	
	//----- Calculate VSH series and output as we go -----
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
						coeffs_VSH_series[0][n][ell] += ( integrand_r[i][j][k] + integrand_r[i][j][k+1] + integrand_r[i][j+1][k] + integrand_r[i][j+1][k+1] + integrand_r[i+1][j][k] + integrand_r[i+1][j][k+1] + integrand_r[i+1][j+1][k] + integrand_r[i+1][j+1][k+1] ) * dV;
						coeffs_VSH_series[1][n][ell] += ( integrand_1[i][j][k] + integrand_1[i][j][k+1] + integrand_1[i][j+1][k] + integrand_1[i][j+1][k+1] + integrand_1[i+1][j][k] + integrand_1[i+1][j][k+1] + integrand_1[i+1][j+1][k] + integrand_1[i+1][j+1][k+1] ) * dV;
						coeffs_VSH_series[2][n][ell] += ( integrand_2[i][j][k] + integrand_2[i][j][k+1] + integrand_2[i][j+1][k] + integrand_2[i][j+1][k+1] + integrand_2[i+1][j][k] + integrand_2[i+1][j][k+1] + integrand_2[i+1][j+1][k] + integrand_2[i+1][j+1][k+1] ) * dV;
					}
				}
			}
			
			if( ell != 0 ){
				coeffs_VSH_series[1][n][ell] /= (double) ell*(ell+1);
				coeffs_VSH_series[2][n][ell] /= (double) ell*(ell+1);
			}
			
			//--- Output and increment coefficient totals (for comparison between code runs) ---
			std::cout <<std::left<<std::setw(4)<< n << std::left<<std::setw(4)<< ell <<":\t"<< std::left<<std::setw(15)<< coeffs_VSH_series[0][n][ell] <<"\t"<< std::left<<std::setw(15)<< coeffs_VSH_series[1][n][ell] <<"\t"<< std::left<<std::setw(15)<< coeffs_VSH_series[2][n][ell] << std::endl;

			coeffs_total_r += coeffs_VSH_series[0][n][ell];
			coeffs_total_1 += coeffs_VSH_series[1][n][ell];
			coeffs_total_2 += coeffs_VSH_series[2][n][ell];
		}
	}
	
	std::cout << "\nCoeffs total:\t" << coeffs_total_r <<"\t"<< coeffs_total_1 <<"\t"<< coeffs_total_2 <<"\t\t"<< coeffs_total_r+coeffs_total_1+coeffs_total_2 << std::endl;
	
	
	//----- Evaluate VSH series at the test point -----
	std::vector<double> A_test_exact { A[0][i_test][j_test][k_test], A[1][i_test][j_test][k_test], A[2][i_test][j_test][k_test] };
	std::vector<double> A_test_calc (3);
	
	//--- n=0 terms ---
	for( int ell=0; ell<=ell_max; ell++ ){
		double factor = sqrt( ( 2.0*ell+1.0 ) / ( 4.0*pi ) );
		A_test_calc[0] += coeffs_VSH_series[0][0][ell] * 0.5 * factor * P0[ell][j_test];
		A_test_calc[1] += coeffs_VSH_series[1][0][ell] * 0.5 * factor * P1[ell][j_test];
		A_test_calc[2] += coeffs_VSH_series[2][0][ell] * 0.5 * factor * P1[ell][j_test];
	}
	//--- n>0 terms ---
	for( int n=1; n<=n_max; n++ ){
		for( int ell=0; ell<=ell_max; ell++ ){
			double factor = sqrt( ( 2.0*ell+1.0 ) / ( 4.0*pi ) );
			A_test_calc[0] += coeffs_VSH_series[0][n][ell] * factor * Tn[n][i_test] * P0[ell][j_test];
			A_test_calc[1] += coeffs_VSH_series[1][n][ell] * factor * Tn[n][i_test] * P1[ell][j_test];
			A_test_calc[2] += coeffs_VSH_series[2][n][ell] * factor * Tn[n][i_test] * P1[ell][j_test];
		}
	}
	
	std::cout << "\nTest fit for (r,theta,phi) = (" << r_list[i_test] <<","<< t_list[j_test] <<","<< p_list[k_test] << "):\n" << std::endl;
	std::cout <<std::left<<std::setw(7)<< "Cmpnt." <<std::left<<std::setw(28)<< "Calc."        <<std::left<<std::setw(28)<< "Exact"         << std::endl;
	std::cout <<std::left<<std::setw(7)<< "r"      <<std::left<<std::setw(28)<< A_test_calc[0] <<std::left<<std::setw(28)<< A_test_exact[0] << std::endl;
	std::cout <<std::left<<std::setw(7)<< "theta"  <<std::left<<std::setw(28)<< A_test_calc[1] <<std::left<<std::setw(28)<< A_test_exact[1] << std::endl;
	std::cout <<std::left<<std::setw(7)<< "phi"    <<std::left<<std::setw(28)<< A_test_calc[2] <<std::left<<std::setw(28)<< A_test_exact[2] << std::endl;
	
	
	//----- Calculate divergence -----
	double r = r_list[i_test];
	double t = t_list[j_test];
	double p = p_list[k_test];
	double div_A_exact     = div_A_function( r, t, p );
	double div_A_numerical = div_spherical ( A_r_function, A_t_function, A_p_function, r, t, p );
	
	double div_A_calc = 0;
		
	//--- V1 ---
	/*for( int ell=0; ell<=ell_max; ell++ ){
		double sum_1 = 0;
		for( int n=1; n<=n_max; n++ ){
			sum_1 += coeffs_VSH_series[0][n][ell] * (double) n * Un[n-1][i_test];
		}
		double sum_2 = 0.5 * coeffs_VSH_series[0][0][ell];
		double sum_3 = 0.5 * coeffs_VSH_series[1][0][ell];
		for( int n=1; n<=n_max; n++ ){
			sum_2 += coeffs_VSH_series[0][n][ell] * Tn[n][i_test];
			sum_3 += coeffs_VSH_series[1][n][ell] * Tn[n][i_test];
		}
		div_A_calc += ( sum_1 + 2.0*pow(r,-1) * sum_2 - (double)(ell*(ell+1)) * pow(r,-1) * sum_3 ) * sqrt( ( 2.0*ell+1.0 ) / ( 4.0*pi ) ) * P0[ell][j_test];
	}*/
	
	//--- V2 (not the shortest but the most readable) ---
	/*for( int ell=0; ell<=ell_max; ell++ ){
		double sum_1 = 0;
		double sum_2 = 0.5 * coeffs_VSH_series[0][0][ell];
		double sum_3 = 0.5 * coeffs_VSH_series[1][0][ell];
		for( int n=1; n<=n_max; n++ ){
			sum_1 += coeffs_VSH_series[0][n][ell] * (double) n * Un[n-1][i_test];
			sum_2 += coeffs_VSH_series[0][n][ell] * Tn[n][i_test];
			sum_3 += coeffs_VSH_series[1][n][ell] * Tn[n][i_test];
		}
		div_A_calc += ( sum_1 + pow(r,-1) * ( 2.0 * sum_2 - (double)(ell*(ell+1)) * sum_3 ) ) * sqrt( ( 2.0*ell+1.0 ) / ( 4.0*pi ) ) * P0[ell][j_test];
	}*/
	
	
	//--- V3 ---
	/*for( int ell=0; ell<=ell_max; ell++ ){
		double term = pow(r,-1) * ( coeffs_VSH_series[0][0][ell] - 0.5 * (double)(ell*(ell+1)) * coeffs_VSH_series[1][0][ell] );
		for( int n=1; n<=n_max; n++ ){
			term += coeffs_VSH_series[0][n][ell] * (double)n * Un[n-1][i_test];
			term += Tn[n][i_test] * pow(r,-1) * ( 2.0 * coeffs_VSH_series[0][n][ell] - (double)(ell*(ell+1)) * coeffs_VSH_series[1][n][ell] );
		}
		div_A_calc += term * sqrt( ( 2.0*ell+1.0 ) / ( 4.0*pi ) ) * P0[ell][j_test];
	}*/
	
	
	//--- V4 (Sum over n separately to give protection against accidentally not accounting for halving the n=0 term) ---
	for( int ell=0; ell<=ell_max; ell++ ){
		double A_r_ell = 0.5 * coeffs_VSH_series[0][0][ell];
		double A_1_ell = 0.5 * coeffs_VSH_series[1][0][ell];
		double A_r_ell_dr = 0;
		for( int n=1; n<=n_max; n++ ){
			double T = cos( (double) n * acos_r_list[i_test] );
			double U = sin( (double) (n) * acos_r_list[i_test] ) / sin( acos_r_list[i_test] );	// U_{n-1}
			A_r_ell    += coeffs_VSH_series[0][n][ell] * T;
			A_1_ell    += coeffs_VSH_series[1][n][ell] * T;
			A_r_ell_dr += coeffs_VSH_series[0][n][ell] * (double) n * U;
		}
		double factor = sqrt( ( 2.0*ell+1.0 ) / ( 4.0*pi ) ) * P0[ell][j_test];
		div_A_calc += factor * ( A_r_ell_dr + 2.0*pow(r,-1) * A_r_ell - (double)(ell*(ell+1))*pow(r,-1) * A_1_ell );
	}
	
	std::cout << "\nDivergence:\t" << div_A_exact <<"\t"<< div_A_numerical <<"\t"<< div_A_calc << std::endl;
	
	
	//----- Calculate curl -----
	std::vector<double> curl_A_exact     = curl_A_function( r, t, p );
	std::vector<double> curl_A_numerical = curl_spherical ( A_r_function, A_t_function, A_p_function, r, t, p );
	std::cout << "\nCurl:" << std::endl;
	print_vector( curl_A_exact    , "Exact\t" );
	print_vector( curl_A_numerical, "Numerical" );
	
	std::vector<double> curl_A_calc (3);
		
	//--- V1 ---
	/*for( int ell=0; ell<=ell_max; ell++ ){
		double sum_r = 0.5 * coeffs_VSH_series[2][0][ell];
		for( int n=1; n<=n_max; n++ ){
			sum_r += coeffs_VSH_series[2][n][ell] * Tn[n][i_test];
		}
		sum_r *= sqrt( ( 2.0*ell+1.0 ) / ( 4.0*pi ) ) * -(double)(ell*(ell+1)) * pow(r,-1) * P0[ell][j_test];
		curl_A_calc[0] += sum_r;
		
		double sum_t_1 = 0;
		double sum_t_2 = 0.5 * coeffs_VSH_series[2][0][ell];
		for( int n=1; n<=n_max; n++ ){
			sum_t_1 += coeffs_VSH_series[2][n][ell] * (double) n * Un[n-1][i_test];
			sum_t_2 += coeffs_VSH_series[2][n][ell] * Tn[n][i_test];
		}
		double sum_t = sum_t_1 + pow(r,-1) * sum_t_2;
		sum_t *= sqrt( ( 2.0*ell+1.0 ) / ( 4.0*pi ) ) * -P1[ell][j_test];
		curl_A_calc[1] += sum_t;
		
		double sum_p_1 = 0.5 * coeffs_VSH_series[1][0][ell];
		double sum_p_2 = 0.5 * coeffs_VSH_series[0][0][ell];
		double sum_p_3 = 0;
		for( int n=1; n<=n_max; n++ ){
			sum_p_1 += coeffs_VSH_series[1][n][ell] * Tn[n][i_test];
			sum_p_2 += coeffs_VSH_series[0][n][ell] * Tn[n][i_test];
			sum_p_3 += coeffs_VSH_series[1][n][ell] * (double) n * Un[n-1][i_test];
		}
		double sum_p = pow(r,-1) * ( sum_p_1 - sum_p_2 ) + sum_p_3;
		sum_p *= sqrt( ( 2.0*ell+1.0 ) / ( 4.0*pi ) ) * P1[ell][j_test];
		curl_A_calc[2] += sum_p;
	}*/
	
	//--- V2 ---
	/*for( int ell=0; ell<=ell_max; ell++ ){
		double factor = sqrt( ( 2.0*ell+1.0 ) / ( 4.0*pi ) );
		
		double sum_r = 0.5 * coeffs_VSH_series[2][0][ell];
		for( int n=1; n<=n_max; n++ ){
			sum_r += coeffs_VSH_series[2][n][ell] * Tn[n][i_test];
		}
		curl_A_calc[0] += sum_r * factor * -(double)(ell*(ell+1)) * pow(r,-1) * P0[ell][j_test];
		
		double sum_t = 0.5 * pow(r,-1) * coeffs_VSH_series[2][0][ell];
		for( int n=1; n<=n_max; n++ ){
			sum_t += coeffs_VSH_series[2][n][ell] * (double) n * Un[n-1][i_test];
			sum_t += pow(r,-1) * coeffs_VSH_series[2][n][ell] * Tn[n][i_test];
		}
		curl_A_calc[1] += sum_t * factor *  -P1[ell][j_test];
		
		double sum_p = 0.5 * pow(r,-1) * ( coeffs_VSH_series[1][0][ell] - coeffs_VSH_series[0][0][ell] );
		for( int n=1; n<=n_max; n++ ){
			sum_p += pow(r,-1) * Tn[n][i_test] * ( coeffs_VSH_series[1][n][ell] - coeffs_VSH_series[0][n][ell] );
			sum_p += coeffs_VSH_series[1][n][ell] * (double) n * Un[n-1][i_test];
		}
		curl_A_calc[2] += sum_p * factor * P1[ell][j_test];
	}*/
	
	//--- V3 ---
	/*for( int ell=0; ell<=ell_max; ell++ ){
		double factor = sqrt( ( 2.0*ell+1.0 ) / ( 4.0*pi ) );
		double sum_r = 0.5 * coeffs_VSH_series[2][0][ell];
		double sum_t = 0.5 * pow(r,-1) * coeffs_VSH_series[2][0][ell];
		double sum_p = 0.5 * pow(r,-1) * ( coeffs_VSH_series[1][0][ell] - coeffs_VSH_series[0][0][ell] );
		for( int n=1; n<=n_max; n++ ){
			double T = Tn[n][i_test];
			double U = (double) n * Un[n-1][i_test];
			sum_r +=   coeffs_VSH_series[2][n][ell] * T;
			sum_t +=   coeffs_VSH_series[2][n][ell] * T * pow(r,-1);
			sum_t +=   coeffs_VSH_series[2][n][ell] * U;
			sum_p +=   coeffs_VSH_series[1][n][ell] * U;
			sum_p += ( coeffs_VSH_series[1][n][ell] - coeffs_VSH_series[0][n][ell] ) * T *  pow(r,-1);
		}
		curl_A_calc[0] += sum_r * factor * -(double)(ell*(ell+1)) * pow(r,-1) * P0[ell][j_test];
		curl_A_calc[1] += sum_t * factor * -P1[ell][j_test];
		curl_A_calc[2] += sum_p * factor *  P1[ell][j_test];
	}*/
	
	//--- V4 ---
	// ell=0 term is always zero. Condense the summation lines. ---
	for( int ell=1; ell<=ell_max; ell++ ){
		double factor = sqrt( ( 2.0*ell+1.0 ) / ( 4.0*pi ) );
		double sum_r = 0.5 * coeffs_VSH_series[2][0][ell];
		double sum_t = 0.5 * pow(r,-1) * coeffs_VSH_series[2][0][ell];
		double sum_p = 0.5 * pow(r,-1) * ( coeffs_VSH_series[1][0][ell] - coeffs_VSH_series[0][0][ell] );
		for( int n=1; n<=n_max; n++ ){
			double T = Tn[n][i_test];
			double U = (double) n * Un[n-1][i_test];
			sum_r += coeffs_VSH_series[2][n][ell] * T;
			sum_t += coeffs_VSH_series[2][n][ell] * T * pow(r,-1) + coeffs_VSH_series[2][n][ell] * U;
			sum_p += ( coeffs_VSH_series[1][n][ell] - coeffs_VSH_series[0][n][ell] ) * T * pow(r,-1) + coeffs_VSH_series[1][n][ell] * U;
		}
		curl_A_calc[0] += sum_r * factor * -(double)(ell*(ell+1)) * pow(r,-1) * P0[ell][j_test];
		curl_A_calc[1] += sum_t * factor * -P1[ell][j_test];
		curl_A_calc[2] += sum_p * factor *  P1[ell][j_test];
	}
	
	print_vector( curl_A_calc, "Calculated" );
	
	
	//----- Output execution time and finish code -----
	auto time_stop = std::chrono::high_resolution_clock::now();
	double exec_time = std::chrono::duration_cast<std::chrono::nanoseconds>( time_stop - time_start ).count() * 1e-9;
	std::cout << "\nExecution time:\t" << floor(exec_time/60) << " m " << exec_time-60*floor(exec_time/60) << " s." << std::endl;
	return 0;
}