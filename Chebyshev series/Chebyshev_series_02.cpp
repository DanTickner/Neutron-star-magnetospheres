/*

cd OneDrive\PhD\Codes\Chebyshev series
g++ Chebyshev_series_02.cpp -o Chebyshev_series_02
Chebyshev_series_02

Read a function f(x) from a header file and calculate its Chebyshev series expansion on [-1,1].

V03: Include calculations of df/dx and d2f/dx2.

*/


#include <iostream>
#include <iomanip>
#include <math.h>
#include <vector>
#include "../Numerical_Integration.h"
#include "../Numerical_Differentiation.h"

const double pi = acos(-1);


//----- Call header file with target function definitions -----
//include "Headers/20230623_chebyshev_header_test.h"
//#include "Headers/20230624_Un_hardcoded.h"
//#include "Headers/20230624_Un.h"
#include "Headers/20230627_Tn.h"



//----- Function declarations -----
double chebyshev_series_integrand( double x, int n );


int main(){
	
	//----- Define variables -----
	int n_max = 7;						// Maximum order for Chebyshev series.
	int n_steps_integration = 1e4;
	int w = 14;							// Fixed width for screen output.
	double x = 0.321;					// Value of x at which to test the series.
	
	std::vector<double> coeffs( n_max+1 );
	
	
	
	
	//----- Fail code if invalid value of x chosen -----
	if( ( x > 1 ) or ( x < -1 ) ){
		std::cout << "Must choose x in [-1,1]. We have x=" << x << ". Stopping the code." << std::endl;
		return 1;
	}
	
	
	
	
	//----- Calculate Chebyshev series coefficients -----
	std::cout << "n" <<"\t"<<std::left<<std::setw(w)<< "coeffs[n]" <<std::left<<std::setw(w)<< "a_n_guess(n)" << std::endl;
	
	for( int n=0; n<=n_max; n++ ){
	
		coeffs[n] = 2.0/pi * rk4( chebyshev_series_integrand, n, 0, pi, n_steps_integration );
		
		if( n == 0 ){
			coeffs[n] *= 0.5;
		}
		
		double coeff_guess = a_n_guess( n );	// Separate to the below in case we need to cout any values during derivation of the analytic coeffs.
		std::cout << n <<"\t"<<std::left<<std::setw(w)<< coeffs[n] <<std::left<<std::setw(w)<< coeff_guess << std::endl;
		
	}
	
	
	
	
	//----- Test the decomposition at the chosen value of x -----
	double f_from_series = 0;
	
	for( int n=0; n<=n_max; n++ ){
		
		double T_n = cos( n * acos( x ) );
		
		f_from_series += coeffs[n] * T_n;
		
	}
	
	double f_from_function = f( x );
	
	
	
	
	//----- Determine Chebyshev series of first derivative -----
	double dfdx_from_numerical = dx( f, x );
	double dfdx_from_function = dfdx( x );
	
	double dfdx_from_series_manual = 0;
	
	double T_0 = 1.0;
	double T_1 = x;
	double T_2 = 2.0*x*x - 1.0;
	double T_3 = 4.0*pow(x,3) - 3.0*x;
	double term_0 = 1.0 * T_0 * ( 1.0*coeffs[1] + 3.0*coeffs[3] + 5.0*coeffs[5] );
	double term_1 = 2.0 * T_1 * ( 2.0*coeffs[2] + 4.0*coeffs[4] + 6.0*coeffs[6] );
	double term_2 = 2.0 * T_2 * ( 3.0*coeffs[3] + 5.0*coeffs[5] );
	double term_3 = 2.0 * T_3 * ( 4.0*coeffs[4] + 6.0*coeffs[6] );
	
	std::cout << std::endl;
	std::cout << 0 << "\t" << term_0 << std::endl;
	std::cout << 1 << "\t" << term_1 << std::endl;
	std::cout << 2 << "\t" << term_2 << std::endl;
	std::cout << 3 << "\t" << term_3 << std::endl;
	
	
	dfdx_from_series_manual += term_0 + term_1 + term_2 + term_3;
	
	
	
	
	
	
	
	
	
	
	double dfdx_from_series = 0;
	
	/*
	//--- V1 ---
	for( int n=0; n<=n_max; n++ ){
		for( int k=n; k<=n_max; k++ ){
			
			if( (k+n)%2 == 1 ){
			
				double T_n = cos( n * acos( x ) );
				
				double term = (double) 2.0 * k * coeffs[k] * T_n;
				
				if( n == 0 ){
					term *= 0.5;
				}
				
				dfdx_from_series += term;
				
			}
			
		}
	}
	*/
	
	
	/*
	//--- V2 ---
	for( int n=0; n<=n_max; n++ ){
		
		double T_n = cos( n * acos( x ) );
		
		double term = 0;
		
		for( int k=n; k<=n_max; k++ ){			
			if( (k+n)%2 == 1 ){
				term += k * coeffs[k];
			}
		}
			
		term *= 2.0 * T_n;
		
		if( n == 0 ){
			term *= 0.5;
		}
		
		dfdx_from_series += term;
		
	}
	*/
	
	
	
	/*
	//--- V3 ---
	for( int n=0; n<=n_max; n++ ){
		
		double T_n = cos( n * acos( x ) );
		
		double term = 0;
		
		for( int k=n+1; k<=n_max; k+=2 ){
			term += k * coeffs[k];
		}
			
		term *= 2.0 * T_n;
		
		if( n == 0 ){
			term *= 0.5;
		}
		
		dfdx_from_series += term;
		
	}
	*/
	
	
	//--- V4 ---
	for( int k=1; k<=n_max; k+=2 ){
		dfdx_from_series += k * coeffs[k];
	}
		
	for( int n=1; n<=n_max; n++ ){
		
		double T_n = cos( n * acos( x ) );
		
		double term = 0;
		
		for( int k=n+1; k<=n_max; k+=2 ){
			term += k * coeffs[k];
		}
			
		term *= 2.0 * T_n;
		
		dfdx_from_series += term;
		
	}
	
	
	
	
	//----- Determine Chebyshev series of second derivative -----
	double d2fdx2_from_function  = d2fdx2( x );
	double d2fdx2_from_numerical = dx2( f, x );
	
	
	double d2fdx2_from_series = 0;
	
	
	/*
	//--- V1 ---
	int infinity = 3.0*n_max;
	
	std::cout << "\nd2f/dx2" << std::endl;
	std::cout << "n\tk\tell\td2f/dx2 so far" << std::endl;
	
	for( int n=0; n<=infinity; n++ ){
		for( int k=0; k<=infinity; k++ ){
			for( int ell=0; ell<=infinity; ell++ ){
				
				bool g1 = ( n <= ell-1 ) and ( (n+ell-1)%2 == 0 );
				bool g2 = ( ell <= k-1 ) and ( (ell+k-1)%2 == 0 );
				bool k_escape = k < coeffs.size();	// Otherwise list overflow and nans.
				
				if( g1 and g2 and k_escape ){
					
					double h_n = 1;
					if( n == 0 ){
						h_n = 0.5;
					}
					
					double T_n = cos( n * acos( x ) );
					
					d2fdx2_from_series += 4.0 * coeffs[k] * k * ell * h_n * T_n;
					std::cout << n <<"\t"<< k <<"\t"<< ell <<"\t"<< d2fdx2_from_series << std::endl;
					
				}
				
			}
		}
	}
	
	std::cout << std::endl;
	*/
	
	
	/*
	//--- V2 ---
	for( int n=0; n<=n_max-2; n++ ){
		for( int k=n+2; k<=n_max; k+=2 ){
			for( int ell=n+1; ell<=k-1; ell+=2 ){
				
				double h_n = 1;
				if( n == 0 ){
					h_n = 0.5;
				}
				
				double T_n = cos( n * acos( x ) );
				
				d2fdx2_from_series += 4.0 * coeffs[k] * k * ell * h_n * T_n;
				
			}
		}
	}
	*/
	
	
	/*
	//--- V3 (separate the n=0 to remove the h_k. Not a major simplification.) ---
	for( int k=2; k<=n_max; k+=2 ){
		for( int ell=1; ell<=k-1; ell+=2 ){
			d2fdx2_from_series += 2.0 * k * ell * coeffs[k];
		}
	}
	
	for( int n=1; n<=n_max-2; n++ ){
		for( int k=n+2; k<=n_max; k+=2 ){
			for( int ell=n+1; ell<=k-1; ell+=2 ){
				double T_n = cos( n * acos( x ) );
				d2fdx2_from_series += 4.0 * k * ell * coeffs[k] * T_n;
			}
		}
	}
	*/
	
	
	/*
	//--- V4 (evaluate the sum over q) ---
	for( int n=0; n<=n_max; n++ ){
		for( int k=n+2; k<=n_max; k++ ){
			if( (n+k)%2==0 ){
				double h_n = 1.0 - 0.5 * ( n == 0 );
				double T_n = cos( n * acos( x ) );
				d2fdx2_from_series += h_n * k * ( k*k - n*n ) * coeffs[k] * T_n;
			}
		}
	}
	*/
	
	
	/*
	//--- V5 ---
	for( int n=0; n<=n_max; n++ ){
		for( int k=n+2; k<=n_max; k+=2 ){
			double h_n = 1.0 - 0.5 * ( n == 0 );
			double T_n = cos( n * acos( x ) );
			d2fdx2_from_series += h_n * k * ( k*k - n*n ) * coeffs[k] * T_n;
		}
	}
	*/
	
	
	//--- V6 (separate n=0 term, no major simplification) ---
	for( int k=2; k<=n_max; k+=2 ){
		d2fdx2_from_series += 0.5 * pow(k,3) * coeffs[k];
	}
	
	for( int n=1; n<=n_max; n++ ){
		for( int k=n+2; k<=n_max; k+=2 ){
			double T_n = cos( n * acos( x ) );
			d2fdx2_from_series += k * ( k*k - n*n ) * coeffs[k] * T_n;
		}
	}
	
	
			
	
	
	
	
	//----- Output results -----
	std::cout << std::endl;
	std::cout << "f(x) from function definition:\t" << f_from_function << std::endl;
	std::cout << "f(x) from Chebyshev series   :\t" << f_from_series   << std::endl;
	
	
	std::cout << std::endl;
	std::cout << "df/dx from function definition :\t" << dfdx_from_function  << std::endl;
	std::cout << "df/dx from numerical derivative:\t" << dfdx_from_numerical << std::endl;
	std::cout << "df/dx from Chebyshev series manual:\t" << dfdx_from_series_manual    << std::endl;
	std::cout << "df/dx from Chebyshev series    :\t" << dfdx_from_series    << std::endl;
	
	
	std::cout << std::endl;
	std::cout << "d2f/dx2 from function definition :\t" << d2fdx2_from_function  << std::endl;
	std::cout << "d2f/dx2 from numerical derivative:\t" << d2fdx2_from_numerical << std::endl;
	std::cout << "d2f/dx2 from Chebyshev series    :\t" << d2fdx2_from_series    << std::endl;
		
	std::cout << "\nDone" << std::endl;
	
	return 0;
}




double chebyshev_series_integrand( double x, int n ){
	return f( cos(x) ) * cos( n * x );
}