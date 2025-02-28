/*

cd OneDrive\PhD\Codes\Chebyshev series
g++ Chebyshev_series_03.cpp -o Chebyshev_series_03
Chebyshev_series_03

Read a function f(x) from a header file and calculate its Chebyshev series expansion on [-1,1].

V03: Remove intermediate versions of functions.

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
//#include "Headers/20230627_Tn.h"
#include "Headers/20230628_one_over_x.h"



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
	double f_from_series_guess = 0;
	
	for( int n=0; n<=n_max; n++ ){
		//if( n%2 == 0 ){ continue; }
		
		double T_n = cos( n * acos( x ) );
		
		f_from_series += coeffs[n] * T_n;
		f_from_series_guess += a_n_guess(n) * T_n;
		
	}
	
	double f_from_function = f( x );
	
	
	
	
	//----- Determine Chebyshev series of first derivative -----
	double dfdx_from_numerical = dx( f, x );
	double dfdx_from_function = dfdx( x );
	
	double dfdx_from_series = 0;
	
	
	for( int k=1; k<=n_max; k+=2 ){
		dfdx_from_series += k * coeffs[k];
	}
		
	for( int n=1; n<=n_max; n++ ){
		double T_n = cos( n * acos( x ) );
		for( int k=n+1; k<=n_max; k+=2 ){
			dfdx_from_series += 2.0 * k * coeffs[k] * T_n;
		}
	}
	
	
	
	
	//----- Determine Chebyshev series of second derivative -----
	double d2fdx2_from_function  = d2fdx2( x );
	double d2fdx2_from_numerical = dx2( f, x );
	
	double d2fdx2_from_series = 0;
	
	for( int k=2; k<=n_max; k+=2 ){
		for( int ell=1; ell<=k-1; ell+=2 ){
			d2fdx2_from_series += 2.0 * k * ell * coeffs[k];
		}
	}
	
	for( int n=1; n<=n_max-2; n++ ){
		double T_n = cos( n * acos( x ) );
		for( int k=n+2; k<=n_max; k+=2 ){
			for( int ell=n+1; ell<=k-1; ell+=2 ){
				d2fdx2_from_series += 4.0 * k * ell * coeffs[k] * T_n;
			}
		}
	}
	
	
			
	
	
	
	
	//----- Output results -----
	std::cout << std::endl;
	std::cout << "f(x) from function definition   :\t" << f_from_function     << std::endl;
	std::cout << "f(x) from Chebyshev series guess:\t" << f_from_series_guess << std::endl;
	std::cout << "f(x) from Chebyshev series      :\t" << f_from_series       << std::endl;
	
	
	std::cout << std::endl;
	std::cout << "df/dx from function definition :\t" << dfdx_from_function  << std::endl;
	std::cout << "df/dx from numerical derivative:\t" << dfdx_from_numerical << std::endl;
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