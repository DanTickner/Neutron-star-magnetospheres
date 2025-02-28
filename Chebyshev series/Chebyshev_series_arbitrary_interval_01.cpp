/*

cd OneDrive\PhD\Codes\Chebyshev series
g++ Chebyshev_series_arbitrary_interval_01.cpp -o Chebyshev_series_arbitrary_interval_01
Chebyshev_series_arbitrary_interval_01

Read a function f(x) from a header file and calculate its Chebyshev series expansion on [x_min,x_max], where x_min and x_max can be freely chosen.

V03: Based on Chebyshev_series_03.cpp.

*/


#include <iostream>
#include <iomanip>
#include <math.h>
#include <vector>
#include "../Numerical_Integration.h"
#include "../Numerical_Differentiation.h"

const double pi = acos(-1);


//----- Call header file with target function definitions -----
#include "Headers/B_r_1.h"



//----- Function declarations -----
double Lambda( double x, double x_min, double x_max );
double Lambda_inverse( double x, double x_min, double x_max );
double T_n( double x, int n, double x_min, double x_max );
double chebyshev_series_integrand( double x, int n, double x_min, double x_max );


int main(){
	
	//----- Define variables -----
	int n_max = 20;						// Maximum order for Chebyshev series.
	int n_steps_integration = 100;
	int w = 14;							// Fixed width for screen output.
	int x_min = 1.0;					// Minimum and maxmium values of x over which the fit should be valid.
	int x_max = 80.0;
	double x = 2.321;					// Value of x at which to test the series.
	
	std::vector<double> coeffs( n_max+1 );
	
	
	
	
	//----- Fail code if invalid value of x chosen -----
	if( ( x > x_max ) or ( x < x_min ) ){
		std::cout << "Must choose x in [" << x_min << "," << x_max << "]. We have x=" << x << ". Stopping the code." << std::endl;
		return 1;
	}
	
	
	
	
	//----- Calculate Chebyshev series coefficients -----
	std::cout << "n" <<"\t"<<std::left<<std::setw(w)<< "coeffs[n]" <<std::left<<std::setw(w)<< "a_n_guess(n)" << std::endl;
	
	for( int n=0; n<=n_max; n++ ){
	
		coeffs[n] = 2.0/pi * rk4( chebyshev_series_integrand, n, x_min, x_max, 0, pi, n_steps_integration );
		
		if( n == 0 ){
			coeffs[n] *= 0.5;
		}
		
		double coeff_guess = a_n_guess( n );	// Separate to the below in case we need to cout any values during derivation of the analytic coeffs.
		std::cout << n <<"\t"<<std::left<<std::setw(w)<< coeffs[n] <<std::left<<std::setw(w)<< coeff_guess << std::endl;
		
	}
	
	
	
	
	//----- Test the decomposition at the chosen value of x -----
	double f_from_series = 0;
	
	for( int n=0; n<=n_max; n++ ){
		f_from_series += coeffs[n] * T_n( x, n, x_min, x_max );
	}
	
	double f_from_function = f( x );
	
	
	
	
	//----- Determine Chebyshev series of first derivative -----
	double dfdx_from_numerical = dx( f, x );
	double dfdx_from_function = dfdx( x );
	
	double dfdx_from_series = 0;
	
	//--- V1 (analytic expression) ---
	for( int n=0; n<=n_max; n++ ){
		for( int k=n+1; k<=n_max; k++ ){
			if( (n+k)%2 == 1 ){
				double h_n = 1.0 - 0.5 * ( n == 0 );
				dfdx_from_series += 4.0/(x_max-x_min) * h_n * k * coeffs[k] * T_n( x, n, x_min, x_max );
			}
		}
	}
	
	
	/*
	//--- V2 (encoded expression) ---
	for( int k=1; k<=n_max; k+=2 ){
		dfdx_from_series += k * coeffs[k];
	}
		
	for( int n=1; n<=n_max; n++ ){
		for( int k=n+1; k<=n_max; k+=2 ){
			dfdx_from_series += 2.0 * k * coeffs[k] * T_n( x, n, x_min, x_max );
		}
	}
	
	dfdx_from_series *= 2.0/(x_max-x_min);
	*/
	
	
	
	
	//----- Determine Chebyshev series of second derivative -----
	double d2fdx2_from_function  = d2fdx2( x );
	double d2fdx2_from_numerical = dx2( f, x );
	
	double d2fdx2_from_series = 0;
	
	/*
	//--- V1 (analytic expression) ---
	for( int n=0; n<=n_max; n++ ){
		
		double term_k = 0;
		for( int k=n+2; k<=n_max; k+=2 ){
			if( (n+k)%2==0 ){
				term_k += k * ( k*k - n*n ) * coeffs[k];
			}
		}
		
		double h_n = 1.0 - 0.5 * ( n == 0 );
		d2fdx2_from_series += 4.0*pow(x_max-x_min,-2) * h_n * term_k * T_n( x, n, x_min, x_max );
	}
	*/
	
	
	
	//--- V2 (encoded expression) ---
	for( int k=2; k<=n_max; k+=2 ){
		d2fdx2_from_series += 0.5 * pow(k,3) * coeffs[k];
	}
	
	for( int n=1; n<=n_max; n++ ){
		for( int k=n+2; k<=n_max; k+=2 ){
			d2fdx2_from_series += k * ( k*k - n*n ) * coeffs[k] * T_n( x, n, x_min, x_max );
		}
	}
	
	d2fdx2_from_series *= 4.0 * pow( x_max-x_min, -2 );
	
	
	
			
	
	
	
	
	//----- Output results -----
	std::cout << std::endl;
	std::cout << "f(x) from function definition:\t" << f_from_function << std::endl;
	std::cout << "f(x) from Chebyshev series   :\t" << f_from_series   << std::endl;
	
	
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




double Lambda( double x, double x_min, double x_max ){
	// Convert x in [-1,1] to x in [x_min,x_max].
	return 0.5 * ( (x_max-x_min)*x + x_min+x_max );
}

double Lambda_inverse( double x, double x_min, double x_max ){
	// Convert x in [x_min,x_max] to x in [-1,1].
	return pow( x_max-x_min, -1 ) * ( 2.0*x - (x_max+x_min) );
}

double T_n( double x, int n, double x_min, double x_max ){
	return cos( n * acos( Lambda_inverse( x, x_min, x_max ) ) );
}

double chebyshev_series_integrand( double x, int n, double x_min, double x_max ){
	return f( Lambda( cos(x), x_min, x_max ) ) * cos( n * x );
}