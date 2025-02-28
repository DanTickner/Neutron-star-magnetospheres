/*

cd OneDrive\PhD\Codes\Chebyshev series
g++ Chebyshev_series_01.cpp -o Chebyshev_series_01
Chebyshev_series_01

Read a function f(x) from a header file and calculate its Chebyshev series expansion on [-1,1].

*/


#include <iostream>
#include <iomanip>
#include <math.h>
#include <vector>
#include "../Numerical_Integration.h"

const double pi = acos(-1);


//----- Call header file with target function definitions -----
//#include "Headers/20230623_chebyshev_header_test.h"
//#include "Headers/20230624_Un_hardcoded.h"
#include "Headers/20230624_Un.h"



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
	
	
	
	
	//----- Output results -----
	std::cout << std::endl;
	std::cout << "f(x) from Chebyshev series   :\t" << f_from_series   << std::endl;
	std::cout << "f(x) from function definition:\t" << f_from_function << std::endl;
		
	std::cout << "\nDone" << std::endl;
	
	return 0;
}




double chebyshev_series_integrand( double x, int n ){
	return f( cos(x) ) * cos( n * x );
}