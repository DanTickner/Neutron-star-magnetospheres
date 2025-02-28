/*

cd OneDrive\PhD\Codes\Chebyshev series
g++ Chebyshev_series_x_times_fx_01.cpp -o Chebyshev_series_x_times_fx_01
Chebyshev_series_x_times_fx_01

Read a function f(x) from a header file and calculate the Chebyshev series expansion of x*f(x)on [-1,1].
Compare to derived expression in terms of the coefficients of the Chebyshev series of f(x).

*/


#include <iostream>
#include <iomanip>
#include <math.h>
#include <vector>
#include "../Numerical_Integration.h"

const double pi = acos(-1);


//----- Call header file with target function definitions -----
#include "Headers/20230623_chebyshev_header_test.h"



//----- Function declarations -----
double chebyshev_series_integrand( double x, int n );
double chebyshev_series_xf_integrand( double x, int n );
double delta( int a, int b );


int main(){
	
	//----- Define variables -----
	int n_max = 7;						// Maximum order for Chebyshev series.
	int n_steps_integration = 1e4;
	int w = 14;							// Fixed width for screen output.
	double x = 0.321;					// Value of x at which to test the series.
	
	std::vector<double> coeffs   ( n_max+1 );
	std::vector<double> coeffs_xf( n_max+1 );
	
	
	
	
	//----- Fail code if invalid value of x chosen -----
	if( ( x > 1 ) or ( x < -1 ) ){
		std::cout << "Must choose x in [-1,1]. We have x=" << x << ". Stopping the code." << std::endl;
		return 1;
	}
	
	
	
	
	//----- Calculate Chebyshev series coefficients -----
	
	//--- f(x) ---
	std::cout << "n" <<"\t"<<std::left<<std::setw(w)<< "coeffs[n]" <<std::left<<std::setw(w)<< "a_n_guess(n)" << std::endl;
	
	for( int n=0; n<=n_max; n++ ){
		coeffs[n] = 2.0/pi * rk4( chebyshev_series_integrand, n, 0, pi, n_steps_integration ) * (1-0.5*(n==0));
		double coeff_guess = a_n_guess( n );
		std::cout << n <<"\t"<<std::left<<std::setw(w)<< coeffs[n] <<std::left<<std::setw(w)<< coeff_guess << std::endl;
	}
	
	//--- x f(x) ---
	std::cout << "\nn" <<"\t"<<std::left<<std::setw(w)<< "coeffs[n]" <<std::left<<std::setw(w)<< "a_n_guess(n)" << std::endl;
	
	for( int n=0; n<=n_max; n++ ){
		coeffs_xf[n] = 2.0/pi * rk4( chebyshev_series_xf_integrand, n, 0, pi, n_steps_integration ) * (1-0.5*(n==0));
		
		/*
		//-- coeff_guess v1 ---
		double coeff_guess = 0;
		if( n == 0 ){
			coeff_guess = 0.5 * coeffs[1];
		}
		else if( n == 1 ){
			coeff_guess = coeffs[0] + 0.5 * coeffs[2];
		}
		else if( n < n_max ){
			coeff_guess = 0.5 * ( coeffs[n-1] + coeffs[n+1] );
		}
		else{
			coeff_guess = 0.5 * coeffs[n_max-1];
		}
		*/
		
		//--- coeff_guess v2 ---
		double coeff_guess = 0.5 * ( ( 1+delta(n,1) ) * coeffs[n-1] + coeffs[n+1] );
		
		std::cout << n <<"\t"<<std::left<<std::setw(w)<< coeffs_xf[n] <<std::left<<std::setw(w)<< coeff_guess << std::endl;
	}
	
	std::cout << coeffs[-1] << std::endl;
	
	
	
	
	//----- Test the decomposition at the chosen value of x -----
	double f_from_series = 0;
	
	for( int n=0; n<=n_max; n++ ){
		f_from_series += coeffs[n] * cos( n * acos( x ) );
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

double chebyshev_series_xf_integrand( double x, int n ){
	return cos(x) * f( cos(x) ) * cos( n * x );
}

double delta( int a, int b ){
	return ( a == b );
}