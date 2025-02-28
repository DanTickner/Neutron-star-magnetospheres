/*

cd OneDrive\PhD\Codes\Chebyshev series
g++ Chebyshev_series_f_times_g_01.cpp -o Chebyshev_series_f_times_g_01 --std=c++11
Chebyshev_series_f_times_g_01

Read two functions f(x) and g(x) from a header file and calculate the Chebyshev series expansion of f(x)*g(x) on [-1,1]
in terms of their individual Chebyshev series.


*/


#include <iostream>
#include <iomanip>
#include <math.h>
#include <vector>
#include "../Numerical_Integration.h"

const double pi = acos(-1);


//----- Call header file with target function definitions -----
#include "Headers/20230628_f_times_g_polynomials.h"



//----- Function declarations -----
double chebyshev_series_f_integrand ( double x, int n );
double chebyshev_series_g_integrand ( double x, int n );
double chebyshev_series_fg_integrand( double x, int n );
double delta( int a, int b );
double T_n( double x, int n );

int main(){
	
	//----- Define variables -----
	int n_max = 4;						// Maximum order for Chebyshev series.
	int n_steps_integration = 1e4;
	int w = 14;							// Fixed width for screen output.
	double x = 0.321;					// Value of x at which to test the series.
	
	std::vector<double> coeffs_f ( n_max+1 );
	std::vector<double> coeffs_g ( n_max+1 );
	std::vector<double> coeffs_fg( 2*n_max+1 );
	
	
	
	
	//----- Fail code if invalid value of x chosen -----
	if( ( x > 1 ) or ( x < -1 ) ){
		std::cout << "Must choose x in [-1,1]. We have x=" << x << ". Stopping the code." << std::endl;
		return 1;
	}
	
	
	
	
	//----- Calculate Chebyshev series coefficients -----
	
	//--- f(x) ---
	std::cout << "n" <<"\t"<<std::left<<std::setw(w)<< "coeffs[n]" <<std::left<<std::setw(w)<< "a_n_guess(n)" << std::endl;
	
	for( int n=0; n<=n_max; n++ ){
		coeffs_f[n] = 2.0/pi * rk4( chebyshev_series_f_integrand, n, 0, pi, n_steps_integration ) * (1-0.5*(n==0));
		double coeff_guess = f_n_guess( n );
		std::cout << n <<"\t"<<std::left<<std::setw(w)<< coeffs_f[n] <<std::left<<std::setw(w)<< coeff_guess << std::endl;
	}
	
	//--- g(x) ---
	std::cout << "\nn" <<"\t"<<std::left<<std::setw(w)<< "coeffs[n]" <<std::left<<std::setw(w)<< "a_n_guess(n)" << std::endl;
	
	for( int n=0; n<=n_max; n++ ){
		coeffs_g[n] = 2.0/pi * rk4( chebyshev_series_g_integrand, n, 0, pi, n_steps_integration ) * (1-0.5*(n==0));
		double coeff_guess = g_n_guess( n );
		std::cout << n <<"\t"<<std::left<<std::setw(w)<< coeffs_g[n] <<std::left<<std::setw(w)<< coeff_guess << std::endl;
	}
	
	//--- f(x)*g(x) ---
	std::cout << "\nn" <<"\t"<<std::left<<std::setw(w)<< "coeffs[n]" <<std::left<<std::setw(w)<< "a_n_guess(n)" << std::endl;
	
	for( int n=0; n<=2*n_max; n++ ){
		
		double coeff_guess = fg_n_guess( n );
		
		
		//---V0 (numerical) ---
		//coeffs_fg[n] = 2.0/pi * rk4( chebyshev_series_fg_integrand, n, 0, pi, n_steps_integration ) * (1-0.5*(n==0));
		
		/*
		//--- V1 (incorrect)---
		std::cout << std::endl;
		for( int q=0; q<=n; q++ ){
			coeffs_fg[n] += 0.5 * coeffs_f[q] * coeffs_g[n-q];
			std::cout << "s1\t" << n <<"\t"<< q <<"\t"<< n-q <<"\t"<< coeffs_f[q] <<"\t"<< coeffs_g[n-q] <<"\t"<< std::endl;
		}
		
		for( int q=n; q<=n_max; q++ ){
			coeffs_fg[n] += 0.5 * coeffs_f[q] * coeffs_g[q-n];
			std::cout << "s2\t" << n <<"\t"<< q <<"\t"<< q-n <<"\t"<< coeffs_f[q] <<"\t"<< coeffs_g[q-n] <<"\t"<< std::endl;
		}
		std::cout << std::endl;
		*/
		
		
		/*
		//--- V2 ---
		if( n == 0 ){
			coeffs_fg[n] = coeffs_f[0] * coeffs_g[0];
			for( int i=1; i<=n_max; i++ ){
				coeffs_fg[n] += 0.5 * coeffs_f[i] * coeffs_g[i];
			}
		}
		else if( n == 2*n_max ){
			coeffs_fg[n] = 0.5 * coeffs_f[n_max] * coeffs_g[n_max];
		}
		*/
		
		/*
		//--- V3 ---
		coeffs_fg[n] += coeffs_f[0] * coeffs_g[n] + coeffs_f[n] * coeffs_g[0];
		
		for( int i=1; i<n; i++ ){
			coeffs_fg[n] += 0.5 * coeffs_f[i] * coeffs_g[ abs( i-n ) ];
			coeffs_fg[n] += 0.5 * coeffs_f[ abs( i-n ) ] * coeffs_g[i];
		}
		*/
		
		
		//---V4 (incorrect) ---
		/*
		if( n == 0 ){
			coeffs_fg[n] = coeffs_f[0] * coeffs_g[0];
			for( int i=1; i<=n_max; i++ ){
				coeffs_fg[n] += 0.5 * coeffs_f[i] * coeffs_g[i];
			}
		}
		
		else {
			for( int n1=0; n1<=2*n_max; n1++ ){
				if( ( n-n1 >= 0 ) and ( n-n1 <= 2*n_max ) ){
					coeffs_fg[n] += 0.5 * coeffs_f[n1] * coeffs_g[n-n1] + coeffs_f[n-n1] * coeffs_g[n1];
				}
				if( ( n+n1 >= 0 ) and ( n+n1 <= 2*n_max ) ){
					coeffs_fg[n] += 0.5 * coeffs_f[n1] * coeffs_g[n+n1] + coeffs_f[n+n1] * coeffs_g[n1];
				}
			}
			if( n%2 == 0 ){
				coeffs_fg[n] -= coeffs_f[n/2] * coeffs_g[n/2];
			}
			//coeffs_fg[n] *= 1.0/pi;
		}
		*/
		
		
		
		
		
		/*
		//--- V5 (numerical integrals to see where error lies in my maths) ---
		if( n == 0 ){
			coeffs_fg[n] = coeffs_f[0] * coeffs_g[0];
			for( int i=1; i<=n_max; i++ ){
				coeffs_fg[n] += 0.5 * coeffs_f[i] * coeffs_g[i];
			}
		}
		
		else {
			for( int n1=0; n1<=n_max; n1++ ){
				for( int n2=n1; n2<=n_max; n2++ ){
					coeffs_fg[n] += 0.5 * ( coeffs_f[n1]*coeffs_g[n2] + coeffs_f[n2]*coeffs_g[n1] ) * delta( n1+n2, n );
					coeffs_fg[n] += 0.5 * ( coeffs_f[n1]*coeffs_g[n2] + coeffs_f[n2]*coeffs_g[n1] ) * delta( n2-n1, n );
					//std::cout << n<<"\t"<< n1 <<"\t"<< n2 <<"\t"<< delta( n1+n2, n ) <<"\t"<< delta( n2-n1, n ) << std::endl;
					if(  delta( n1+n2, n ) ){
						std::cout << n <<"\t"<< n1 <<"\t"<< n2 << "\tdelta(n1+n2,n)=1" << std::endl;
					}
				}
				coeffs_fg[n] -= 0.5 * coeffs_f[n1]*coeffs_g[n1] * delta( 2*n1, n );
			}
		}
		*/
		
		
		
		//--- V6 ---
		
		if( n == 0 ){
			coeffs_fg[n] = coeffs_f[0] * coeffs_g[0];
			for( int i=1; i<=n_max; i++ ){
				coeffs_fg[n] += 0.5 * coeffs_f[i] * coeffs_g[i];
			}
		}
		
		else {
			for( int n1=0; n1<=n_max; n1++ ){
				if( not( n-n1 <= n_max ) ){
					std::cout << n <<"\t"<< n1 <<"\t"<< n-n1 <<"\t"<< n_max <<"\t!(n-n1<=n_max)" << std::endl;
				}
				if( ( n-n1 >= n1 ) and ( n-n1 <= n_max ) ){
					coeffs_fg[n] += 0.5 * ( coeffs_f[n1]*coeffs_g[n-n1] + coeffs_f[n-n1]*coeffs_g[n1] ) * delta( n1+n-n1, n );
				}
				if( ( n+n1 >= n1 ) and ( n+n1 <= n_max ) ){
					coeffs_fg[n] += 0.5 * ( coeffs_f[n1]*coeffs_g[n+n1] + coeffs_f[n+n1]*coeffs_g[n1] ) * delta( n+n1-n1, n );
				}
				coeffs_fg[n] -= 0.5 * coeffs_f[n1]*coeffs_g[n1] * delta( 2*n1, n );
			}
		}
		
			
		
				
		
		
		
		std::cout << n <<"\t"<<std::left<<std::setw(w)<< coeffs_fg[n] <<std::left<<std::setw(w)<< coeff_guess << std::endl;
	}
	
	
	
	
	//----- Test the decomposition at the chosen value of x -----
	double f_from_series  = 0;
	double g_from_series  = 0;
	double fg_from_series = 0;
	
	for( int n=0; n<=n_max; n++ ){
		f_from_series  += coeffs_f [n] * cos( n * acos( x ) );
		g_from_series  += coeffs_g [n] * cos( n * acos( x ) );
	}
	for( int n=0; n<=2*n_max; n++ ){
		fg_from_series += coeffs_fg[n] * cos( n * acos( x ) );
	}
	
	double f_from_function  = f( x );
	double g_from_function  = g( x );
	double fg_from_function = fg( x );
	double fg_from_function_explicit = fg_explicit( x );
	
	double f_from_exact_coeffs  = 0;
	double g_from_exact_coeffs  = 0;
	double fg_from_exact_coeffs = 0;
	
	for( int n=0; n<=n_max; n++ ){
		f_from_exact_coeffs  += f_n_guess ( n ) * T_n( x, n );
		g_from_exact_coeffs  += g_n_guess ( n ) * T_n( x, n );
	}
	for( int n=0; n<=2*n_max; n++ ){
		fg_from_exact_coeffs += fg_n_guess( n ) * T_n( x, n );
	}
	
	
	//----- Guess fg by multiplying the series -----
	double fg_from_multiplying = 0;
	
	//-- V1 --- 
	/*
	for( int n=0; n<=n_max; n++ ){
		for( int m=n; m<=n_max; m++ ){
			fg_from_multiplying += ( coeffs_f[n]*coeffs_g[m] + coeffs_f[m]*coeffs_g[n] ) * T_n(x,n) * T_n(x,m);
		}
		fg_from_multiplying -= coeffs_f[n] * coeffs_g[n] * pow( T_n(x,n), 2 );	// Avoid double-counting.
	}
	*/
	
	
	//--- V2 ---
	/*
	for( int n=0; n<=n_max; n++ ){
		for( int m=n; m<=n_max; m++ ){
			fg_from_multiplying += ( coeffs_f[n]*coeffs_g[m] + coeffs_f[m]*coeffs_g[n] ) * 0.5 * ( T_n(x,n+m) + T_n(x,n-m) );
		}
		fg_from_multiplying -= coeffs_f[n] * coeffs_g[n] * 0.5 * ( T_n(x,2*n) + T_n(x,0) );
	}
	*/
	
	//--- V3 ---
	for( int n=0; n<=n_max; n++ ){
		for( int m=n; m<=n_max; m++ ){
			fg_from_multiplying += ( coeffs_f[n]*coeffs_g[m] + coeffs_f[m]*coeffs_g[n] ) * ( T_n(x,n+m) + T_n(x,m-n) );
		}
		fg_from_multiplying -= coeffs_f[n] * coeffs_g[n] * ( T_n(x,2*n) + T_n(x,0) );
	}
	fg_from_multiplying *= 0.5;
	
	
	//----- Output results -----
	std::cout << std::endl;
	std::cout << "f(x) from Chebyshev series   :\t" << f_from_series       << std::endl;
	std::cout << "f(x) from function definition:\t" << f_from_function     << std::endl;
	std::cout << "f(x) from exact coefficients :\t" << f_from_exact_coeffs << std::endl;
	
	std::cout << std::endl;
	std::cout << "g(x) from Chebyshev series   :\t" << g_from_series       << std::endl;
	std::cout << "g(x) from function definition:\t" << g_from_function     << std::endl;
	std::cout << "g(x) from exact coefficients :\t" << g_from_exact_coeffs << std::endl;
	
	std::cout << std::endl;
	std::cout << "f(x)*g(x) from Chebyshev series   :\t" << fg_from_series            << std::endl;
	std::cout << "f(x)*g(x) from function definition:\t" << fg_from_function          << std::endl;
	std::cout << "f(x)*g(x) from function explicit  :\t" << fg_from_function_explicit << std::endl;
	std::cout << "f(x)*g(x) from exact coefficients :\t" << fg_from_exact_coeffs      << std::endl;
	std::cout << "f(x)*g(x) from multiplying series :\t" << fg_from_multiplying       << std::endl;
		
	std::cout << "\nDone" << std::endl;
	
	return 0;
}




double chebyshev_series_f_integrand( double x, int n ){
	return f( cos(x) ) * cos( n * x );
}

double chebyshev_series_g_integrand( double x, int n ){
	return g( cos(x) ) * cos( n * x );
}

double chebyshev_series_fg_integrand( double x, int n ){
	return f( cos(x) ) * g( cos(x) ) * cos( n * x );
}

double chebyshev_series_Tn_Tm_integrand( double x, int n, int m ){
	return cos( n * x ) * cos( m * x );
}

double delta( int a, int b ){
	return ( a == b );
}

double T_n( double x, int n ){
	return cos( n * acos( x ) );
}