/*

cd OneDrive\PhD\Codes\Chebyshev series
g++ Chebyshev_series_f_times_g_04.cpp -o Chebyshev_series_f_times_g_04 --std=c++11
Chebyshev_series_f_times_g_04

Read two functions f(x) and g(x) from a header file and calculate the Chebyshev series expansion of f(x)*g(x) on [-1,1]
in terms of their individual Chebyshev series.

V04: Try again.
     Only cout the coefficients for fg, since we know by now that f and g are separately correct.


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
double chebyshev_series_Tn_Tm_integrand( double x, int n, int m );
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
	
	for( int n=0; n<=n_max; n++ ){
		coeffs_f[n] = 2.0/pi * rk4( chebyshev_series_f_integrand, n, 0, pi, n_steps_integration ) * (1-0.5*(n==0));
	}
	
	//--- g(x) ---
	
	for( int n=0; n<=n_max; n++ ){
		coeffs_g[n] = 2.0/pi * rk4( chebyshev_series_g_integrand, n, 0, pi, n_steps_integration ) * (1-0.5*(n==0));
	}
	
	//--- f(x)*g(x) ---
	std::cout << "\nn" <<"\t"<<std::left<<std::setw(w)<< "coeffs[n]" <<std::left<<std::setw(w)<< "a_n_guess(n)" << std::endl;
	
	//--- n=0 is unambiguous at this point ---
	coeffs_fg[0] = coeffs_f[0] * coeffs_g[0];
	for( int n1=1; n1<=n_max; n1++ ){
		coeffs_fg[0] += 0.5 * coeffs_f[n1] * coeffs_g[n1];
	}
	
	//--- n>0 terms --
	
	/*
	//--- V14 ---
	for( int n=1; n<=2*n_max; n++ ){
		for( int n1=0; n1<=n_max; n1++ ){
			for( int n2=0; n2<=n_max; n2++ ){
				coeffs_fg[n] += pow(pi,-1) * coeffs_f[n1] * coeffs_g[n2] * rk4( chebyshev_series_Tn_Tm_integrand, n1+n2, n, 0, pi, n_steps_integration );
				coeffs_fg[n] += pow(pi,-1) * coeffs_f[n1] * coeffs_g[n2] * rk4( chebyshev_series_Tn_Tm_integrand, n1-n2, n, 0, pi, n_steps_integration );
			}
		}
	}
	*/
	
	
	/*
	//--- V15 ---
	for( int n=1; n<=2*n_max; n++ ){
		
		for( int n1=0; n1<=n_max; n1++ ){
			for( int n2=0; n2<=n_max; n2++ ){
				
				//--- First delta function ---
				//if( n1+n2 == n ){
					//coeffs_fg[n] += 0.5 * coeffs_f[n1] * coeffs_g[n2];
				//}
				coeffs_fg[n] += 0.5 * coeffs_f[n1] * coeffs_g[n2] * delta(n1+n2,n);
				
				//--- Second delta function ---
				//if( n1-n2 == n ){
					//coeffs_fg[n] += 0.5 * coeffs_f[n1] * coeffs_g[n2];
				//}
				//if( n1-n2 == -n ){
					//coeffs_fg[n] += 0.5 * coeffs_f[n1] * coeffs_g[n2];
				//}
				coeffs_fg[n] += 0.5 * coeffs_f[n1] * coeffs_g[n2] * ( delta(n1-n2,n) + delta(n2-n1,n) );
				
			}
		}
	}
	*/
	
	/*
	//--- V16 ---
	for( int n=1; n<=2*n_max; n++ ){
		for( int n1=0; n1<=n_max; n1++ ){
			for( int n2=0; n2<=n_max; n2++ ){
				coeffs_fg[n] += 0.5 * coeffs_f[n1] * coeffs_g[n2] * ( delta(n1+n2,n) + delta(n1-n2,n) + delta(n2-n1,n) );
			}
		}
	}
	*/
	
	/*
	//--- V17 ---
	for( int n=1; n<=2*n_max; n++ ){
		for( int n1=0; n1<=n_max; n1++ ){
			for( int n2=0; n2<=n_max; n2++ ){
				coeffs_fg[n] += 0.5 * coeffs_f[n1] * coeffs_g[n2] * delta(n1+n2,n);
			}
		}
	}
	
	for( int n=1; n<=2*n_max; n++ ){
		for( int n1=0; n1<=n_max; n1++ ){
			for( int n2=0; n2<=n_max; n2++ ){
				coeffs_fg[n] += 0.5 * coeffs_f[n1] * coeffs_g[n2] * delta(n1-n2,n);
			}
		}
	}
	
	for( int n=1; n<=2*n_max; n++ ){
		for( int n1=0; n1<=n_max; n1++ ){
			for( int n2=0; n2<=n_max; n2++ ){
				coeffs_fg[n] += 0.5 * coeffs_f[n1] * coeffs_g[n2] * delta(n2-n1,n);
			}
		}
	}
	*/
	
	
	/*
	//--- V18 ---
	for( int n=1; n<=2*n_max; n++ ){
		
		for( int n1=0; n1<=n_max; n1++ ){
			int n2 = n-n1;
			//std::cout << n1 <<"\t"<< n2 << std::endl;
			if( n2 <= n_max ){
				coeffs_fg[n] += 0.5 * coeffs_f[n1] * coeffs_g[n2];
			}
			if( n2 > n_max ){
				std::cout << "index overflow\t" << n2 << std::endl;
			}
		}
	}
		
	for( int n=1; n<=2*n_max; n++ ){
		for( int n1=n; n1<=n_max; n1++ ){
			int n2 = n1-n;
			if( n2 > n_max ){
				std::cout << "index overflow" << std::endl;
			}
			coeffs_fg[n] += 0.5 * coeffs_f[n1] * coeffs_g[n2];
		}
	}
	
	for( int n=1; n<=2*n_max; n++ ){
		for( int n2=n; n2<=n_max; n2++ ){
			int n1 = n2-n;
			if( n1 > n_max ){
				std::cout << "index overflow" << std::endl;
			}
			coeffs_fg[n] += 0.5 * coeffs_f[n1] * coeffs_g[n2];
		}
	}
	*/
	
	//--- V19 ---
	for( int n=1; n<=2*n_max; n++ ){
		
		for( int n1=0; n1<=n_max; n1++ ){
			int n2 = n-n1;
			if( n2 <= n_max ){
				coeffs_fg[n] += 0.5 * coeffs_f[n1] * coeffs_g[n2];
			}
		}
		
		for( int n1=n; n1<=n_max; n1++ ){
			coeffs_fg[n] += 0.5 * ( coeffs_f[n1] * coeffs_g[n1-n] + coeffs_f[n1-n] * coeffs_g[n1] );
		}
		
	}
	
	
	
	
	
	
	
	//--- Output ---
	for( int n=0; n<=2*n_max; n++ ){
		std::cout << n <<"\t"<<std::left<<std::setw(w)<< coeffs_fg[n] <<std::left<<std::setw(w)<< fg_n_guess( n ) << std::endl;
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
	
	for( int n1=0; n1<=n_max; n1++ ){
		for( int n2=0; n2<=n_max; n2++ ){
			fg_from_multiplying += 0.5 * coeffs_f[n1] * coeffs_g[n2] * ( T_n(x,n1+n2) + T_n(x,n1-n2) );
		}
	}
	
	
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