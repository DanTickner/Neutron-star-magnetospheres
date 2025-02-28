/*
cd OneDrive\PhD\Codes\Check mathematical expressions
g++ 20230624_Chebyshev_polynomial_second_kind_explicit_power_series.cpp -o 20230624_Chebyshev_polynomial_second_kind_explicit_power_series
20230624_Chebyshev_polynomial_second_kind_explicit_power_series


U_n(x) = sum_{k=0}^floor(n/2) (-1)^k (n-k)! / ( k! (n-2k)! ) 2^(n-2k) x^(n-2k)

Book 19, p4 and Eq (1.10) of
https://hal.science/hal-01705040/document
which is quoting from pp.432-433 of
D. Zwillinger, CRC Standard Mathematical Tables and Formulae, Thirty-first edition, CRC Press, Boca Raton, FL, 2003.

*/

#include <iostream>
#include <iomanip>
#include <math.h>

double Un( double x, int n );
double Un_hardcoded( double x, int n );
double Un_powerseries( double x, int n );


int main(){
	
	//----- Define variables -----
	int n_max = 12;
	double x = 0.345;
	int w = 16;			// Fixed width.
	
	
	std::cout <<std::left<<std::setw(w)<< "n" <<std::left<<std::setw(w)<< "LHS" <<std::left<<std::setw(w)<< "RHS" <<std::left<<std::setw(w)<< "Difference" << std::endl;
	
	
	for( int n=0; n<=n_max; n++ ){
		
		double LHS = Un( x, n );
		double RHS = Un_powerseries( x, n );
		
		std::cout <<std::left<<std::setw(w)<< n <<std::left<<std::setw(w)<< LHS <<std::left<<std::setw(w)<< RHS <<std::left<<std::setw(w)<< abs( LHS - RHS ) << std::endl;
		
	}
	
	
	return 0;
	
}




double Un( double x, int n ){
	
	double endpoint_tolerance = 1e-5;	// Freely adjustable.
	
	if( abs( x - 1 ) < endpoint_tolerance ){
		//U_n(1) = n+1
		return (double) n+1.0;
	}
	
	if( abs( x - (-1) ) < endpoint_tolerance ){
		// U_n(-1) = (-1)^n (n+1).
		return pow(-1,n) * (double) n+1.0;
	}
	
	return sin( (double)(n+1.0) * acos(x) ) / sqrt( 1.0 -x*x );
	
}




double Un_hardcoded( double x, int n ){
	
	switch( n ){
		case 0:
			return 1.0;
			break;
		case 1:
			return 2.0*x;
			break;
		case 2:
			return 4.0*x*x - 1.0;
			break;
		case 3:
			return 8.0*pow(x,3) - 4.0*x;
			break;
		
		
		default:
			return 0;
	}
	
}




double Un_powerseries( double x, int n ){
	
	double Un = 0;
	
	for( int k=0; k<=n/2; k++ ){
		double a_k = pow(-1,k) * pow( 2, n-2*k ) * tgamma( n-k +1 ) / ( tgamma( k +1 ) * tgamma( n-2*k +1 ) );
		Un += a_k * pow( x, n-2*k );
	}
	
	return Un;
}