// cd OneDrive\PhD\Codes\Test mathematical expressions
// g++ Test_cos_to_the_n.cpp -o Test_cos_to_the_n

#include <iostream>
#include <complex>
#include <math.h>

/*
Test the following equation
cos^n(x) = Re( 2^-m sum_{k=0}^n (n choose k) exp( i ( 2k-m ) x ) )
         = 2^-m sum_{k=0}^n (n choose k) cos( (2k-m) x )
*/

int main(){
	
	//----- Define variables -----
	int n = 5;
	double x_min = 0;
	double x_max = 6.28;	// 2pi
	int n_terms = 10;
	double h = ( x_max - x_min ) / ( (double) n_terms - 1 );
	double x = 0;
	
	for( int i=0; i<n_terms; i++){
		
		double LHS = pow( cos(x), n );
		
		//std::complex<double> RHS = 0;
		double RHS = 0;
		
		for( int k=0; k<=n; k++ ){
			//RHS += tgamma( n + 1 ) / ( tgamma( k + 1 ) * tgamma( n - k + 1 ) ) * std::complex<double> ( cos( (2*k-n)*x ), sin( (2*k-n)*x ) );
			//RHS += tgamma( n + 1 ) / ( tgamma( k + 1 ) * tgamma( n - k + 1 ) ) * exp( std::complex<double> ( 0, (2*k-n)*x ) );
			RHS += tgamma( n + 1 ) / ( tgamma( k + 1 ) * tgamma( n - k + 1 ) ) * cos( (2*k-n)*x );
		}
		
		RHS *= pow( 2, -n );
		
		std::cout << x << "\t" << LHS << "\t" << RHS << "\t" << LHS - RHS << std::endl;
		
		x += h;
		
	}
	
	return 0;
}