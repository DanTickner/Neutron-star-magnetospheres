// cd OneDrive\PhD\Codes\Test mathematical expressions
// g++ Test_fraction_of_gamma_functions_3.cpp -o Test_fraction_of_gamma_functions_3

// Gamma( n + 1/2 ) / ( 2 Gamma( +1/2) ) = (2n)! / ( 2^(2n+2) n! )

#include <iostream>
#include <math.h>


//----- Global variables and function declarations -----
const double pi = 3.14159265358979323846;

double LHS( int n );
double RHS( int n );


int main(){
	
	//----- Define variables -----
	int n_max = 10;	// Maximum value of n
	
	
	//----- Output results for only the half-integers -----
	for( int n=1; n<=n_max; n+=2 ){
		std::cout << n+0.5 << ":\t" << LHS( n ) << "\t" << RHS( n ) << std::endl;
	}
	
	
	return 0;
}

double LHS( int n ){
	return tgamma( n + 0.5 ) / ( 2 * tgamma( 0.5 ) );
}

double RHS( int n ){
	return tgamma(2*n+1) / ( pow( 2, 2*n+1 ) * tgamma( n+1 ) );
}