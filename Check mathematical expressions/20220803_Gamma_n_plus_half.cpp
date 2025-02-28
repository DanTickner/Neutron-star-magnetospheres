// cd OneDrive\PhD\Codes\Test mathematical expressions
// g++ Gamma_n_plus_half.cpp -o Gamma_n_plus_half

// Gamma( n + 1/2 ) = sqrt(pi) * (2n)! / ( 2^(2n) n! )

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
		std::cout << n-0.5 << ":\t" << LHS( n ) << "\t" << RHS( n ) << std::endl;
	}
	
	
	return 0;
}




double LHS( int n ){
	return tgamma( n + 0.5 );
}

double RHS( int n ){
	//return sqrt( pi ) * tgamma(2*n+1) / ( pow( 2, 2*n ) * tgamma(n+1) );
	return tgamma( n - 0.5 ) * ( n - 0.5 );
}