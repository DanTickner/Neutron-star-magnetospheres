/*
cd OneDrive\PhD\Codes\Test mathematical expressions
g++ 20230310_Legendre_polynomial_expression.cpp -o 20230310_Legendre_polynomial_expression
20230310_Legendre_polynomial_expression

P_ell(x) = 2^ell sum_{k=0}^ell x^k ( ell choose k ) ( (ell+k-1)/2 choose ell )
*/

#include <iostream>
#include <math.h>
#include "../Hardcoded_Legendre_Polynomials.h"


//----- Global variables and function declarations -----
const double pi = 3.14159265358979323846;

double choose( int n, int r );
double choose_generalised( double alpha, int r );
double p_expression( double x, int ell );


int main(){
	
	//----- Define variables -----
	double x = 0.4;	// x in (-1,1) for which to test the relation.
	int n_max = 10;	// Maximum value of n
	
	
	//----- Prove that choose() and choose_generalised() functions work -----
	std::cout << choose( 5, 3 ) <<"\t"<< 10 << std::endl;
	std::cout << choose_generalised( 1, 3 ) <<"\t"<< 0 << std::endl;
	std::cout << choose_generalised( 1.5, 3 ) <<"\t"<< -1.0/16.0 << std::endl;
	std::cout << std::endl;
	
	
	
	
	//----- Output results for only the half-integers -----
	std::cout << p_expression( x, 0 ) <<"\t"<< p0( x ) << std::endl;
	std::cout << p_expression( x, 1 ) <<"\t"<< p1( x ) << std::endl;
	std::cout << p_expression( x, 2 ) <<"\t"<< p2( x ) << std::endl;
	std::cout << p_expression( x, 3 ) <<"\t"<< p3( x ) << std::endl;
	std::cout << p_expression( x, 4 ) <<"\t"<< p4( x ) << std::endl;
	std::cout << p_expression( x, 5 ) <<"\t"<< p5( x ) << std::endl;
	std::cout << p_expression( x, 6 ) <<"\t"<< p6( x ) << std::endl;
	std::cout << p_expression( x, 7 ) <<"\t"<< p7( x ) << std::endl;
	std::cout << p_expression( x, 8 ) <<"\t"<< p8( x ) << std::endl;
	std::cout << p_expression( x, 9 ) <<"\t"<< p9( x ) << std::endl;
	
	return 0;
}




double p_expression( double x, int ell ){
	double ret = 0;
	for( int k=0; k<=ell; k++ ){
		
		// V1: Usual definition of choose. Incorrect because second choose has non-integer values
		//ret += pow( x, k ) * choose( ell, k ) * choose( 0.5*(ell+k-1), ell );
		
		// V2: Generalised binomial coefficient in denominator
		ret += pow( x, k ) * choose( ell, k ) * choose_generalised( 0.5*(ell+k-1), ell );
	}
	return pow( 2, ell ) * ret;
}


double choose( int n, int r ){
	// Binomial coefficient.
	return tgamma( n + 1 ) / ( tgamma( r + 1 ) * tgamma( n - r + 1 ) );
}

double choose_generalised( double alpha, int r ){
	// Generalised binomial coefficient.
	// https://en.wikipedia.org/wiki/Binomial_coefficient#Generalization_and_connection_to_the_binomial_series
	
	double ret = 1;
	
	for( int i=0; i<=r-1; i++ ){
		ret *= alpha - (double) i;
	}
	
	return ret / tgamma( r + 1 );
}