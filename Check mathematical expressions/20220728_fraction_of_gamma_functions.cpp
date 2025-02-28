// cd OneDrive\PhD\Codes\Test mathematical expressions
// g++ Test_fraction_of_gamma_functions.cpp -o Test_fraction_of_gamma_functions

// - Gamma( ell/2 ) / Gamma( - ell/2 ) = ell / 2^(2ell-5) * ( (ell-2)! / ((ell-3)/2)! )^2 x { -1 if ell=1 mod 4, or 1 if ell=3 mod 4 }. But if ell=1 then instead the answer is 1/2.

#include <iostream>
#include <math.h>


//----- Global variables and function declarations -----
const double pi = 3.14159265358979323846;

double gamma_half_integer( int n );


int main(){
	
	//----- Define variables -----
	int n_max = 10;	// Maximum value of n
	
	
	//----- Output results for only the half-integers -----
	for( int n=1; n<=n_max; n+=2 ){
		std::cout << 0.5*n << ":\t" << - tgamma( 0.5*n ) / tgamma( -0.5*n ) << "\t" << gamma_half_integer( n ) << std::endl;
	}
	
	
	return 0;
}

double gamma_half_integer( int n ){
	/*// V1. See Book 6, p51.
	double ret_positive = sqrt(pi) * tgamma(n-2+1) / ( pow(2,n-2) * tgamma(0.5*(n-3)+1) );
	double ret_negative = sqrt(pi) * pow(2,n-1) * pow(n,-1) * tgamma(0.5*(n-3)+1) / tgamma(n-2+1);
	if( n % 4 == 1 ){
		ret_negative *= -1;
	}
	return - ret_positive / ret_negative;*/
	
	/*// V2. Temporarily recombine the +- cases into a sine function for easier function manipulation.
	double ret_positive = sqrt(pi) * tgamma(n-2+1) / ( pow(2,n-2) * tgamma(0.5*(n-3)+1) );
	double ret_negative = - sin(n*0.5*pi) * sqrt(pi) * pow(2,n-1) * pow(n,-1) * tgamma(0.5*(n-3)+1) / tgamma(n-2+1);
	return - ret_positive / ret_negative;*/
	
	/*// V3. Write Gamma(-n/2) in terms of Gamma(n/2).
	double ret_positive = sqrt(pi) * tgamma(n-2+1) / ( pow(2,n-2) * tgamma(0.5*(n-3)+1) );
	double ret_negative = - 2 * pi / ( sin(n*0.5*pi) * ret_positive * n );
	return - ret_positive / ret_negative;*/
	
	/*// V4. See whiteboard photo 20220729.
	double ret_positive = sqrt(pi) * tgamma(n-2+1) / ( pow(2,n-2) * tgamma(0.5*(n-3)+1) );
	return pow( ret_positive, 2 ) * n * sin(0.5*pi*n) / ( 2 * pi );*/
	
	// V5.
	//return ( sin(0.5*n*pi) * n / (2*pi) ) * pow( sqrt(pi)/pow(2,n-2) * tgamma(n-2) / tgamma(0.5*(n-3)), 2 );
	//return pow( sqrt(pi) * tgamma(n-2+1) / ( pow(2,n-2) * tgamma(0.5*(n-3)+1) ), 2 ) * n * sin(0.5*pi*n) / ( 2 * pi );
	//return pow( tgamma(n-2+1) / ( pow(2,n-2) * tgamma(0.5*(n-3)+1) ), 2 ) * n * sin(0.5*pi*n) / ( (double) 2 );
	//return pow( tgamma(n-2+1) / ( tgamma(0.5*(n-3)+1) ), 2 ) * pow( pow(2,n-2), -2 ) * n * sin(0.5*pi*n) / ( (double) 2 );
	//return pow( tgamma(n-2+1) / ( tgamma(0.5*(n-3)+1) ), 2 ) * pow(2,-(2*n-4)) * n * sin(0.5*pi*n) / ( (double) 2 );
	if( n == 1 ){
		return 0.5;
	}
	return pow( tgamma(n-2+1) / ( tgamma(0.5*(n-3)+1) ), 2 ) * n * sin(0.5*pi*n) / pow(2,2*n-3);
}