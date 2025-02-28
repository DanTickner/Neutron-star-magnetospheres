// cd OneDrive\PhD\Codes\Test mathematical expressions
// g++ Test_gamma_negative_half_integer.cpp -o Test_gamma_negative_half_integer

// Gamma( -n/2 ) = sqrt(pi) * 2^(n-1)/n * ((n-3)/2)! / (n-2)! x { -1 if n=1 mod 4, or 1 if n=3 mod 4 }.

#include <iostream>
#include <math.h>


//----- Global variables and function declarations -----
const double pi = 3.14159265358979323846;

int double_factorial( int n );
double gamma_half_integer( int n );
int double_factorial( int n );


int main(){
	
	//----- Define variables -----
	int n_max = 10;	// Maximum value of n
	
	
	//----- Output results for only the half-integers -----
	for( int n=1; n<=n_max; n+=2 ){
		std::cout << 0.5*n << ":\t" << tgamma( -0.5*n ) << "\t" << gamma_half_integer( n ) << std::endl;
	}
	
	
	return 0;
}

double gamma_half_integer( int n ){
	// V1. See Book 6, p48, p50.
	//return pi/sin(pi * -0.5*n) / tgamma(1+0.5*n);
	//return pi/sin(pi * -0.5*n) / ( 0.5*n * tgamma(0.5*n) );
	//return - pi/sin(pi*0.5*n) * (2.0/((double)n)) / (  tgamma(0.5*n) );
	return -2*pi / ( sin(0.5*n*pi) * n * tgamma(0.5*n) );
	
	// V2. Absorb the gamma function, using that sin(pi*n/2) is +1 if n=1,5,9... or -1 if n=3,7,11,..., and that we already have a minus sign due to sin(x) being odd.
	/*double ret = sqrt(pi) * pow(2,n-1) * pow(n,-1) * tgamma(0.5*(n-3)+1) / tgamma(n-2+1);
	if( n % 4 == 1 ){
		return -ret;
	}
	if( n % 4 == 3 ){
		return ret;
	}
	else{
		std::cout << "There is a mistake." << std::endl;
		return 111111l;
	}*/
	
	
	//return sqrt(pi) * double_factorial( n - 2 ) / pow( 2, 0.5*(n-1) );
	
	// V2. n must be odd. Evaluate the gamma function at n-2.
	//return sqrt(pi) * tgamma(n-2+1) / ( pow(2,0.5*(n-1)) * pow(2,0.5*(n-3)) * tgamma(0.5*(n-3)+1) );
	
	// V3. Cleanup.
	//return sqrt(pi) * tgamma(n-2+1) / ( pow(2,n-2) * tgamma(0.5*(n-3)+1) );
}




int double_factorial( int n ){
	if( n % 2 == 0 ){
		return pow( 2, 0.5*n ) * tgamma( 0.5*n + 1 );
	}
	return tgamma( n + 1 ) / ( pow( 2, 0.5*(n-1) ) * tgamma( 0.5*(n+1) ) );
}