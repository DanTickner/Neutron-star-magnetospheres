// cd OneDrive\PhD\Codes\Test mathematical expressions
// g++ 20221213_Chebyshev_and_derivatives_in_terms_cos_and_sin.cpp -o 20221213_Chebyshev_and_derivatives_in_terms_cos_and_sin

/* To save having to generate the coefficients of Chebyshev polynomials, test the following direct expressions:
T_n(x) = cos( n acos(x) )
U_n(x) = sin( (n+1) acos(x) ) / sin( acos(x) ) = sin( (n+1) acos(x) ) / sqrt( 1-x^2 )
d/dx T_n(x) = n U_{n-1}(x) = n sin( n acos(x) ) / sqrt( 1-x^2 )
*/


#include <iostream>
#include <math.h>
#include "../Generate_Chebyshev_Polynomials.h"


//----- Global variables and function declarations -----
double LHS( int n );
double RHS( int n );


int main(){
	
	//----- Define variables -----
	int    n = 3;
	double x = 0.4;
	
	
	//----- Generate basis functions -----
	generate_coeffs_chebyshev( n );
	
	
	//----- Evaluate the functions -----
	double Tn_exact = polynomial_chebyshev_1( x, n );
	double Tn_guess = cos( n * acos(x) );
	
	double Un_exact   = polynomial_chebyshev_2( x, n );
	double Un_guess_1 = sin( (n+1) * acos(x) ) / sin( acos(x) );
	double Un_guess_2 = sin( (n+1) * acos(x) ) / sqrt( 1-x*x );
	
	double Tn_dx_exact = polynomial_chebyshev_1_dr( x, n );
	double Tn_dx_guess = n * sin( n * acos(x) ) / sqrt( 1-x*x );
	
	std::cout << "Tn   :\t" << Tn_exact    <<"\t"<< Tn_guess      << std::endl;
	std::cout << "Un   :\t" << Un_exact    <<"\t"<< Un_guess_1    <<"\t"<< Un_guess_2 << std::endl;
	std::cout << "Tn_dx:\t" << Tn_dx_exact <<"\t"<< Tn_dx_guess   << std::endl;
	
	
	return 0;
}