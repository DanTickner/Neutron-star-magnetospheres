// cd OneDrive\PhD\Codes\Test mathematical expressions
// g++ Test_cos_identity.cpp -o Test_cos_identity

#include <math.h>
#include <iostream>

/* Test the equation
cos(ax) cos(bx) = 0.5 ( cos( (a+b)x ) + cos( (a-b)x )
*/


int main(){
	
	//----- Define variables -----
	const double pi = 3.14159265358979323846;
	int a_max = 7;
	int b_max = 7;
	double x = 0.32;	// One value of x.
	
	for( int a=1; a<=a_max; a++ ){
		for( int b=1; b<=b_max; b++ ){
			double LHS = cos( a*x ) * cos( b*x );
			double RHS = 0.5 * ( cos( (a+b)*x ) + cos( (a-b)*x ) );
			std::cout << "a = " << a << "\tb = " << b << "\tLHS = " << LHS << "\tRHS = " << RHS << "\tDiff = " << LHS - RHS << std::endl;
		}
		std::cout << std::endl;
	}
	
	return 0;
}