// cd OneDrive\PhD\Codes\Test mathematical expressions
// g++ Test_cos_5x.cpp -o Test_cos_5x

#include <math.h>
#include <iostream>

/* Test the equation
cos(5x) = cos^5(x) - 10 cos^3(x)sin^2(x) + 5 cos(x)sin^4(x)
*/


int main(){
	
	//----- Define variables -----
	const double pi = 3.14159265358979323846;
	double x_min = 0;
	double x_max = 2*pi;
	int n_steps = 10;
	double h = ( x_max - x_min ) / ( (double) n_steps - 1 );
	double x = x_min;
	
	//----- Loop through values of x, testing identity for each x -----
	for( int i=0; i<n_steps; i++ ){
		double LHS = cos( 5 * x );
		double RHS = pow( cos(x), 5 ) - 10 * pow( cos(x), 3 ) * pow( sin(x), 2 ) + 5 * cos(x) * pow( sin(x), 4 );
		std::cout << "x = " << x << "\tLHS = " << LHS << "\tRHS = " << RHS << "\tDiff = " << LHS - RHS << std::endl;
		x += h;
	}
	
	return 0;
}