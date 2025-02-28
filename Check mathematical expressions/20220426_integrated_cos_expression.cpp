// cd OneDrive\PhD\Codes\Test mathematical expressions
// g++ Test_integrated_cos_expression.cpp -o Test_integrated_cos_expression

#include <math.h>
#include <iostream>

/* Test the equation
0.5 * (  [ (-1)^(a+b) - 1 ] / [ a + b ] + [ (-1)^(a-b) - 1 ] / [ a - b ] ) = 0 for all a, b
*/


int main(){
	
	//----- Define variables -----
	const double pi = 3.14159265358979323846;
	int a_max = 7;
	int b_max = 7;
	
	for( int a=1; a<=a_max; a++ ){
		for( int b=1; b<=b_max; b++ ){
			if( a != b ) {
				//double x = 0.5 * ( ( pow( -1, a+b ) - 1 ) / ( (double) a + b ) + ( pow( -1, a-b ) - 1 ) / ( (double) a - b ) );
				double x = 0.5 * ( ( cos((a+b)*pi) - cos((a+b)*0) ) / ( (double) a + b ) + ( cos((a-b)*pi) - cos((a-b)*0) ) / ( (double) a - b ) );
				std::cout << "a = " << a << "\tb = " << b << "\tx = " << x << std::endl;
			}
		}
		std::cout << std::endl;
	}
	
	return 0;
}