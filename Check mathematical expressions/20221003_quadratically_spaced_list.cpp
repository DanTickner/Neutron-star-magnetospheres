// cd OneDrive\PhD\Codes\Test mathematical expressions
// g++ 20221003_quadratically_spaced_list.cpp -o 20221003_quadratically_spaced_list

/*
N quadratically spaced values of x between x_min and x_max, inclusive.
x_i = a + b i + c i^2
We have that
a = x_min
c = ( x_max - x_min ) / (N-1)^2 - b / (N-1)
and b is a free parameter. In this code, we test values of b.
*/

#include <math.h>
#include <vector>
#include <iostream>


int main(){
	
	//----- Define variables -----
	int    N     = 10;
	double x_min = 0;
	double x_max = 1;
	
	std::vector<double> x (N);
	
	//----- Compile the list -----
	//--- Version for any b to try ---
	/*for( int i=0; i<N; i++ ){
		double a = x_min;
		double b = ( x_max - x_min ) * pow( N-1, -2 );
		//double c = ( x_max - x_min ) * pow( N-1, -2 ) - b * pow( N-1, -1 );
		double c = ( x_max - x_min ) * ( N - 2 ) * pow( N-1, -3 );
		x[i] = a + b * i + c * pow(i,2);
	}*/
	
	//--- Version where b can be any value within limits (see book 8, p8) ---
	/*for( int i=0; i<N; i++ ){
		double a = x_min;
		double b = ( x_max - x_min ) * pow( N-1, -1 ) * 0;	// Multiply by any number between 0 and 1, not inclusive. b=1 reduces to linear fit; b=0 reduces to purely squared fit.
		double c = ( x_max - x_min ) * pow( N-1, -2 ) - b * pow( N-1, -1 );
		x[i] = a + b * i + c * pow(i,2);
	}*/
	
	//--- Version with b halfway through above interval
	for( int i=0; i<N; i++ ){
		double a = x_min;
		double b = ( x_max - x_min ) * pow( N-1, -1 ) * 0.5;	// Multiply by any number between 0 and 1, not inclusive. b=1 reduces to linear fit; b=0 reduces to purely squared fit.
		double c = ( x_max - x_min ) * pow( N-1, -2 ) * 0.5;
		//x[i] = a + b * i + c * pow(i,2);
		
		x[i] = x_min + (x_max-x_min)*pow(N-1,-1)*0.5 * i + (x_max-x_min)*pow(N-1,-2)*0.5 * i*i;
	}
	
	//----- Output the list -----
	for( int i=0; i<N; i++ ){
		std::cout << i <<"\t"<< x[i] << std::endl;
	}
	
	return 0;
}