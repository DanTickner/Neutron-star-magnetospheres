// cd OneDrive\PhD\Codes\Test mathematical expressions
// g++ 20221122_A_dot_curl_A.cpp -o 20221122_A_dot_curl_A

// Calculate A dot ( curl(A) ), in Cartesian coordinates.

#include <iostream>
#include <math.h>
#include <vector>
#include "../Numerical_Differentiation.h"
#include "../Vector_operations.h"

using namespace std;


//----- Function declarations -----
double A_x_function( double x, double y, double z );
double A_y_function( double x, double y, double z );
double A_z_function( double x, double y, double z );


int main(){
	
	//----- Define variables -----
	double x = 0.123;
	double y = 0.543;
	double z = 0.678;
	
	vector<double> A { x, y, z };
	vector<double> curl_A ( 3 );
	
	
	//----- Calculate curl(A) -----
	curl_A[0] = partial_x2( A_z_function, x, y, z ) - partial_x3( A_y_function, x, y, z );
	curl_A[1] = partial_x3( A_x_function, x, y, z ) - partial_x1( A_z_function, x, y, z );
	curl_A[2] = partial_x1( A_y_function, x, y, z ) - partial_x2( A_x_function, x, y, z );
	
	
	//----- Calculate A dot curl(A) -----
	double A_dot_curl_A = dot_product_cart( A, curl_A ); 
	
	
	//----- Calculate guess for A dot curl(A) -----
	double A_dot_curl_A_guess = A[0] * ( partial_x2( A_z_function, x, y, z ) - partial_x3( A_y_function, x, y, z ) )
	                          + A[1] * ( partial_x3( A_x_function, x, y, z ) - partial_x1( A_z_function, x, y, z ) )
							  + A[2] * ( partial_x1( A_y_function, x, y, z ) - partial_x2( A_x_function, x, y, z ) );
	
	
	print_vector( A, "A\t" );
	print_vector( curl_A, "curl(A)\t" );
	cout << "A dot curl(A):\t" << A_dot_curl_A <<"\t"<< A_dot_curl_A_guess << endl;
	
	
	return 0;
}


double A_x_function( double x, double y, double z ){
	return pow(x,3) * y + sin(z);
}

double A_y_function( double x, double y, double z ){
	return x * y * z;
}

double A_z_function( double x, double y, double z ){
	return cos(x) + y * z;
}