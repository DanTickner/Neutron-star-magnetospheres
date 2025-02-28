// cd OneDrive\PhD\Codes\Test mathematical expressions
// g++ 20221122_curl_of_vector.cpp -o 20221122_curl_of_vector

// Calculate the curl of a vector in Cartesian coordinates, ahead of testing a vector calculus identity.

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
	double x_test = 0.123;
	double y_test = 0.543;
	double z_test = 0.678;
	
	vector<double> A_test { x_test, y_test, z_test };
	vector<double> curl_A_test ( 3 );
	
	//----- Calculate curl(A) -----
	curl_A_test[0] = partial_x2( A_z_function, x_test, y_test, z_test ) - partial_x3( A_y_function, x_test, y_test, z_test );
	curl_A_test[1] = partial_x3( A_x_function, x_test, y_test, z_test ) - partial_x1( A_z_function, x_test, y_test, z_test );
	curl_A_test[2] = partial_x1( A_y_function, x_test, y_test, z_test ) - partial_x2( A_x_function, x_test, y_test, z_test );
	
	print_vector( A_test, "A" );
	print_vector( curl_A_test, "curl(A)" );
	
	
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