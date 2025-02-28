/*
cd C:\Users\Dan\OneDrive\PhD\Codes\Check mathematical expressions
g++ test_Gaussian_elimination_20240413.cpp -o test_Gaussian_elimination_20240413 --std=c++11
test_Gaussian_elimination_20240413

Test the Gaussian elimination code by solving a system of equations with known result.

2x + y -  z = 5
3x - y + 2z = -1
 x - y -  z = 0
 
x = 1
y = 2
z = -1
*/


#include <vector>
#include <iostream>

#include "../Gaussian_elimination.h"


int main(){
	
	
	//----- Define variables -----
	int n_variables = 3;
	
	std::vector< std::vector<double> > A ( n_variables, std::vector<double> ( n_variables ) );
	std::vector<double> b ( n_variables );
	
	//----- Build matrix A and vector b -----
	
	A[0][0] =   2.0;
	A[0][1] =   1.0;
	A[0][2] =  -1.0;
	
	A[1][0] =  3.0;
	A[1][1] = -1.0;
	A[1][2] =  2.0;
	
	A[2][0] =  1.0;
	A[2][1] = -1.0;
	A[2][2] = -1.0;
	
	b[0] =  5.0;
	b[1] = -1.0;
	b[2] =  0.0;
	
	//----- Perform Gaussian elimination -----
	
	std::vector<double> x = gaussian_elimination( A, b );
	
	//----- Output results -----
	for( int i=0; i<n_variables; i++ ){
		std::cout << "x_" << i << ":\t" << x[i] << std::endl;
	}
	
	
	std::cout << "Code finished";
	
	return 0;
}