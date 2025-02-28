// cd OneDrive\PhD\Codes\Test mathematical expressions
// g++ minus_one_to_the_minus_n.cpp -o minus_one_to_the_minus_n

// (-1)^-n = (-1)^n

#include <iostream>
#include <math.h>


//----- Global variables and function declarations -----
const double pi = 3.14159265358979323846;

double LHS( int n );
double RHS( int n );


int main(){
	
	//----- Define variables -----
	int n_min = -10;	// Minimum value of n
	int n_max =  10;	// Maximum value of n
	
	
	//----- Output results for increasing n -----
	for( int n=1; n<=n_max; n++ ){
		std::cout << n <<":\t"<< LHS( n ) <<"\t"<< RHS( n ) <<"\t"<< LHS( n ) - RHS( n ) << std::endl;
	}
	
	
	return 0;
}




double LHS( int n ){
	return pow( -1, -n );
}

double RHS( int n ){
	return pow( -1, n );
}