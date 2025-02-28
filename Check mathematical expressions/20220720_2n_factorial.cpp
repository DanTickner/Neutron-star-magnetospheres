// cd OneDrive\PhD\Codes\Test mathematical expressions
// g++ Test_2n_factorial.cpp -o Test_2n_factorial

// (2n)! = 2^n n!
// This expression has been proven incorrect by the code.

#include <iostream>
#include <math.h>


//----- Declare functions -----
int double_factorial( int n );


int main(){
	
	//----- Define variables -----
	int n_max = 10;	// Maximum value of n
	
	
	//----- Output results for odd integers -----
	for( int n=1; n<=11; n++ ){
		std::cout << n << ":\t" << tgamma( 2*n+1 ) << "\t" << pow( 2, n ) * tgamma( n+1 ) << std::endl;
	}
	
	
	return 0;
}