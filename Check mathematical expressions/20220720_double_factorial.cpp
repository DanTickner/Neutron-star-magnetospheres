/*
cd OneDrive\PhD\Codes\Test mathematical expressions
g++ Test_double_factorial.cpp -o Test_double_factorial
Test_double_factorial

n!! = 2^(n/2) (n/2)!                  if n even
n!! = n! / ( 2^((n-1)/2) ((n-1)/2)! ) if n odd

Test the recursive relation
(n-1)!! / n!! = ( 1 - 1/n ) * (n-3)!! / (n-2)!!
*/

#include <iostream>
#include <math.h>


//----- Declare functions -----
int double_factorial( int n );


int main(){
	
	//----- Define variables -----
	int df_odd  = 1;	// Starting point for odd  numbers 1!!
	int df_even = 1;	// Starting point for even numbers 2!!
	
	
	//----- Output results for odd integers -----
	std::cout << "Double factorials of odd integers:" << std::endl;
	for( int n=1; n<=11; n+=2 ){
		df_odd *= n;
		std::cout << n << ":\t" << double_factorial( n ) << "\t" << df_odd << std::endl;
	}
	
	
	//----- Output results for even integers -----
	std::cout << "\nDouble factorials of even integers:" << std::endl;
	for( int n=2; n<=12; n+=2 ){
		df_even *= n;
		std::cout << n << ":\t" << double_factorial( n ) << "\t" << df_even << std::endl;
	}
	
	
	//----- Output recursion relation -----
	std::cout <<"\nRecursion relation for (n-1)!! / n!!" << std::endl;
	for( int n=0; n<=12; n++ ){
		std::cout << n << ":\t" << double_factorial( n-1 ) / ( (double) double_factorial( n ) )
		                <<  "\t" << ( 1.0 - 1.0/n ) * double_factorial( n-3 ) / ( (double) double_factorial( n-2 ) ) << std::endl;
	}
	
	
	return 0;
}

int double_factorial( int n ){
	if( n % 2 == 0 ){
		return pow( 2, 0.5*n ) * tgamma( 0.5*n +1 );
	}
	return tgamma( n +1 ) / ( pow( 2, 0.5*(n-1) ) * tgamma( 0.5*(n-1) +1 ) );
}