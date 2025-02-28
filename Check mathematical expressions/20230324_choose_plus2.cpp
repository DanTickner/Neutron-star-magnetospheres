/*
cd OneDrive\PhD\Codes\Test mathematical expressions
g++ 20230324_choose_plus2.cpp -o 20230324_choose_plus2
20230324_choose_plus2

( n+2 choose r ) = (n-r)(n-r-1) / ( (r+2)(r+1) ) ( n choose r )

( (n+(k+2)-1)/2 choose n ) = (n+k+1)/(n-k+1) ( (n+k-1)/2 choose )

Book 14, p49.
*/

#include <iostream>
#include <iomanip>
#include <math.h>


//----- Global variables and function declarations -----
double choose( double n, double r );


int main(){
	
	//----- Define variables -----
	int n = 14;
	int r = 5;	// Must be at most n-2.
	
	double exact_1 = choose( n, r+2 );
	double guess_1 = (n-r)*(n-r-1.0) / ( (r+2.0)*(r+1.0) ) * choose( n, r );
	
	std::cout <<std::left<<std::setw(15)<< exact_1 <<std::left<<std::setw(15)<< guess_1 << std::endl;
	
	
	r = 12;
	double exact_2 = choose( 0.5*(n+(r+2)-1), n );
	double guess_2 = (n+r+1.0)/(n-r+1.0) * choose( 0.5*(n+r-1), n );
	
	std::cout <<std::left<<std::setw(15)<< exact_2 <<std::left<<std::setw(15)<< guess_2 << std::endl;
	
	
	
	for( int k=0; k<=n; k++ ){
		std::cout <<std::left<<std::setw(15)<< k <<std::left<<std::setw(15)<< 0.5*(n+k-1)
		          <<std::left<<std::setw(15)<< choose( 0.5*(n+k-1), n ) << std::endl;
	}
	
	return 0;
}


double choose( double n, double r ){
	return tgamma( n +1 ) / ( tgamma( r +1 ) * tgamma( n-r +1 ) );
}