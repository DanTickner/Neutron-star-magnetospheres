/*
cd OneDrive\PhD\Codes\Test mathematical expressions
g++ 20230324_choose_plus2_02.cpp -o 20230324_choose_plus2_02
20230324_choose_plus2_02

( ell choose k ) = (ell-k+2)(ell-k+1)/(k(k-1)) ( ell choose k-2 )

( (ell+k-1)/2 choose ell = (ell+k-1)/(-ell+k-1) ( (ell+[k-2]-1)/2 choose ell )

( ell choose k ) = ell(ell-k+1)/(k(k-1)) ( ell-1 choose k-2 )

( (ell-1)/2 choose ell ) = 1/4 ( 1/ell - 1 ) ( ([ell-2]-1)/2 choose ell-2 )

A4 paper, photographed 20230324.
*/

#include <iostream>
#include <iomanip>
#include <math.h>


//----- Global variables and function declarations -----
double choose( double ell, double k );


int main(){
	
	//----- Define variables -----
	int ell = 3;
	int k = 4;	// Must be at most ell-2.
	
	double exact_1 = choose( ell, k );
	double guess_1 = (ell-k+2.0)*(ell-k+1.0) / ( (double) k*(k-1) ) * choose( ell, k-2 );
	
	std::cout <<std::left<<std::setw(15)<< exact_1 <<std::left<<std::setw(15)<< guess_1 << std::endl;
	
	
	
	double exact_2 = choose( 0.5*(ell+k-1), ell );
	double guess_2 = (ell+k-1.0) / ( (double) -ell+k-1.0 ) * choose( 0.5*(ell+(k-2)-1), ell );
	
	std::cout <<std::left<<std::setw(15)<< exact_2 <<std::left<<std::setw(15)<< guess_2 << std::endl;
	
	
	
	double exact_tot = choose( ell, k ) * choose( 0.5*(ell+k-1), ell );
	double guess_tot = (ell-k+2.0)*(ell-k+1.0)*(ell+k-1.0) / ( (double) k*(k-1)*(-ell+k-1.0) ) * choose( ell, k-2 ) * choose( 0.5*(ell+(k-2)-1), ell );
	
	std::cout <<std::left<<std::setw(15)<< exact_tot <<std::left<<std::setw(15)<< guess_tot << std::endl;
	
	
	double exact_3 = choose( ell, k );
	double guess_3 = ell*(ell-k+1.0) / ( (double) k*(k-1.0) ) * choose( ell-1, k-2 );
	
	std::cout <<std::left<<std::setw(15)<< exact_3 <<std::left<<std::setw(15)<< guess_3 << std::endl;
	
	double exact_4 = choose( 0.5*(ell-1), ell );
	double guess_4 = 0.25 * ( 1.0/ell - 1.0 ) * choose( 0.5*((ell-2)-1), ell-2 );
	
	std::cout <<std::left<<std::setw(15)<< exact_4 <<std::left<<std::setw(15)<< guess_4 << std::endl;
	
	return 0;
}


double choose( double ell, double k ){
	return tgamma( ell +1 ) / ( tgamma( k +1 ) * tgamma( ell-k +1 ) );
}