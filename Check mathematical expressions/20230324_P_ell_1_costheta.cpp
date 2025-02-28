/*
cd OneDrive\PhD\Codes\Test mathematical expressions
g++ 20230324_P_ell_1_costheta.cpp -o 20230324_P_ell_1_costheta
20230324_P_ell_1_costheta

P_ell^1[cos(theta)] = -2^ell sin(theta) sum_{k=1}^ell k {ell choose k} {(ell+k-1)/2 choose ell} cos^{k-1}(theta)
On paper, the choose-terms are zero for some k. In codes, they return nan, so we must skip over them. Hence the k=ell%2 and k+=2.
*/

#include <iostream>
#include <iomanip>
#include <math.h>
#include "../Generate_Associated_Legendre_Functions.h"


//----- Global variables and function declarations -----
double choose( double n, double r );


int main(){
	
	//----- Define variables -----
	double theta = 0.321 * pi;
	
	int ell_max = 30;	// Maximum value of ell
	
	double x = cos( theta );
	
	
	//----- Generate expression for P_ell^1[cos(theta)] -----
	
	generate_coeffs_associated_legendre( ell_max );
	
	std::cout <<std::left<<std::setw( 5)<< "ell" <<std::left<<std::setw(15)<< "P_test" <<std::left<<std::setw(15)<< "P_exact"
                  <<std::left<<std::setw(15)<< "P_test/P_exact" << std::endl;
	
	for( int ell=1; ell<=ell_max; ell++ ){
		double P_test = 0;
		for( int k=ell%2; k<=ell; k+=2 ){
			P_test += k * choose( ell, k ) * choose( 0.5*(ell+k-1), ell ) * pow( x, k-1 );
			//std::cout << choose( ell, k ) <<"\t"<< choose( 0.5*(ell+k-1), ell ) <<"\t"<< pow( x, k-1 ) <<"\t"<< P_test << std::endl;
		}
		
		P_test *= - pow(2,ell) * sin(theta);
		double P_exact = associated_legendre_function( x, ell, 1 );
		
		std::cout <<std::left<<std::setw( 5)<< ell <<std::left<<std::setw(15)<< P_test <<std::left<<std::setw(15)<< P_exact
                  <<std::left<<std::setw(15)<< P_test/P_exact << std::endl;
	}
	/*
	
	
	
	std::vector<double> sequence_explicit ( ell-1 );
	std::vector<double> sequence_recursion( ell-1 );
	
	
	//----- Calculate the sequence by the explicit method -----
	for( int k=0; k<=k_max; k++ ){
		sequence_explicit[k] = k * choose( ell, k ) * choose( 0.5*(ell+k-1), ell );
	}
	
	
	//----- Output results -----
	std::cout <<std::left<<std::setw( 5)<< "ell" <<std::left<<std::setw(15)<< "Recursion" <<std::left<<std::setw(15)<< "Explicit"
              <<std::left<<std::setw(15)<< "Ratio"<< std::endl;
	
	for( int k=0; k<=k_max; k++ ){
		
		double P_explicit = -pow(2,ell)
		
		std::cout <<std::left<<std::setw( 5)<< k <<std::left<<std::setw(15)<< sequence_recursion[k] <<std::left<<std::setw(15)<< sequence_explicit[k]
                  <<std::left<<std::setw(15)<< sequence_recursion[k] / sequence_explicit[k] << std::endl;
	}
	
	std::cout << "Done" << std::endl;
	*/
	return 0;
}


double choose( double n, double r ){
	return tgamma( n +1 ) / ( tgamma( r +1 ) * tgamma( n-r +1 ) );
}