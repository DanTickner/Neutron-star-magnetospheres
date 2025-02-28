/*
cd OneDrive\PhD\Codes\Check mathematical expressions
g++ 20230726_Plm_over_sintheta.cpp -o 20230726_Plm_over_sintheta
20230726_Plm_over_sintheta


P_L^m[cos(theta)] * 1/sin(theta) = -1/(2m) P_{L+1}^{m+1} - (L-m+1)(L-m+2)/(2m) P_{L+1}^{m-1}

OneNote Force-Free Documentation Proofs -> Spherical harmonics
*/

#include <iostream>
#include <iomanip>
#include <math.h>

#include "../Generate_Associated_Legendre_Functions.h"
#include "../Hardcoded_Legendre_Functions.h"


int main(){
	
	//----- Define variables -----
	int n_max = 12;
	double x = 0.345;
	int w = 40;			// Fixed width.
	
	int L_max = 6;
	
	generate_coeffs_associated_legendre( L_max );
	//std::cout << ylm( 0.345, 0.432, 3, 2 ) <<"\t"<< y32( 0.345, 0.432 ) << std::endl;
	
	double theta = 0.345;
	double phi   = 0.432;
	
	std::cout << "L" <<"\t"<< "m" <<"\t|"<<std::left<<std::setw(w)<< "LHS" <<std::left<<std::setw(w)<< "RHS" <<std::left<<std::setw(w)<< "Difference" << std::endl;
	
	for( int L=0; L<=L_max+1; L++ ){
		// L_max+1 because otherwise RHS cannot be evaluated for L=L_max.
		for( int m=-L; m<=L; m++ ){
			
			double LHS = associated_legendre_function( cos(theta), L, m ) / sin( theta );
			
			double RHS = -1.0/(2.0*m) * associated_legendre_function( cos(theta), L+1, m+1 ) - (L-m+1.0)*(L-m+2.0)/(2.0*m) * associated_legendre_function( cos(theta), L+1, m-1 );
			
			std::cout << L <<"\t"<< m <<"\t|"<<std::left<<std::setw(w)<< LHS <<std::left<<std::setw(w)<< RHS <<std::left<<std::setw(w)<< abs( LHS - RHS ) << std::endl;
			
		}
		std::cout << std::endl;
	}
	
	std::cout << "\ndone" << std::endl;
	
	
	return 0;
	
}