/*
cd OneDrive\PhD\Codes\Check mathematical expressions
g++ 20230726_Ylm_dphi_over_sintheta.cpp -o 20230726_Ylm_dphi_over_sintheta
20230726_Ylm_dphi_over_sintheta


Y_L^m * 1/sin(theta) = -1/(2m) sqrt( (2L+1)/(2L+3) ) ( sqrt( (L+m+2)(L+m+1) ) Y_{L+1}^{m-1} e^{-i*phi} + sqrt( (L-m+2)(L-m+1) ) Y_{L+1}^{m-1} e^{i*phi} )

OneNote Force-Free Documentation Proofs -> Spherical harmonics
*/

#include <iostream>
#include <iomanip>
#include <math.h>

#include "../Generate_Spherical_Harmonics.h"
#include "../Hardcoded_Spherical_Harmonics.h"


int main(){
	
	//----- Define variables -----
	int n_max = 12;
	double x = 0.345;
	int w = 40;			// Fixed width.
	
	int L_max = 6;
	
	generate_coeffs_spherical_harmonic( L_max );
	//std::cout << ylm( 0.345, 0.432, 3, 2 ) <<"\t"<< y32( 0.345, 0.432 ) << std::endl;
	
	double theta = 0.345;
	double phi   = 0.432;
	double h = 1.0e-5;		// Small positive number for numerical differentiation.
	
	std::cout << "L" <<"\t"<< "m" <<"\t|"<<std::left<<std::setw(w)<< "LHS" <<std::left<<std::setw(w)<< "RHS" <<std::left<<std::setw(w)<< "Difference" << std::endl;
	
	for( int L=0; L<=L_max+1; L++ ){
		// L_max+1 because otherwise RHS cannot be evaluated for L=L_max.
		for( int m=-L; m<=L; m++ ){
			
			std::complex<double> ylm_dphi = ( ylm( theta, phi+h, L, m ) - ylm( theta, phi-h, L, m ) ) / ( 2.0*h );
			std::complex<double> LHS = ylm_dphi / sin( theta );
			
			/*
			std::complex<double> RHS_a = sqrt( (L+m+2.0)*(L+m+1.0) ) * ylm( theta, phi, L+1, m+1 ) * exp( std::complex<double>{0,-phi} );
			std::complex<double> RHS_b = sqrt( (L-m+2.0)*(L-m+1.0) ) * ylm( theta, phi, L+1, m-1 ) * exp( std::complex<double>{0, phi} );
			std::complex<double> RHS = -0.5 * std::complex<double>{0,1} * sqrt( (2.0*L+1.0)/((double)(2.0*L+3.0)) ) * ( RHS_a + RHS_b );
			*/
			
			std::complex<double> RHS = ylm_dphi_over_sintheta( theta, phi, L, m );
			
			
			std::cout << L <<"\t"<< m <<"\t|"<<std::left<<std::setw(w)<< LHS <<std::left<<std::setw(w)<< RHS <<std::left<<std::setw(w)<< abs( LHS - RHS ) << std::endl;
			
		}
		std::cout << std::endl;
	}
	
	std::cout << "\ndone" << std::endl;
	
	
	return 0;
	
}