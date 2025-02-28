/*

cd OneDrive\PhD\Codes\VSH series
g++ Wigner_3j_test_20230802.cpp -o Wigner_3j_test_20230802 --std=c++11
Wigner_3j_test_20230802

Test values of the Wigner 3j symbol.
Use values found by the calculator at:
https://www-stone.ch.cam.ac.uk/wigner.shtml

The Wigner 3j symbols are always real. Here, I allow them to be complex while doing code development. When taking values and expressions from this code, you can safely force them to be real.

*/


#include <iostream>
#include <iomanip>
#include <math.h>
#include <vector>
#include <chrono>

#include "../Generate_Spherical_Harmonics.h"	// Includes "../Vector_Operations.h", which we also need here."
#include "../Generate_Associated_Legendre_Functions.h"
#include "../Numerical_Integration.h"




int main(){
	
	int L1 = 1;
	int L2 = 2;
	int L3 = 3;
	int M1 = 1;
	int M2 = 0;
	int M3 = -1;	// Must have M1+M2+M3=0, so must be M3=-M1-M2 for the symbol to be nonzero.
	
	
	std::cout << wigner_3j( 1, 2, 3, 1, 0, -1 ) << std::endl;
	std::cout << std::endl;
	
	std::cout << sqrt(2.0/35.0) << std::endl;
	
	
	std::cout << wigner_3j( L1, L2, L3, M1, M2, M3 ) <<"\t"<< wigner_3j( L2, L3, L1, M2, M3, M1 ) <<"\t"<< wigner_3j( L3, L1, L2, M3, M1, M2 ) << std::endl;
	
	
	//----- Test integral relation for three spherical harmonics -----
	
	generate_coeffs_spherical_harmonic( L1+L2+L3 );
	
	
	int n_points_theta = 100;
	int n_points_phi   = 100;
	
	std::vector<double> theta_list ( n_points_theta );
	std::vector<double> phi_list   ( n_points_phi   );
	
	double delta_theta = 1.0*pi / ( (double) n_points_theta - 1.0 );
	double delta_phi   = 2.0*pi / ( (double) n_points_phi   - 1.0 );
	
	for( int i=0; i<n_points_theta; i++ ){
		theta_list[i] = (double) i * delta_theta;
	}
	
	for( int j=0; j<n_points_phi; j++ ){
		phi_list[j] = (double) j * delta_phi;
	}
	
	std::vector< std::vector< std::complex< double > > > integrand_array ( n_points_theta, std::vector< std::complex< double > > ( n_points_phi ) );
	
	for( int i=0; i<n_points_theta; i++ ){
		for( int j=0; j<n_points_phi; j++ ){
			double t = theta_list[i];
			double p = phi_list  [j];
			integrand_array[i][j] = ylm(t,p,L1,M1) * ylm(t,p,L2,M2) * ylm(t,p,L3,M3);
		}
	}
	
	std::complex<double> integral_of_three_SHs_LHS = integral_spherical_surface( integrand_array, theta_list, phi_list );
	
	std::complex<double> integral_of_three_SHs_RHS  = sqrt( (2.0*L1+1.0)*(2.0*L2+1.0)*(2.0*L3+1.0) / ( 4.0*pi ) ) * wigner_3j( L1, L2, L3, 0, 0, 0 ) * wigner_3j( L1, L2, L3, M1, M2, M3 );
	
	std::cout << "Integral of three SHs:\t" << integral_of_three_SHs_LHS <<"\t"<< integral_of_three_SHs_RHS << std::endl;
	
	
	//----- Test integral relation for three Legendre polynomials -----
	
	//--- Triangle condition ---
	// https://mathworld.wolfram.com/TriangleCondition.html
	
	bool triangle_condition = false;
	
	for( int triangle_condition_L_comparison = abs( L1-L2 ); triangle_condition_L_comparison <= L1+L2; triangle_condition_L_comparison++ ){
		if( L3 == triangle_condition_L_comparison ){
			triangle_condition = true;
		}
	}
	
	std::cout << "Triangle condition:\t" << triangle_condition << std::endl;
	
	
	//--- Build lists of variables and integrands ---
	int n_points_x = 100;
	double delta_x = ( 1.0 - (-1.0) ) / ( (double) n_points_x - 1.0 );
	
	std::vector<double> x_list ( n_points_x );
	
	for( int i=0; i<n_points_x; i++ ){
		x_list[i] = -1.0 + (double) i * delta_x;
	}
	
	std::vector<double> integrand_list_three_LPs_x ( n_points_x );
	for( int i=0; i<n_points_x; i++ ){
		
		double x = x_list[i];
		integrand_list_three_LPs_x[i] = p_ell(x,L1) * p_ell(x,L2) * p_ell(x,L3);
		
		/*
		double P1 = x;
		double P2 = 0.5 * ( 3.0*x*x - 1.0 );
		double P3 = 0.5 * ( 5.0*x*x*x - 3.0*x );
		integrand_list_three_LPs_x[i] = P1 * P2 * P3;
		*/
		
		
	}
	
	std::vector<double> integrand_list_three_LPs_costheta ( n_points_theta );
	for( int i=0; i<n_points_theta; i++ ){
		double t = theta_list[i];
		integrand_list_three_LPs_costheta[i] = p_ell(cos(t),L1) * p_ell(cos(t),L2) * p_ell(cos(t),L3) * sin(t);
	}
	
	//--- Perform and output integrals and expected values ---
	std::complex<double> integral_of_three_LPs_LHS_x = trapezium( integrand_list_three_LPs_x, x_list );				// Should be real but leave as complex until I know it is working.
	std::complex<double> integral_of_three_LPs_LHS_costheta = trapezium( integrand_list_three_LPs_costheta, theta_list );
	std::complex<double> integral_of_three_LPs_RHS = 2.0 * pow( wigner_3j( L1, L2, L3, 0, 0, 0 ), 2 );
	
	std::cout << "Integral of three LPs:\t" << integral_of_three_LPs_LHS_x <<"\t"<< integral_of_three_LPs_LHS_costheta <<"\t"<< integral_of_three_LPs_RHS << std::endl;
	
	
	
	
	//----- Product of two associated Legendre functions -----
	// The Overlap Integral of Three Associated Legendre Polynomials, Dong and Lemus, 2002, Applied Mathematics Letters, 15, Eqs. (2) and (3).
	
	generate_coeffs_associated_legendre( L1+L2+L3 );
	
	double costheta = cos( theta_list[43] );
	double product_of_two_ALFs_LHS = associated_legendre_function(costheta,L1,M1) * associated_legendre_function(costheta,L2,M2);
	
	double product_of_two_ALFs_RHS = 0;
	
	for( int ell=abs(L1-L2); ell<=L1+L2; ell++ ){
		double G = ( pow( -1, -M1+M2 ) * ( 2.0*ell+1.0 ) * wigner_3j( L1, L2, ell, 0, 0, 0 ) * wigner_3j( L1, L2, ell, -M1, M2, M1-M2 ) ).real();
		double sqrt_factor = sqrt( tgamma( ell-M2+M1 +1 ) / tgamma( ell+M2-M1 +1 ) );
		double ALF_ell = associated_legendre_function( costheta, ell, -M1+M2 );
		product_of_two_ALFs_RHS += G * sqrt_factor * ALF_ell;
	}
	product_of_two_ALFs_RHS *= pow( -1, M1 ) * sqrt( tgamma( L1+M1 +1 ) * tgamma( L2+M2 +1 ) / ( tgamma( L1-M1 + 1 ) * tgamma( L2-M2 + 1 ) ) );
	
	
	std::cout << "Product of two ALFs\t" << product_of_two_ALFs_LHS <<"\t"<< product_of_two_ALFs_RHS << std::endl;
	
	
	
	
	//-----
	
	
	
	
	
	
	
	
	
			
	
	
	
	std::cout << "done" << std::endl;
	return 0;
}