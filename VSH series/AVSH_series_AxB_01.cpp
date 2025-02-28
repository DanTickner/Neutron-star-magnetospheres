/*

cd OneDrive\PhD\Codes\VSH series
g++ AVSH_series_AxB_01.cpp -o AVSH_series_AxB_01 --std=c++11
AVSH_series_AxB_01

Read two vector functions
A(r,theta) = A_r(r,theta) e_r + A_theta(r,theta) e_theta + A_phi(r,theta) e_phi
B(r,theta) = B_r(r,theta) e_r + B_theta(r,theta) e_theta + B_phi(r,theta) e_phi
and calculate the AVSH (axisymmetric vector spherical harmonic) series coefficients of AxB at a given value of r.

Based on AVSH_series_AxB_02.cpp.

*/


#include <iostream>
#include <iomanip>
#include <math.h>
#include <vector>
#include <chrono>

#include "../Generate_Associated_Legendre_Functions.h"	// Includes "../Vector_Operations.h", which we also need here."
#include "../Generate_Spherical_Harmonics.h"			// Wigner_3j
#include "../Numerical_Integration.h"



//----- Call header file with target function definitions -----
#include "Headers/AVSH_series_test_AxB_01.h"

double I_011( int L1, int L2, int L3 );


int main(){
	
	//----- Define variables -----
	int ell_max = 10;
	int n_steps_integration = 1e4;
	int w = 14;							// Fixed width for screen output.
	
	double test_r     = 0.321;				// Value of x at which to test the series.
	double test_theta = 0.432;
	
	
	std::vector<double> integrand_values_A_r             ( n_steps_integration );
	std::vector<double> integrand_values_A_1             ( n_steps_integration );
	std::vector<double> integrand_values_A_2             ( n_steps_integration );
	std::vector<double> integrand_values_B_r             ( n_steps_integration );
	std::vector<double> integrand_values_B_1             ( n_steps_integration );
	std::vector<double> integrand_values_B_2             ( n_steps_integration );
	std::vector<double> integrand_values_AxB_r          ( n_steps_integration );
	std::vector<double> integrand_values_AxB_1          ( n_steps_integration );
	std::vector<double> integrand_values_AxB_2          ( n_steps_integration );
	std::vector<double> integrand_values_AxB_r_from_AB  ( n_steps_integration );
	std::vector<double> integrand_values_AxB_12_from_AB ( n_steps_integration );
	
	std::vector<double> coeffs_A_r            ( ell_max+1 );
	std::vector<double> coeffs_A_1            ( ell_max+1 );
	std::vector<double> coeffs_A_2            ( ell_max+1 );
	std::vector<double> coeffs_B_r            ( ell_max+1 );
	std::vector<double> coeffs_B_1            ( ell_max+1 );
	std::vector<double> coeffs_B_2            ( ell_max+1 );
	std::vector<double> coeffs_AxB_r         ( ell_max+1 );
	std::vector<double> coeffs_AxB_1         ( ell_max+1 );
	std::vector<double> coeffs_AxB_2         ( ell_max+1 );
	std::vector<double> coeffs_AxB_r_from_AB ( ell_max+1 );
	std::vector<double> coeffs_AxB_1_from_AB ( ell_max+1 );
	std::vector<double> coeffs_AxB_2_from_AB ( ell_max+1 );
	
	std::vector<double> A_r_values        ( n_steps_integration );
	std::vector<double> A_theta_values    ( n_steps_integration );
	std::vector<double> A_phi_values      ( n_steps_integration );
	std::vector<double> B_r_values        ( n_steps_integration );
	std::vector<double> B_theta_values    ( n_steps_integration );
	std::vector<double> B_phi_values      ( n_steps_integration );
	std::vector<double> AxB_r_values     ( n_steps_integration );
	std::vector<double> AxB_theta_values ( n_steps_integration );
	std::vector<double> AxB_phi_values   ( n_steps_integration );
	
	

	
	//----- Start timing -----
	auto time_start = std::chrono::high_resolution_clock::now();	
	
	//----- Build lists of coordinates and precalculate vector values -----
	double delta_theta = pi / ( (double) n_steps_integration - 1.0 );
	
	std::vector<double> theta_list ( n_steps_integration );
	
	for( int i=0; i<n_steps_integration; i++ ){
		
		theta_list[i] = (double) i * delta_theta;
		
		A_r_values    [i] = A_r_function    ( test_r, theta_list[i] );
		A_theta_values[i] = A_theta_function( test_r, theta_list[i] );
		A_phi_values  [i] = A_phi_function  ( test_r, theta_list[i] );
		
		B_r_values    [i] = B_r_function    ( test_r, theta_list[i] );
		B_theta_values[i] = B_theta_function( test_r, theta_list[i] );
		B_phi_values  [i] = B_phi_function  ( test_r, theta_list[i] );
		
		AxB_r_values    [i] = A_theta_values[i]*B_phi_values  [i] - A_phi_values  [i]*B_theta_values[i];
		AxB_theta_values[i] = A_phi_values  [i]*B_r_values    [i] - A_r_values    [i]*B_phi_values  [i];
		AxB_phi_values  [i] = A_r_values    [i]*B_theta_values[i] - A_theta_values[i]*B_r_values    [i];
		
	}
	
	
	//----- Calculate VSH series coefficients -----
	
	std::cout << "ell" << "\t"<< "0" <<"\t|\t"<< std::left<<std::setw(w)<< "AxB^{r,ell}_m" <<std::left<<std::setw(w)<< "" <<"\t|\t"<< std::left<<std::setw(w)<< "AxB^{(1),ell}_m" <<std::left<<std::setw(w)<< "" <<"\t|\t"<< std::left<<std::setw(w)<< "AxB^{(2),ell}_m" <<std::left<<std::setw(w)<< "" << std::endl;
	
	for( int ell=0; ell<=ell_max; ell++ ){
		
		//--- Build lists of integrand values ---
		for( int i=0; i<n_steps_integration; i++ ){
				
			double theta = theta_list[i];
			
			double P_ell_0 = p_ell  ( cos(theta), ell );
			double P_ell_1 = p_ell_1( cos(theta), ell );
			
			integrand_values_A_r[i] = A_r_values    [i] * P_ell_0 * sin( theta );
			integrand_values_A_1[i] = A_theta_values[i] * P_ell_1 * sin( theta );
			integrand_values_A_2[i] = A_phi_values  [i] * P_ell_1 * sin( theta );
			
			integrand_values_B_r[i] = B_r_values    [i] * P_ell_0 * sin( theta );
			integrand_values_B_1[i] = B_theta_values[i] * P_ell_1 * sin( theta );
			integrand_values_B_2[i] = B_phi_values  [i] * P_ell_1 * sin( theta );
			
			integrand_values_AxB_r[i] = AxB_r_values    [i] * P_ell_0 * sin( theta );
			integrand_values_AxB_1[i] = AxB_theta_values[i] * P_ell_1 * sin( theta );
			integrand_values_AxB_2[i] = AxB_phi_values  [i] * P_ell_1 * sin( theta );
			
		}
		
		//--- Integrate and put values into array ---
		coeffs_A_r[ell] = trapezium( integrand_values_A_r, theta_list ) * sqrt( (2.0*ell+1.0)*pi );
		coeffs_A_1[ell] = trapezium( integrand_values_A_1, theta_list ) * sqrt( (2.0*ell+1.0)*pi ) / ( (double)ell*(ell+1.0) );
		coeffs_A_2[ell] = trapezium( integrand_values_A_2, theta_list ) * sqrt( (2.0*ell+1.0)*pi ) / ( (double)ell*(ell+1.0) );
		
		coeffs_B_r[ell] = trapezium( integrand_values_B_r, theta_list ) * sqrt( (2.0*ell+1.0)*pi );
		coeffs_B_1[ell] = trapezium( integrand_values_B_1, theta_list ) * sqrt( (2.0*ell+1.0)*pi ) / ( (double)ell*(ell+1.0) );
		coeffs_B_2[ell] = trapezium( integrand_values_B_2, theta_list ) * sqrt( (2.0*ell+1.0)*pi ) / ( (double)ell*(ell+1.0) );
		
		coeffs_AxB_r[ell] = trapezium( integrand_values_AxB_r, theta_list ) * sqrt( (2.0*ell+1.0)*pi );
		coeffs_AxB_1[ell] = trapezium( integrand_values_AxB_1, theta_list ) * sqrt( (2.0*ell+1.0)*pi ) / ( (double)ell*(ell+1.0) );
		coeffs_AxB_2[ell] = trapezium( integrand_values_AxB_2, theta_list ) * sqrt( (2.0*ell+1.0)*pi ) / ( (double)ell*(ell+1.0) );
		
		
		if( ell == 0 ){
			coeffs_A_1  [0] = 0;
			coeffs_A_2  [0] = 0;
			coeffs_B_1  [0] = 0;
			coeffs_B_2  [0] = 0;
			coeffs_AxB_1[0] = 0;
			coeffs_AxB_2[0] = 0;
		}
		
	}
	
	
	
	
	//----- Calculate coefficients of AxB in terms of coefficients of A and B -----
	
	
	for( int ell=0; ell<=ell_max; ell++ ){	
		
		//-- Build lists of integrand values --
		for( int L1=0; L1<=ell_max; L1++ ){
				for( int L2=0; L2<=ell_max; L2++ ){
				
				coeffs_AxB_r_from_AB[ell] += ( coeffs_A_1[L1]*coeffs_B_2[L2] - coeffs_A_2[L2]*coeffs_B_1[L1] ) * I_011( ell, L1, L2  );
				coeffs_AxB_1_from_AB[ell] += ( coeffs_A_2[L2]*coeffs_B_r[L1] - coeffs_A_r[L1]*coeffs_B_2[L2] ) * I_011( L1 , L2, ell );
				coeffs_AxB_2_from_AB[ell] += ( coeffs_A_r[L1]*coeffs_B_1[L2] - coeffs_A_1[L2]*coeffs_B_r[L1] ) * I_011( L1 , L2, ell );
				
			}
		}
		
		coeffs_AxB_1_from_AB[ell] /= ( (double) ell*(ell+1) );
		coeffs_AxB_2_from_AB[ell] /= ( (double) ell*(ell+1) );
		
		if( ell == 0 ){
			coeffs_AxB_1_from_AB[ell] = 0;
			coeffs_AxB_2_from_AB[ell] = 0;
		}
		
		
		
		std::cout << ell << "\t"<< 0 <<"\t|\t"<< std::left<<std::setw(w)<< coeffs_AxB_r[ell] <<std::left<<std::setw(w)<< coeffs_AxB_r_from_AB[ell]<<"\t|\t"<< std::left<<std::setw(w)<< coeffs_AxB_1[ell] <<std::left<<std::setw(w)<< coeffs_AxB_1_from_AB[ell]<<"\t|\t"<< std::left<<std::setw(w)<< coeffs_AxB_2[ell] <<std::left<<std::setw(w)<< coeffs_AxB_2_from_AB[ell] << std::endl;
		//std::cout << ell << "\t"<< 0 <<"\t|\t"<< std::left<<std::setw(w)<< coeffs_A_r[ell] <<std::left<<std::setw(w)<< coeffs_A_1[ell] << std::left<<std::setw(w)<< coeffs_A_2[ell] << "(coeffs A for checking)" << std::endl;
		//std::cout << ell << "\t"<< 0 <<"\t|\t"<< std::left<<std::setw(w)<< coeffs_B_r[ell] <<std::left<<std::setw(w)<< coeffs_B_1[ell] << std::left<<std::setw(w)<< coeffs_B_2[ell] << "(coeffs B for checking)" << std::endl;
	}
	
	
	
	
	//----- Evaluate VSH series at the chosen gridpoint -----
	
	std::vector<double> A_exact  { A_r_function(test_r,test_theta), A_theta_function(test_r,test_theta), A_phi_function(test_r,test_theta) };
	std::vector<double> B_exact  { B_r_function(test_r,test_theta), B_theta_function(test_r,test_theta), B_phi_function(test_r,test_theta) };
	std::vector<double> AxB_exact = cross_product( A_exact, B_exact );
	std::vector<double> AxB_expression { AxB_r_function(test_r,test_theta), AxB_theta_function(test_r,test_theta), AxB_phi_function(test_r,test_theta) };
	
	std::vector<double> A_series           { 0, 0, 0 };
	std::vector<double> B_series           { 0, 0, 0 };
	std::vector<double> AxB_series         { 0, 0, 0 };
	std::vector<double> AxB_series_from_AB { 0, 0, 0 };
	
	for( int ell=0; ell<=ell_max; ell++ ){
		
		double sqrt_factor = sqrt( (2.0*ell+1.0) / ( 4.0 * pi ) );
		
		double P_ell_0 = p_ell  ( cos(test_theta), ell );
		double P_ell_1 = p_ell_1( cos(test_theta), ell );
		
		A_series[0] += sqrt_factor * P_ell_0 * coeffs_A_r[ell];
		A_series[1] += sqrt_factor * P_ell_1 * coeffs_A_1[ell];
		A_series[2] += sqrt_factor * P_ell_1 * coeffs_A_2[ell];
		
		B_series[0] += sqrt_factor * P_ell_0 * coeffs_B_r[ell];
		B_series[1] += sqrt_factor * P_ell_1 * coeffs_B_1[ell];
		B_series[2] += sqrt_factor * P_ell_1 * coeffs_B_2[ell];
		
		AxB_series[0] += sqrt_factor * P_ell_0 * coeffs_AxB_r[ell];
		AxB_series[1] += sqrt_factor * P_ell_1 * coeffs_AxB_1[ell];
		AxB_series[2] += sqrt_factor * P_ell_1 * coeffs_AxB_2[ell];
		
		AxB_series_from_AB[0] += sqrt_factor * P_ell_0 * coeffs_AxB_r_from_AB[ell];
		AxB_series_from_AB[1] += sqrt_factor * P_ell_1 * coeffs_AxB_1_from_AB[ell];
		AxB_series_from_AB[2] += sqrt_factor * P_ell_1 * coeffs_AxB_2_from_AB[ell];
		
	}
	
	std::cout << std::endl;
	cout_vector_setw( A_exact , w, "A_exact " );
	cout_vector_setw( A_series, w, "A_series" );
	std::cout << std::endl;
	cout_vector_setw( B_exact , w, "B_exact " );
	cout_vector_setw( B_series, w, "B_series" );
	std::cout << std::endl;
	cout_vector_setw( AxB_exact         , w, "AxB_exact         " );
	cout_vector_setw( AxB_expression    , w, "AxB_expression    " );
	cout_vector_setw( AxB_series        , w, "AxB_series        " );
	cout_vector_setw( AxB_series_from_AB, w, "AxB_series_from_AB" );
	
	
	
	
	//----- Output execution time and finish code -----
	auto   time_stop   = std::chrono::high_resolution_clock::now();
	double exec_time   = std::chrono::duration_cast<std::chrono::nanoseconds>( time_stop - time_start ).count() * 1e-9;
	int    exec_time_h = exec_time / 3600;
	int    exec_time_m = exec_time / 60 - exec_time_h * 60;
	int    exec_time_s = exec_time - exec_time_h * 3600 - exec_time_m * 60;
	std::cout << "\nExecution time (h:mm:ss):\t" << exec_time_h;
	if( exec_time_m < 10 ){
		std::cout << ":0" << exec_time_m;
	} else {
		std::cout << ":"  << exec_time_m;
	}
	if( exec_time_s < 10 ){
		std::cout << ":0" << exec_time_s << std::endl;
	} else {
		std::cout << ":"  << exec_time_s << std::endl;
	}
	
	
	std::cout << "done" << std::endl;
	return 0;
}




double I_011( int L1, int L2, int L3 ){
	// 1/(4sqrt(pi)) sqrt(2L1+1) sqrt(L2(L2+1)(2L2+1)) sqrt(L3(L3+1)(2L3+1)) int_0^pi P_{L1}^0 P_{L2}^1 P_{L3}^1 sin(theta) dtheta
	// ../Check integrals/20230807_P_L1_0_P_L2_1_P_L_1.cpp
	
	if( L3 > L1+L2 ){
		return 0;
	}
	
	double ret = 0;
	
	for( int L12=std::max(abs(L1-L2),1); L12<=L1+L2; L12++ ){
		if( (L1+L2+L12)%2 == 0 ){
		
			double G12 = - ( 2*L12+1 ) * wigner_3j_000( L1, L2, L12 ) * wigner_3j( L1, L2, L12, 0, 1, -1 );
			
			if( (L12+L3)%2 == 0 ){
				for( int L123=std::max(abs(L12-L3),2); L123<=L12+L3; L123++ ){
					if( L123 %2 == 0 ){
						double G123          = ( 2*L123+1 ) * wigner_3j_000( L12, L3, L123 ) * wigner_3j( L12, L3, L123, 1, 1, -2 );
						double sqrt_fraction = pow( (L123+2) * (L123+1) * L123 * (L123-1), -0.5 );
						
						ret += G12 * G123 * sqrt_fraction;
					}
				}
			}
			
		}
	}
	
	return ret * sqrt( (2*L1+1) * L2*(L2+1)*(2*L2+1) * L3*(L3+1)*(2*L3+1) / pi );
}