/*
cd OneDrive\PhD\Codes\Check vector expressions
g++ Check_parallel_and_perpendicular_decomposition.cpp -o Check_parallel_and_perpendicular_decomposition --std=c++11
Check_parallel_and_perpendicular_decomposition

Given two vectors A and B, check the decomposition of A into components parallel and perpendicular to B.

A_para = ( A dot B ) / B^2 * B
A_perp = - 1/B^2 * [ ( A x B ) x B ]
A_para + A_perp = A
*/

#include <iostream>
#include <iomanip>
#include <math.h>
#include <vector>
#include <chrono>

#include "../Generate_Spherical_Harmonics.h"
#include "../Generate_Associated_Legendre_Functions.h"


int main(){
	
	//----- Define variables -----
	
	double Ax = 1.234;
	double Ay = 2.345;
	double Az = 4.321;
	double Bx = 5.321;
	double By = 6.324;
	double Bz = 0.432;
	
	std::vector<double> A { Ax, Ay, Az };
	std::vector<double> B { Bx, By, Bz };
	
	double A_mag = vector_magnitude( A );
	double B_mag = vector_magnitude( B );
	
	double one_over_Bsq = pow( vector_magnitude( B ), -2 );
	
	
	//----- Parallel component -----
	double para_factor = dot_product( A, B ) * one_over_Bsq;
	std::vector<double> A_para = scalar_times_vector( para_factor, B );
	double A_para_mag = vector_magnitude( A_para );
	
	double A_para_dot_B = dot_product( A_para, B );
	std::vector<double> A_para_cross_B = cross_product( A_para, B );
	
	std::cout << "A_para dot B:\t" << A_para_dot_B <<"\t"<< A_para_mag * B_mag << std::endl;
	cout_vector( A_para_cross_B, "A_para x B" );
	
	
	//----- Perpendicular component -----
	double perp_factor = - one_over_Bsq;
	std::vector<double> A_cross_B = cross_product( A, B );
	std::vector<double> A_cross_B_cross_B = cross_product( A_cross_B, B );
	std::vector<double> A_perp = scalar_times_vector( perp_factor, A_cross_B_cross_B );
	double A_perp_mag = vector_magnitude( A_perp );
	
	double A_perp_dot_B = dot_product( A_perp, B );
	std::vector<double> A_perp_cross_B = cross_product( A_perp, B );
	
	std::cout << "\nA_perp dot B:\t" << A_perp_dot_B << std::endl;
	cout_vector( A_perp_cross_B, "A_perp x B" );
	
	
	//----- Vector sum -----
	std::vector<double> A_para_plus_A_perp = vector_sum( A_para, A_perp );
	
	cout_vector( A_para_plus_A_perp, "\nA_para + A_perp" );
	cout_vector( A, "A\t\t" );
	
	
	
	
	std::cout << "\ndone" << std::endl;
	
	
	return 0;
	
}