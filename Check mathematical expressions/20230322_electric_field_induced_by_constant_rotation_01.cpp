/*
cd OneDrive\PhD\Codes\Test mathematical expressions
g++ 20230322_electric_field_induced_by_constant_rotation_01.cpp -o 20230322_electric_field_induced_by_constant_rotation_01
20230322_electric_field_induced_by_constant_rotation_01

Let an object with magnetic field B = B_r e_r + B_theta e_theta + B_phi e_phi
be rotating with angular velocity Omega = Omega_r e_r + Omega_theta e_theta + Omega_phi e_phi.
Then, it induces an electric field E = E_r e_r + E_theta e_theta + e_phi e_phi
given by
E = -( Omega cross r ) cross B
This code tests whether the RHS of that expression becomes
  [ ( Omega_r theta - Omega_theta r ) B_theta + ( Omega_r phi - Omega_phi r ) B_phi ] e_r
+ [ ( Omega_theta phi - Omega_phi theta ) B_phi + ( Omega_theta r - Omega_r theta ) B_r ] e_theta
+ [ ( Omega_phi r - Omega_r phi ) B_r + ( Omega_phi theta - Omega_theta phi ) B_theta ] e_phi
*/


#include <iostream>
#include <iomanip>
#include <math.h>
#include <vector>
#include "../Vector_Operations.h"


int main(){
	
	
	//----- Define the vectors -----
	double Omega_r = 0.123;
	double Omega_t = 0.321;
	double Omega_p = 0.432;
	double r = 0.567;
	double t = 0.765;
	double p = 0.654;
	double B_r = 0.445;
	double B_t = 0.543;
	double B_p = 0.234;
	
	
	
	//----- Build vectors -----
	std::vector<double> vector_Omega { Omega_r, Omega_t, Omega_p };
	std::vector<double> vector_r { r, t, p };
	std::vector<double> vector_B { B_r, B_t, B_p };
	
	
	//----- Calculate correct expression -----
	std::vector<double> Omega_cross_r = cross_product( vector_Omega, vector_r );
	std::vector<double> Omega_cross_r_all_cross_B = cross_product( Omega_cross_r, vector_B );
	std::vector<double> E_LHS { -Omega_cross_r_all_cross_B[0], -Omega_cross_r_all_cross_B[1], -Omega_cross_r_all_cross_B[2] };
	
	
	//----- Calculate guessed expression -----
	double E_RHS_r = ( Omega_r * t - Omega_t * r ) * B_t + ( Omega_r * p - Omega_p * r ) * B_p;
	double E_RHS_t = ( Omega_t * p - Omega_p * t ) * B_p + ( Omega_t * r - Omega_r * t ) * B_r;
	double E_RHS_p = ( Omega_p * r - Omega_r * p ) * B_r + ( Omega_p * t - Omega_t * p ) * B_t;
	
	std::vector<double> E_RHS { E_RHS_r, E_RHS_t, E_RHS_p };
	
	
	//----- Output results -----
	cout_vector( E_LHS, "Correct" );
	cout_vector( E_RHS, "Guessed" );
	
	
	
	//----- Special case 1: Omega is azimuthal (ignore expressions for Omega_r and Omega_t) -----
	std::vector<double> vector_Omega_sc1 { 0, 0, Omega_p };
	
	std::vector<double> Omega_cross_r_sc1 = cross_product( vector_Omega_sc1, vector_r );
	std::vector<double> Omega_cross_r_all_cross_B_sc1 = cross_product( Omega_cross_r_sc1, vector_B );
	std::vector<double> E_LHS_sc1 { -Omega_cross_r_all_cross_B_sc1[0], -Omega_cross_r_all_cross_B_sc1[1], -Omega_cross_r_all_cross_B_sc1[2] };
	
	double E_RHS_sc1_r = Omega_p * -r * B_p;
	double E_RHS_sc1_t = Omega_p * -t * B_p;
	double E_RHS_sc1_p = Omega_p * ( r * B_r + t * B_t );
	std::vector<double> E_RHS_sc1 { E_RHS_sc1_r, E_RHS_sc1_t, E_RHS_sc1_p };
	
	std::cout << std::endl;
	cout_vector( E_LHS_sc1, "Correct SC1" );
	cout_vector( E_RHS_sc1, "Guessed SC1" );
	
	
	
	
	//----- Special case 2: Omega is azimuthal and B has no azimuthal component -----
	std::vector<double> vector_B_sc2 { B_r, B_t, 0 };
	
	std::vector<double> Omega_cross_r_all_cross_B_sc2 = cross_product( Omega_cross_r_sc1, vector_B_sc2 );
	std::vector<double> E_LHS_sc2 { -Omega_cross_r_all_cross_B_sc2[0], -Omega_cross_r_all_cross_B_sc2[1], -Omega_cross_r_all_cross_B_sc2[2] };
	
	std::vector<double> E_RHS_sc2 { 0, 0, E_RHS_sc1_p };
	
	std::cout << std::endl;
	cout_vector( E_LHS_sc2, "Correct SC2" );
	cout_vector( E_RHS_sc2, "Guessed SC2" );
	
	return 0;
}