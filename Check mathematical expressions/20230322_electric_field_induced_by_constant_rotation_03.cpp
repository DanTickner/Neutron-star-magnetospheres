/*
cd OneDrive\PhD\Codes\Check mathematical expressions
g++ 20230322_electric_field_induced_by_constant_rotation_03.cpp -o 20230322_electric_field_induced_by_constant_rotation_03
20230322_electric_field_induced_by_constant_rotation_03

Let an object with magnetic field B = B_r e_r + B_theta e_theta + B_phi e_phi
be rotating with angular velocity Omega = Omega_r e_r + Omega_theta e_theta + Omega_phi e_phi.
Then, it induces an electric field E = E_r e_r + E_theta e_theta + e_phi e_phi
given by
E = -( Omega cross r ) cross B
This code tests whether the RHS of that expression becomes
  [ ( Omega_r theta - Omega_theta r ) B_theta + ( Omega_r phi - Omega_phi r ) B_phi ] e_r
+ [ ( Omega_theta phi - Omega_phi theta ) B_phi + ( Omega_theta r - Omega_r theta ) B_r ] e_theta
+ [ ( Omega_phi r - Omega_r phi ) B_r + ( Omega_phi theta - Omega_theta phi ) B_theta ] e_phi

V02: Previous expressions for Omega were incorrect. We actually have
     vector Omega = Omega e_z = Omega cos(t) e_r - Omega sin(t) e_t
	 where e_z = cos(t) e_r - sin(t) e_t.
	 Also remove generic expressions for vector components and input the exact expressions.

*/


#include <iostream>
#include <iomanip>
#include <math.h>
#include <vector>
#include "../Vector_Operations.h"


int main(){
	
	
	//----- Define variables -----
	double x = 0.567;
	double y = 0.543;
	double z = 0.765;
	
	double B_0 = 0.456;
	double Omega = 0.432;
	
	
	
	//----- Build vectors -----
	double r = sqrt( x*x + y*y + z*z );
	std::vector<double> vector_Omega { 0, 0, Omega };
	std::vector<double> vector_r { x, y, z };
	std::vector<double> vector_B { B_0*pow(r,-5)*3.0*x*z, B_0*pow(r,-5)*3.0*y*z, B_0*pow(r,-5)*(3.0*z*z-r*r) };
	
	
	//----- Calculate correct expression -----
	std::vector<double> Omega_cross_r = cross_product( vector_Omega, vector_r );
	std::vector<double> Omega_cross_r_all_cross_B = cross_product( Omega_cross_r, vector_B );
	std::vector<double> E_LHS { -Omega_cross_r_all_cross_B[0], -Omega_cross_r_all_cross_B[1], -Omega_cross_r_all_cross_B[2] };
	
	
	//----- Calculate guessed expression -----
	std::vector<double> Omega_cross_r_guess { -Omega*y, Omega*x, 0 };
	
	double E_RHS_0 = B_0 * Omega ;
	
	//std::vector<double> E_RHS { E_RHS_0*pow(r,-5)*(-3.0*x*z*z+x*r*r), E_RHS_0*pow(r,-5)*(-3.0*y*z*z+y*r*r), E_RHS_0*pow(r,-5)*(3.0*y*y*z+3.0*x*x*z) };
	
	std::vector<double> E_RHS {
		                        E_RHS_0 * ( x*pow(r,-3) - 3.0*x*z*z*pow(r,-5) ),
								E_RHS_0 * ( y*pow(r,-3) - 3.0*y*z*z*pow(r,-5) ),
								//E_RHS_0 * ( 3.0*z*(x*x+y*y)*pow(r,-5) )
								E_RHS_0 * ( 3.0*z*pow(r,-3) - 3.0*pow(z,3)*pow(r,-5) )
	};
	
	/*std::vector<double> E_RHS {
		                        E_RHS_0*pow(r,-3) * x * ( 1.0 - 3.0*z*z*pow(r,-2) ),
								E_RHS_0*pow(r,-3) * y * ( 1.0 - 3.0*z*z*pow(r,-2) ),
								E_RHS_0*pow(r,-3) * 3.0*z * ( 1.0 - z*z*pow(r,-2) )
	};*/
	
	
	//----- Output results -----
	cout_vector( Omega_cross_r      , "Omega x r correct" );
	cout_vector( Omega_cross_r_guess, "Omega x r guess  " );
	std::cout << std::endl;
	
	cout_vector( E_LHS, "Correct" );
	cout_vector( E_RHS, "Guessed" );
	
	return 0;
}