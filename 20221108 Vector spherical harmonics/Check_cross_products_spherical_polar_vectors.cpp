// cd OneDrive\PhD\Codes\20221108 Vector spherical harmonics
// g++ Check_cross_products_spherical_polar_vectors.cpp -o Check_cross_products_spherical_polar_vectors

/*
Test whether A x B = ( A_tB_p -A_pB_t ) e_r + ( A_pB_r - A_rB_p ) e_t + ( A_rB_t - A_tB_r ) e_p
*/

#include <iostream>
#include <vector>
#include <math.h>
#include "../Vector_Operations.h"						// print_vector()

int main(){
	
	//----- Define variables -----
	const double pi = 3.14159265358979323846;
	
	double t = 0.7 * pi;	// theta
	double p = 1.2 * pi;	// phi
	
	std::vector<double> A_spher = { 0.5, 0.2*pi, 1.1*pi };
	std::vector<double> B_spher = { 0.3, 0.7*pi, 1.3*pi };
	
	std::vector<double> e_r_spher = { 1, 0, 0 };
	std::vector<double> e_t_spher = { 0, 1, 0 };
	std::vector<double> e_p_spher = { 0, 0, 1 };
	std::vector<double> e_r_cart  = { sin(t)*cos(p), sin(t)*sin(p),  cos(t) };
	std::vector<double> e_t_cart  = { cos(t)*cos(p), cos(t)*sin(p), -sin(t) };
	std::vector<double> e_p_cart  = {       -sin(p),        cos(p),      0  };
	
	
	//----- Calculate cross product -----
	std::vector<double> A_cross_B_spher_test = { A_spher[1]*B_spher[2]-A_spher[2]*B_spher[1], A_spher[2]*B_spher[0]-A_spher[0]*B_spher[2], A_spher[0]*B_spher[1]-A_spher[1]*B_spher[0] };
	
	print_vector( A_cross_B_spher_test, "Spher. test\t" );
	
	std::vector<double> A_cart = spher_to_cart( A_spher );
	std::vector<double> B_cart = spher_to_cart( B_spher );
	std::vector<double> A_cross_B_cart = cross_product_cart( A_cart, B_cart );
	std::vector<double> A_cross_B_spher_fromcart = cart_to_spher( A_cross_B_cart );
	
	print_vector( A_cross_B_spher_fromcart, "Spher. from cart" );
	
	
	//----- Calculate cross product with e_r-----
	std::vector<double> A_cross_e_r_spher_test = { A_spher[1]*e_r_spher[2]-A_spher[2]*e_r_spher[1], A_spher[2]*e_r_spher[0]-A_spher[0]*e_r_spher[2], A_spher[0]*e_r_spher[1]-A_spher[1]*e_r_spher[0] };
	
	print_vector( A_cross_B_spher_test, "\nA x e_r, Spher. test\t" );
	
	std::vector<double> e_r_cross_A_cart = cross_product_cart( e_r_cart, A_cart );
	std::vector<double> e_r_cross_A_spher_fromcart = cart_to_spher( e_r_cross_A_cart );
	
	print_vector( e_r_cross_A_spher_fromcart, "e_r x A, Spher. from cart" );
	
	return 0;
}