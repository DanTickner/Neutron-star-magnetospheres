// cd OneDrive\PhD\Codes\20221108 Vector spherical harmonics
// g++ Check_cross_products_spherical_polar_unit_vectors.cpp -o Check_cross_products_spherical_polar_unit_vectors

/*
Check that e_i x e_j = delta_{ij} epsilon_{ijk} with i,j in { r, theta, phi }.
*/

#include <iostream>
#include <vector>
#include <math.h>
#include "../Vector_Operations.h"						// print_vector()

int main(){
	
	//----- Define variables -----
	const double pi = 3.14159265358979323846;
	
	double t = 0;//0.7 * pi;	// theta
	double p = 0;//1.2 * pi;	// phi
	
	
	/*//----- Check cross product works (it was defined at the same time this code was written) -----
	std::vector<double> a_test = { 2,  3,  5 };
	std::vector<double> b_test = { 7, 11, 13 };
	std::vector<double> test = cross_product_cart( a_test, b_test );
	std::cout << "a x b:\t" << test[0] <<"\t"<< test[1] <<"\t"<< test[2] << std::endl;	// Correct result is ( -16, 9, 1 ).*/
	
	
	//----- Define spherical polar unit vectors in terms of i,j,k -----
	std::vector<double> e_r = { sin(t)*cos(p), sin(t)*sin(p),  cos(t) };
	std::vector<double> e_t = { cos(t)*cos(p), cos(t)*sin(p), -sin(t) };
	std::vector<double> e_p = {       -sin(p),        cos(p),      0  };
	
	print_vector( e_r, "e_r" );
	print_vector( e_t, "e_t" );
	print_vector( e_p, "e_p", "\n" );
	
	
	//----- Cycle through the combinations of vectors -----
	print_vector( cross_product_cart( e_r, e_r ), "e_r x e_r" );
	print_vector( cross_product_cart( e_r, e_t ), "e_r x e_t" );
	print_vector( cross_product_cart( e_r, e_p ), "e_r x e_p", "\n" );
	
	print_vector( cross_product_cart( e_t, e_r ), "e_t x e_r" );
	print_vector( cross_product_cart( e_t, e_t ), "e_t x e_t" );
	print_vector( cross_product_cart( e_t, e_p ), "e_t x e_p", "\n" );
	
	print_vector( cross_product_cart( e_p, e_r ), "e_p x e_r" );
	print_vector( cross_product_cart( e_p, e_t ), "e_p x e_t" );
	print_vector( cross_product_cart( e_p, e_p ), "e_p x e_p", "\n" );
	
	
	return 0;
}