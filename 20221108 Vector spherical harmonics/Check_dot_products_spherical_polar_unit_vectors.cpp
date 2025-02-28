// cd OneDrive\PhD\Codes\20221108 Vector spherical harmonics
// g++ Check_dot_products_spherical_polar_unit_vectors.cpp -o Check_dot_products_spherical_polar_unit_vectors

/*
Check that e_i cdot e_j = delta_{ij}.
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
	
	
	//----- Check cross product works (it was defined at the same time this code was written) -----
	std::vector<double> a_test = { 2,  3,  5 };
	std::vector<double> b_test = { 7, 11, 13 };
	double test = dot_product_cart( a_test, b_test );
	std::cout << "a dot b:\t" << test << std::endl;		// Correct result is 112
	
	
	//----- Define spherical polar unit vectors in terms of i,j,k -----
	std::vector<double> e_r = { sin(t)*cos(p), sin(t)*sin(p),  cos(t) };
	std::vector<double> e_t = { cos(t)*cos(p), cos(t)*sin(p), -sin(t) };
	std::vector<double> e_p = {       -sin(p),        cos(p),      0  };
	
	print_vector( e_r, "e_r" );
	print_vector( e_t, "e_t" );
	print_vector( e_p, "e_p", "\n" );
	
	
	//----- Cycle through the combinations of vectors -----
	std::cout << "e_r dot e_r\t" << dot_product_cart( e_r, e_r )         << std::endl;
	std::cout << "e_r dot e_t\t" << dot_product_cart( e_r, e_t )         << std::endl;
	std::cout << "e_r dot e_p\t" << dot_product_cart( e_r, e_p ) << "\n" << std::endl;
	
	std::cout << "e_t dot e_r\t" << dot_product_cart( e_t, e_r )         << std::endl;
	std::cout << "e_t dot e_t\t" << dot_product_cart( e_t, e_t )         << std::endl;
	std::cout << "e_t dot e_p\t" << dot_product_cart( e_t, e_p ) << "\n" << std::endl;
	
	std::cout << "e_p dot e_r\t" << dot_product_cart( e_p, e_r )         << std::endl;
	std::cout << "e_p dot e_t\t" << dot_product_cart( e_p, e_t )         << std::endl;
	std::cout << "e_p dot e_p\t" << dot_product_cart( e_p, e_p ) << "\n" << std::endl;	
	
	return 0;
}