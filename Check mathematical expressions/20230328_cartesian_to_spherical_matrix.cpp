/*
cd OneDrive\PhD\Codes\Check mathematical expressions
g++ 20230328_cartesian_to_spherical_matrix.cpp -o 20230328_cartesian_to_spherical_matrix
20230328_cartesian_to_spherical_matrix

Check that the inverse matrix going from Cartesian to spherical polar coordinates is correct.
Book 14, p60.
*/

#include <iostream>
#include <iomanip>
#include <math.h>
#include "../Matrix_Operations.h"


//----- Function declarations -----
const double pi = acos( -1 );


int main(){
	
	//----- Choose the coordinates -----
	double t = 0.123 * pi;	// theta
	double p = 0.432 * pi;	// phi
	
	double A_x = 0.321;	// Cartesian vector components.
	double A_y = 0.543;
	double A_z = 0.345;
	
	//----- Define matrix -----
	std::vector< std::vector<double> > A ( 3, std::vector<double> ( 3 ) );
	// The matrix calculating Cartesian components from spherical components. This is known to be correct.
	A[0] = {  sin(t)*cos(p),  cos(t)*cos(p), -sin(p) };
	A[1] = {  sin(t)*sin(p),  cos(t)*sin(p),  cos(p) };
	A[2] = {  cos(t)       , -sin(t)       ,  0      };
	
	cout_matrix_3x3( A );
	
	
	double det_A_guess = 1;
	
	
	std::vector< std::vector<double> > A_inverse_guess ( 3, std::vector<double> ( 3 ) );
	A_inverse_guess[0] = {  sin(t)*cos(p),  sin(t)*sin(p),  cos(t) };
	A_inverse_guess[1] = {  cos(t)*cos(p),  cos(t)*sin(p), -sin(t) };
	A_inverse_guess[2] = {  -sin(p)      ,  cos(p)       ,  0      };
	
	
	//----- Determinant -----
	std::cout << "det(A):\t" << det_3x3( A ) <<"\t"<< det_A_guess << std::endl;
	
	
	//----- Inverse -----
	std::vector< std::vector<double> > A_inverse_calc = inverse_3x3( A );
	cout_matrix_3x3( A_inverse_calc );
	cout_matrix_3x3( A_inverse_guess );
	
	
	//----- Calculate spherical polar vector components -----
	// Now that we know A^{-1} is correct, hardcode the transformation to spherical polar coordinates.
	
	double A_r = sin(t)*cos(p)*A_x + sin(t)*sin(p)*A_y + cos(t)*A_z;
	double A_t = cos(t)*cos(p)*A_x + cos(t)*sin(p)*A_y - sin(t)*A_z;
	double A_p =       -sin(p)*A_x +        cos(p)*A_y             ;
	
	int w = 15;
	std::cout << "Spherical polar components:\t" <<std::left<<std::setw(w)<< A_r <<std::left<<std::setw(w)<< A_t <<std::left<<std::setw(w)<< A_p << std::endl;
	
	double A_x2 = sin(t)*cos(p)*A_r + cos(t)*cos(p)*A_t - sin(p)*A_p;
	double A_y2 = sin(t)*sin(p)*A_r + cos(t)*sin(p)*A_t + cos(p)*A_p;
	double A_z2 = cos(t)       *A_r - sin(t)       *A_t             ;
	
	std::cout << "Recalc. Cart. components:\t" <<std::left<<std::setw(w)<< A_x2 <<std::left<<std::setw(w)<< A_y2 <<std::left<<std::setw(w)<< A_z2 << std::endl;
	
	
	
	return 0;
}