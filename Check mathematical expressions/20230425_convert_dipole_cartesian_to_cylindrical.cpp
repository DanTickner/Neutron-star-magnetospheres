/*
cd OneDrive\PhD\Codes\Check mathematical expressions
g++ 20230425_convert_dipole_cartesian_to_cylindrical.cpp -o 20230425_convert_dipole_cartesian_to_cylindrical
20230425_convert_dipole_cartesian_to_cylindrical


B = B_0 1/r^5 ( 3xz ihat - 3yz jhat + (3z^2-r^2) khat )
  = B_0 1/(varpi^2+z^2)^3/2 ( 3 varpi z cos(2phi) e_varpi - 3 varpi z sin(2phi) e_phi + (2z^2-varpi^2) e_z )

Book 16, p19.
*/

#include <iostream>
#include <iomanip>
#include <math.h>
#include <vector>

#include "../Vector_Operations.h"	// cout_vector()


int main(){
	
	//----- Define variables -----
	double B_0 = 0.987;
	double x = 0.234;
	double y = 0.345;
	double z = 0.456;
	
	int w = 12;				// Fixed width for text output.
	
	
	//----- Calculated variables -----
	double r = sqrt( x*x + y*y + z*z );
	
	double varpi = sqrt( x*x + y*y );
	double phi = atan( y / x );
	
	
	//----- Calculate the translation of the Cartesian unit vectors to cylindrical unit vectors -----
	std::vector<double> ihat_to_cyl { cos(phi), -sin(phi), 0   };
	std::vector<double> jhat_to_cyl { sin(phi),  cos(phi), 0   };
	std::vector<double> khat_to_cyl {        0,         0, 1.0 };
	
	
	//----- Calculate B in cartesian coordinates -----
	double B_cart_0 = B_0 * pow(r,-5);
	
	std::vector<double> B_cart {
		B_cart_0 *  3.0*x*z,
		B_cart_0 * -3.0*y*z,
		B_cart_0 * ( 3.0*z*z - r*r )
	};
	
	
	//----- Calculate B in cylindrical coordinates -----
	std::vector<double> B_cyl_from_cart {
		B_cart[0]*ihat_to_cyl[0] + B_cart[1]*jhat_to_cyl[0] + B_cart[2]*khat_to_cyl[0],
		B_cart[0]*ihat_to_cyl[1] + B_cart[1]*jhat_to_cyl[1] + B_cart[2]*khat_to_cyl[1],
		B_cart[0]*ihat_to_cyl[2] + B_cart[1]*jhat_to_cyl[2] + B_cart[2]*khat_to_cyl[2]
	};
	
	double B_cyl_0 = B_0 * pow(varpi*varpi+z*z,-2.5);
	
	std::vector<double> B_cyl_guess {
		
		/*
		//----- V1 -----
		B_cyl_0 * ( 3.0*varpi*cos(phi)*z* cos(phi) - 3.0*varpi*sin(phi)*z*sin(phi) ),
		B_cyl_0 * ( 3.0*varpi*cos(phi)*z*-sin(phi) - 3.0*varpi*sin(phi)*z*cos(phi) ),
		B_cyl_0 * ( 3.0*z*z - ( varpi*varpi+z*z ) )
		*/
		
		/*
		//----- V2 -----
		B_cyl_0 * 3.0*varpi*z * ( cos(phi)* cos(phi) - sin(phi)*sin(phi) ),
		B_cyl_0 * 3.0*varpi*z * ( cos(phi)*-sin(phi) - sin(phi)*cos(phi) ),
		B_cyl_0 * ( 2.0*z*z - varpi*varpi )
		*/
		
		//----- V3 -----
		B_cyl_0 *  3.0*varpi*z * cos(2.0*phi),
		B_cyl_0 * -3.0*varpi*z * sin(2.0*phi),
		B_cyl_0 * ( 2.0*z*z - varpi*varpi )
	};
	
	std::cout << "Check r:\t" << B_cart_0 <<"\t"<< B_cyl_0 << std::endl;
	
	//----- Output results -----
	cout_vector_setw( B_cart         , w, "B_cart         " );
	
	cout_vector_setw( B_cyl_from_cart, w, "B_cyl_from_cart" );
	cout_vector_setw( B_cyl_guess    , w, "B_cyl_guess    " );
	
	
	std::cout << "\nDone" << std::endl;
	
	return 0;
}