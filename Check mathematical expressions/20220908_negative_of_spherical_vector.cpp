// cd OneDrive\PhD\Codes\Test mathematical expressions
// g++ 20220908_negative_of_spherical_vector.cpp -o 20220908_negative_of_spherical_vector

// -x^{2n+2} - (2n+2)! sum_{k=0}^n (-1)^k/(2k)! x^{2k} ?=? (2n+2)! sum_{k=0}^{n+1} (-1)^k/(2k)! x^{2k}

#include <iostream>
#include <math.h>

int main(){
	
	//----- Define variables -----
	const double pi = 3.14159265358979323846;
	std::cout.precision( 3 );
	
	double A_r     = 0.5;
	double A_theta = 0.2 * pi;
	double A_phi   = 1.2 * pi;
	
	
	//----- Calculate A in Cartesian coordinates -----
	double A_x = A_r * sin(A_theta) * cos(A_phi );
	double A_y = A_r * sin(A_theta) * sin(A_phi );
	double A_z = A_r * cos(A_theta);
	
	//----- Use this to calculate mA=-A in Cartesian coordinates -----
	double mA_x = - A_x;
	double mA_y = - A_y;
	double mA_z = - A_z;
	
	//----- Use this to calculate mA in spherical coordinates -----
	double mA_r     = sqrt ( mA_x*mA_x + mA_y*mA_y + mA_z*mA_z );
	double mA_theta = acos ( mA_z / mA_r );
	double mA_phi   = atan2( mA_y, mA_x );
	
	//----- Guessed mA in spherical coordinates -----
	double mA_r_guess     = A_r;
	double mA_theta_guess = pi - A_theta;
	double mA_phi_guess   = A_phi - pi;
	/*if( A_y >= 0 ){
		mA_phi_guess = A_phi - pi;
	}
	else {
		mA_phi_guess = A_phi + pi;
	}*/
	
	
	//----- Output results -----
	std::cout << " A in spherical coordinates:" <<"\t"<<  A_r  <<"\t"<<  A_theta <<"\t"<<  A_phi << std::endl;
	std::cout << " A in Cartesian coordinates:" <<"\t"<<  A_x  <<"\t"<<  A_y     <<"\t"<<  A_z   << std::endl;
	std::cout << "-A in Cartesian coordinates:" <<"\t"<< mA_x  <<"\t"<< mA_y     <<"\t"<< mA_z   << std::endl;
	std::cout << "-A in spherical coordinates:" <<"\t"<< mA_r  <<"\t"<< mA_theta <<"\t"<< mA_phi << std::endl;
	
	std::cout << "\nGuessed -A in spherical coordinates:" <<"\t"<< mA_r_guess  <<"\t"<< mA_theta_guess <<"\t"<< mA_phi_guess << std::endl;
	
	return 0;
}