// cd OneDrive\PhD\Codes\Test mathematical expressions
// g++ 20220908_atan2_sinphi_cosphi.cpp -o 20220908_atan2_sinphi_cosphi

// atan2( sin(phi), cos(phi) ) = ?

#include <iostream>
#include <math.h>

int main(){
	
	//----- Define variables -----
	const double pi = 3.14159265358979323846;
	
	double phi_1 = 0;
	double phi_2 = pi * 0.3;
	double phi_3 = pi * 0.5;
	double phi_4 = pi * 0.8;
	double phi_5 = pi;
	double phi_6 = pi * 1.3;
	double phi_7 = pi * 1.5;
	double phi_8 = pi * 1.8;
	double phi_9 = pi * 2.0;
	
	
	std::cout.precision( 3 );
	
	
	//----- atan2( y, x ) -----
	std::cout << atan2( sin(phi_1), cos(phi_1) ) <<"\t"<< atan2( sin(phi_2), cos(phi_2) ) <<"\t"<< atan2( sin(phi_3), cos(phi_3) ) <<"\t"<<
	             atan2( sin(phi_4), cos(phi_4) ) <<"\t"<< atan2( sin(phi_5), cos(phi_5) ) <<"\t"<< atan2( sin(phi_6), cos(phi_6) ) <<"\t"<<
				 atan2( sin(phi_7), cos(phi_7) ) <<"\t"<< atan2( sin(phi_8), cos(phi_8) ) <<"\t"<< atan2( sin(phi_9), cos(phi_9) )       << std::endl;
						 
	std::cout << phi_1 <<"\t"<<  phi_2 <<"\t"<<  phi_3 <<"\t"<<  phi_4 <<"\t"<<  phi_5 <<"\t"<<  phi_6 - 2*pi <<"\t"<<  phi_7 - 2*pi <<"\t"<<  phi_8 - 2*pi <<"\t"<< phi_9 - 2*pi << "\n" << std::endl;
	
	
	//----- Continuous range -----
	double phi_min = pi * -2.0;
	double phi_max = pi *  2.0;
	int    n       = 41;
	double delta   = ( phi_max - phi_min ) / ( (double) n - 1 );
	
	double phi = phi_min;
	for( int i=0; i<n; i++ ){
		//std::cout << phi << "\t" << atan2( sin(phi), cos(phi) ) << "\t" << phi - floor( ( abs(phi) ) / ( 2.0 * pi ) ) << std::endl; // Undeveloped continuous function
		
		double atan2_guess;
		if( ( phi >= -pi ) && ( phi <= pi ) ){
			atan2_guess = phi;
		}
		else if( ( phi >= -2.0*pi ) && ( phi <= -pi ) ){
			atan2_guess = phi + 2.0*pi;
		}
		else if( ( phi >= pi ) && ( phi <= 3.0*pi ) ){
			atan2_guess = phi - 2.0*pi;
		}
		std::cout << phi << "\t" << atan2( sin(phi), cos(phi) ) << "\t" << atan2_guess << std::endl;
		
		//std::cout << phi << "\t" << atan2( sin(phi), cos(phi) ) << "\t" << fmod( phi + pi, 2*pi ) - pi << std::endl; // Undeveloped continuous function
		
		phi += delta;
	}

	
	return 0;
}