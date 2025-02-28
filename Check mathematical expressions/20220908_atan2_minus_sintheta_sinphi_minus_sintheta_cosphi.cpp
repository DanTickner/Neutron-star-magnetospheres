// cd OneDrive\PhD\Codes\Test mathematical expressions
// g++ 20220908_atan2_minus_sintheta_sinphi_minus_sintheta_cosphi.cpp -o 20220908_atan2_minus_sintheta_sinphi_minus_sintheta_cosphi

// atan2( sin(phi), cos(phi) ) = ?

#include <iostream>
#include <math.h>

int main(){
	
	//----- Define variables -----
	const double pi         = 3.14159265358979323846;
	double      theta_min   = 0.1;														// My derived result should mean that atan2(...) is undefined if theta=0.
	double      theta_max   = pi;
	double      phi_min     = 0;
	double      phi_max     = pi*2.0;
	int         n_theta     = 16;
	int         n_phi       = 16;
	double      delta_theta = ( theta_max - theta_min ) / ( (double) n_theta - 1 );
	double      delta_phi   = ( phi_max   - phi_min   ) / ( (double) n_phi   - 1 );
	
	std::cout.precision( 3 );
	
	//----- Cycle through all theta,phi and evaluate the function -----
	double theta = theta_min;
	double phi   = phi_min;
	for( int i=0; i<n_theta; i++ ){
		phi = phi_min;
		for( int j=0; j<n_phi; j++ ){
			
			double atan2_exact = atan2( -sin(theta)*sin(phi) , -sin(theta)*cos(phi) );
			
			/*double atan2_guess;
			if( theta == 0 ){
				atan2_guess = 0;
			} else {
				atan2_guess = phi - pi;
			}*/
			double atan2_guess = phi - pi;
			
			std::cout << theta <<", "<< phi <<":\t" << atan2_exact <<"\t"<< atan2_guess <<"\t\t"<< atan2_exact - atan2_guess << std::endl;
			phi += delta_phi;
		}
		theta += delta_theta;
	}

	
	return 0;
}