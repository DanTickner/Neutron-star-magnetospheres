// cd OneDrive\PhD\Codes\Test mathematical expressions
// g++ 20220908_negative_of_spherical_vector_fulltest.cpp -o 20220908_negative_of_spherical_vector_fulltest

// -x^{2n+2} - (2n+2)! sum_{k=0}^n (-1)^k/(2k)! x^{2k} ?=? (2n+2)! sum_{k=0}^{n+1} (-1)^k/(2k)! x^{2k}

#include <iostream>
#include <math.h>

int main(){
	
	//----- Define variables -----
	const double pi = 3.14159265358979323846;
	std::cout.precision( 3 );
	
	double A_r_min     = 1e-3;	// Recall that there will be a blowing-up problem at the origin.
	double A_r_max     = 1;
	double A_theta_min = 0;
	double A_theta_max = pi;
	double A_phi_min   = 0;
	double A_phi_max   = 2*pi;
	
	int n_r     = 1;
	int n_theta = 5;
	int n_phi   = 16;
	
	double delta_A_r     = ( A_r_max     - A_r_min     ) / ( (double) n_r     - 1 );
	double delta_A_theta = ( A_theta_max - A_theta_min ) / ( (double) n_theta - 1 );
	double delta_A_phi   = ( A_phi_max   - A_phi_min   ) / ( (double) n_phi   - 1 );
	
	
	//----- Cycle through all A_r,A_theta,A_phi and calculate -A -----
	double A_r     = A_r_min;
	double A_theta = A_theta_min;
	double A_phi   = A_phi_min;
	
	
	for( int i=0; i<n_r; i++ ){
		A_theta = A_theta_min;
		for( int j=0; j<n_theta; j++ ){
			A_phi = A_phi_min;
			for( int k=0; k<n_phi; k++ ){
				
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
				
				/*double mA_phi_guess;
				if ( A_theta != 0 ){
					mA_phi_guess = A_phi - phi;
				}
				else {
					if ( A_phi < 0.5*pi ){
						mA_phi_guess = -pi;
					}
					else if( A_phi < 1.5*pi ){
						mA_phi_guess = 0;
					}
					else {
						mA_phi_guess = pi;
					}
				}*/
				
				//double mA_phi_guess = atan2( -sin(A_theta)*sin(A_phi), -sin(A_theta)*cos(A_phi) );
				double mA_phi_guess = atan2( -sin(A_theta)*sin(A_phi), -sin(A_theta)*cos(A_phi) );
				
				
				//----- Output difference between guessed and correct result -----
				std::cout << "(" << A_r            <<"\t"<<  A_theta                <<"\t"<<  A_phi              <<"\t):\t("<<
				                   mA_r            <<", "<< mA_theta                <<", "<< mA_phi              <<")\t("   <<
								   mA_r_guess      <<", "<< mA_theta_guess          <<", "<< mA_phi_guess        <<")\t"    <<
								   mA_r_guess-mA_r <<"\t"<< mA_theta_guess-mA_theta <<"\t"<< mA_phi_guess-mA_phi            << std::endl;
				
				A_phi += delta_A_phi;
			}
			std::cout << std::endl;
			A_theta += delta_A_theta;
		}
		std::cout << std::endl; std::cout << std::endl; std::cout << std::endl;
		A_r += delta_A_r;
	}
	
	std::cout << "Code finished" << std::endl;
	
	return 0;
}