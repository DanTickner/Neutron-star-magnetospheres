/*
cd OneDrive\PhD\Codes\20240311 Time evolution code\Test codes
g++ Test_FF_region_expression_dipole.cpp -o Test_FF_region_expression_dipole --std=c++11
Test_FF_region_expression_dipole

In the report, we showed that a magnetic dipole rotating with angular velocity Omega satisfies the force-free condition
B^2 - E^2 > 0
for all points satisfying
r > R_LC csc(theta)
where R_LC = 1 / Omega.
Thereforce, r * sin(theta) < R_LC.
Test this expression by calculating B^2 and E^2 at all gridpoints and comparing it to r * sin(theta).
*/

#include <iostream>
#include <iomanip>
#include <vector>
#include <math.h>
#include <fstream>
#include <chrono>

#include "../Header_Initial_Conditions_08.h"
#include "../Header_Time_Evolution_40.h"


//----- Declare functions -----
std::vector<double> cross_product( std::vector<double> a, std::vector<double> b );



//----- Main code -----
int main(){
	
	
	//----- Extra variables useful only for this testing code -----
	std::string output_filename_test_base = "20240412_Test_FF_region_expression_dipole";	// Don't add ".csv" as the code will append _r, _t, etc to the filename for each csv.
	
	std::string output_filename_test_csv = "../CSV/"  + output_filename_test_base + "_function_values.csv";
	
	std::ofstream output_file_test_csv;
	
	double Bsq_minus_Esq_tol = 1e-5;	// Tolerance above and below zero at which we accept that B^2-E^2 = 0.
	bool   output_to_csv = true;		// Option to output all gridpoints to a CSV file.
	int    n_gridpoints_failsafe = 1e5;	// Don't output to CSV if the total number of gridpoints exceeds this.
	
	std::vector< std::vector<double> > r_sin_theta_minus_R_LC     ( n_r, std::vector<double> ( n_t ) );
	
	
	
	//----- Calculate the gridpoints using the function as it appears in the evolution code -----
	calculate_time_and_rotation_values();
	calculate_gridpoints();
	apply_initial_field_values_for_B();
	
	
	//----- Calculate E for a constantly rotating dipole -----
	for( int i=0; i<n_r; i++ ){
		for( int j=0; j<n_t; j++ ){
			
			std::vector<double> vector_r {
				r[i],
				0,
				0
			};
			
			std::vector<double> vector_Omega {
				Omega.back() * cos(t[j]),
				Omega.back() * sin(t[j]),
				0
			};
			
			std::vector<double> vector_B {
				B[0][i][j],
				B[1][i][j],
				B[2][i][j]
			};
			
			std::vector<double> Omega_cross_r             = cross_product( vector_Omega , vector_r );
			std::vector<double> Omega_cross_r_all_cross_B = cross_product( Omega_cross_r, vector_B );
			
			for( int f=0; f<3; f++ ){
				E[f][i][j] = - Omega_cross_r_all_cross_B[f];
			}
			
		}
	}
	
	
	//----- Calculate B^2 - E^2 at all gridpoints -----
	calculate_dot_products();
	
	
	//----- Calculate r sin(theta) - R_LC -----
	// All points where this is greater than zero satisfy the FF condition.
	
	for( int i=0; i<n_r; i++ ){
		for( int j=0; j<n_t; j++ ){
			r_sin_theta_minus_R_LC[i][j] = r[i] * st[j] - R_LC.back();
		}
	}
	
	
	//----- Output values to CSV -----
	int n_gridpoints = n_r * n_t;
	
	std::cout << "n_gridpoints         :\t" << n_gridpoints          << std::endl;
	std::cout << "n_gridpoints_failsafe:\t" << n_gridpoints_failsafe << std::endl;
	
	if( not( output_to_csv ) ){
		std::cout << "CSV file not created, as requested." << std::endl;
	}
	
	else if( n_gridpoints > n_gridpoints_failsafe ){
		std::cout << "n_gridpoints > n_gridpoints_failsafe. CSV file not created." << std::endl;
	}
	
	else {
		
		output_file_test_csv.open( output_filename_test_csv );
		output_file_test_csv << "i,j,r,t,E_dot_B,B_squared_minus_E_squared,r_sin_theta_minus_R_LC,in_FF_region_based_on_FF_condition,in_FF_region_based_on_expression,are_conditions_equal,x,r_sin_theta\n";
		
		for( int i=0; i<n_r; i++ ){
			for( int j=0; j<n_t; j++ ){
				
				bool in_FF_region_based_on_FF_condition = ( B_squared_minus_E_squared[i][j] > 0 );
				bool in_FF_region_based_on_expression   = ( r[i] * st[j] < R_LC.back() );
				bool are_conditions_equal               = ( in_FF_region_based_on_expression == in_FF_region_based_on_expression );
				
				output_file_test_csv << i <<","<< j <<","<< r[i] <<","<< t[j] <<","<< E_dot_B[i][j] <<","<< B_squared_minus_E_squared[i][j] <<","<< r_sin_theta_minus_R_LC[i][j]
				                     <<","<< in_FF_region_based_on_FF_condition <<","<< in_FF_region_based_on_expression <<","<< are_conditions_equal <<","<< x[i][j] <<","<< r[i]*st[j]<< "\n";
			}			
		}
		
		std::cout << "\nCSV file saved:\t" << output_filename_test_csv << std::endl;
		
	}
	
	std::cout << "\nCode finished" << std::endl;
	
	return 0;
	
}


//----- Define functions -----
std::vector<double> cross_product( std::vector<double> a, std::vector<double> b ){
	return std::vector<double> {
		a[1]*b[2] - a[2]*b[1],
		a[2]*b[0] - a[0]*b[2],
		a[0]*b[1] - a[1]*b[0]
	};
}