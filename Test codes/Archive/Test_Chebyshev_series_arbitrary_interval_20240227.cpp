/*
cd OneDrive\PhD\Codes\20230810 Time evolution with updated equations\Test codes
g++ Test_Chebyshev_series_arbitrary_interval.cpp -o Test_Chebyshev_series_arbitrary_interval --std=c++11
Test_Chebyshev_series_arbitrary_interval

Test the Chebyshev series of a function f(r) on an arbitrary interval [ r_min, r_max ].
*/

#include <iostream>
#include <iomanip>
#include <vector>
#include <math.h>
#include <fstream>
#include <chrono>

#include "../../Vector_operations.h"
#include "../Header_Initial_Conditions_04.h"
#include "../Header_Vector_B_Rotating_Dipole_Ramp.h"
#include "../Header_Time_Evolution_31.h"
#include "Trapezium_rule.h"


//----- Declare functions -----
double f_function( double r );



//----- Main code -----
int main(){
	
	
	//----- Extra variables useful only for this testing code -----
	std::string output_filename_test_base = "../CSV/20240227_Test_Chebyshev_series_arbitrary_interval";	// Don't add ".csv" as the code will append _r, _t, etc to the filename for each csv.
	
	std::string output_filename_test_r = output_filename_test_base + "_r.csv";
	std::string output_filename_test_t = output_filename_test_base + "_t.csv";
	
	std::ofstream output_file_test_r;
	std::ofstream output_file_test_t;
	
	std::vector<double> r_lambda_inverse ( n_points_r );
	
	std::vector<double> chebyshev_series_integrand            ( n_points_r );
	std::vector<double> chebyshev_series_integration_variable ( n_points_r );
	
	double delta_csiv = pi / ( (double) n_points_r - 1.0 );
	
	int n_max = 10;	// Maximum Chebyshev series index.
	
	std::vector<double> chebyshev_series( n_max + 1 );
	
	
	
	//----- Calculate the gridpoints using the function as it appears in the evolution code -----
	
	calculate_gridpoints();
	
	
	
	//----- Calculate transformed gridpoints -----
	for( int i=0; i<n_points_r; i++ ){
		r_lambda_inverse[i] = ( 2.0 * r[i] - ( r_max + r_min ) ) / ( r_max - r_min );
		//std::cout << i <<"\t"<< r[i] <<"\t"<< r_lambda_inverse[i] << std::endl;
	}
	
	
	
	
	//----- Calculate Chevyshev series coefficients -----
	
	for( int i=0; i<n_points_r; i++ ){
		chebyshev_series_integration_variable[i] = (double) i * delta_csiv;
	}
	
	int n = 3;	// Put this in loop over n later.
	
	//--- Create list of integrand values ---
	for( int i=0; i<n_points_r; i++ ){
		chebyshev_series_integrand[i] = f_function( r[i] ) * cos( n * chebyshev_series_integration_variable[i] );
	}
	
	chebyshev_series[n] = ( integral_trapezium_1D_Cartesian( chebyshev_series_integrand, chebyshev_series_integrand ) );
	std::cout << chebyshev_series[n] << std::endl;
		
	
	
	
	
	//----- Code finished -----
	std::cout << "\nCode finished" << std::endl;
	
	return 0;
	
}




//----- Define functions -----
double f_function( double r ){
	return r * r;
}