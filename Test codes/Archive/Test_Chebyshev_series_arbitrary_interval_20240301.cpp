/*
cd OneDrive\PhD\Codes\20230810 Time evolution with updated equations\Test codes
g++ Test_Chebyshev_series_arbitrary_interval.cpp -o Test_Chebyshev_series_arbitrary_interval --std=c++11
Test_Chebyshev_series_arbitrary_interval

Test the Chebyshev series of a function f(r) on an arbitrary interval [ r_min, r_max ].
NOTE: In the writeup, we denote the integration variable for Chebyshev series coefficients by theta.
Here, we instead use R to avoid confusion with the polar coordinate and to highlight the fact that it's related to r, being merely a transformed radial coordinate.
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
	std::string output_filename_test_base = "20240301_Test_Chebyshev_series_arbitrary_interval";	// Don't add ".csv" as the code will append _r, _t, etc to the filename for each csv.
	
	std::string output_filename_test_csv_coordinates = "../CSV/" + output_filename_test_base + "_coordinates.csv";
	std::string output_filename_test_log             = "../Logs/" + output_filename_test_base + "_log.txt";
	
	std::ofstream output_file_test_csv_coordinates;
	std::ofstream output_file_test_log;
	
	n_points_r = 100;	// Change from value set in header files for easier testing.
	delta_r = ( r_max - r_min ) / ( (double) n_points_r - 1.0 );	// Must recalculate delta_r if n_points_r changed.
	
	std::vector<double> chebyshev_series_integrand ( n_points_r );
	std::vector<double> lambda_inverse             ( n_points_r ); // Undo the coordinate transform to check that (r_min,r_max) maps back to (-1,1).
	std::vector<double> R                          ( n_points_r ); // Integration variable for Chebyshev series coefficients.
	std::vector<double> f                          ( n_points_r ); // Function evaluated at the radial gridpoints.
	
	double delta_R = pi / ( (double) n_points_r - 1.0 );
	
	int n_max = 10;	// Maximum Chebyshev series index.
	
	std::vector<double> chebyshev_series( n_max + 1 );
	
	
	
	//----- Calculate the gridpoints using the function as it appears in the evolution code -----
	
	calculate_gridpoints();
	
	
	
	
	//----- Expected results -----
	std::vector<double> chebyshev_series_exact( n_max+1 );
	
	chebyshev_series_exact[0] = 0.5;// 0.25  * ( pow( r_min, 2 ) + pow( r_max, 2 ) );
	chebyshev_series_exact[1] = 0;// 0.125 * ( pow( r_max, 2 ) - pow( r_min, 2 ) );
	chebyshev_series_exact[2] = 0.5;// 0.125 * ( pow( r_min, 2 ) + pow( r_max, 2 ) - 2.0 * r_min * r_max );
	
	
	
	//----- Evaluate Chebyshev series integration variable and function and at the discrete radial gridpoints -----
	for( int i=0; i<n_points_r; i++ ){
		lambda_inverse[i] = ( 2.0 * r[i] - ( r_max + r_min ) ) / ( r_max - r_min );
		R             [i] = acos( lambda_inverse[i] );
		f             [i] = f_function( r[i] );
	}
	
	
	
	
	//----- Calculate Chevyshev series coefficients -----
	output_file_test_log.open( output_filename_test_log );
	
	std::cout            << "\nn\ta_n" << std::endl;
	output_file_test_log << "\nn\ta_n" << "\n";
	
	for( int n=0; n<=n_max; n++ ){
	
		//--- Create list of integrand values ---
		for( int i=0; i<n_points_r; i++ ){
			chebyshev_series_integrand[i] = f_function( r[i] ) * cos( n * R[i] );
			//std::cout << n <<"\t"<< i <<"\t"<< f_function( r[i] ) <<"\t"<< cos( n * R[i] ) <<"\t"<< chebyshev_series_integrand[i] << std::endl;
		}
		
		//--- Calculate and output Chebyshev series coefficient ---
		//chebyshev_series[n] = ( integral_trapezium_1D_Cartesian_with_cout( chebyshev_series_integrand, chebyshev_series_integrand ) );
		chebyshev_series[n] = ( integral_trapezium_1D_Cartesian( chebyshev_series_integrand, chebyshev_series_integrand ) );
		
		std::cout            << n <<"\t"<< chebyshev_series[n] <<"\t"<< chebyshev_series_exact[n] << std::endl;
		output_file_test_log << n <<"\t"<< chebyshev_series[n] <<"\t"<< chebyshev_series_exact[n] << "\n";
		
	}
	
	std::cout << "Log file saved:\t" << output_filename_test_log << std::endl;
	
	
	
	
	//----- Build Chebyshev polynomials -----
	std::vector< std::vector<double> > T_n ( n_max+1, std::vector<double> ( n_points_r ) );
	
	for( int n=0; n<=n_max; n++ ){
		for( int i=0; i<n_points_r; i++ ){
			T_n[n][i] = cos( n * R[i] );
		}
	}
	
	
	//----- Output coordinates and Chebyshev polynomials to CSV -----
	output_file_test_csv_coordinates.open( output_filename_test_csv_coordinates );
	output_file_test_csv_coordinates << "i,r,Lambda_inverse,R";
	for( int n=0; n<=n_max; n++ ){
		output_file_test_csv_coordinates << ",T_" << n;
	}
	output_file_test_csv_coordinates << "\n";
	
	for( int i=0; i<n_points_r; i++ ){
		
		output_file_test_csv_coordinates << i <<","<< r[i] <<","<< lambda_inverse[i] <<","<< R[i];
		
		for( int n=0; n<=n_max; n++ ){
			output_file_test_csv_coordinates <<"," << T_n[n][i];
		}
		output_file_test_csv_coordinates << "\n";
	
	}
	
	std::cout << "\nCSV file saved:\t"<< output_filename_test_csv_coordinates << std::endl;
	
	
	
	//----- Output expected Chebyshev series evaluated -----
	
	for( int i=0; i<n_points_r; i++ ){
		
		double result = 0;
		
		for( int n=0; n<=n_max; n++ ){
			result += chebyshev_series_exact[n] * T_n[n][i];
		}
		
		//std::cout << i <<"\t"<< r[i] <<"\t"<< f_function( r[i] ) <<"\t"<< result << std::endl;
		
	}
			
		
	
	
	
	
	
	//----- Code finished -----
	std::cout << "\nCode finished" << std::endl;
	
	return 0;
	
}




//----- Define functions -----
double f_function( double r ){
	return r * r;
}