/*
cd OneDrive\PhD\Codes\20240311 Time evolution code\Test codes
g++ Test_curl_from_VSH_series_axisymmetric.cpp -o Test_curl_from_VSH_series_axisymmetric --std=c++11
Test_curl_from_VSH_series_axisymmetric

Test the calculation of the divergence of a vector given its VSH series.
For ease of using the main code, the vector is referred to as B.
Set the known exact value of div(B) in the function double B_div_exact( double r, double t ).
The accuracy of the VSH series itself is already tested in Test_VSH_series_axisymmetric.cpp.

Output a CSV with the calculated and exact divergence of B for each gridpoint.

Calculate the stdev over all gridpoints.
*/

#include <iostream>
#include <iomanip>
#include <vector>
#include <math.h>
#include <fstream>
#include <chrono>

#include "../../Vector_operations.h"
#include "../Header_Initial_Conditions_08.h"
#include "../Header_Time_Evolution_38.h"




//----- Declare functions -----
std::vector<double> B_curl_exact_function( double r, double t );




//----- Main code -----
int main(){
	
	
	//----- Extra variables useful only for this testing code -----
	
	int test_type = 3;	// 1 = Original; 2 = Divergenceless with finite differencing; 3 = Divergenceless with Chebyshev decomposition; 
	
	std::string output_filename_test_base = "20240329_Test_VSH_series_axisymmetric_divergence";
	
	switch( test_type ){
		case 1:
			output_filename_test_base += "_Original";
			std::cout << "Test type:\tOriginal" << std::endl;
			break;
		case 2:
			output_filename_test_base += "_Divergenceless_FD";
			std::cout << "Test type:\tDivergenceless with finite differencing" << std::endl;
			break;
		case 3:
			output_filename_test_base += "_Divergenceless_Chebyshev_nmax_" + std::to_string( n_max );
			std::cout << "Test type:\tDivergenceless with Chebyshev series" << std::endl;
			break;
		default:
			std::cout << "Test type not recognised. Stopping the code.";
			return 1;
	}
	
	std::string output_filename_test_values = "../CSV/"  + output_filename_test_base + "_values.csv";
	std::string output_filename_test_log    = "../Logs/" + output_filename_test_base + ".txt";
	
	std::ofstream output_file_test_values;
	std::ofstream output_file_test_log;
	
	std::vector< std::vector< double> > B_div_exact    ( n_r, std::vector<double> ( n_t ) );	// Exact values of div(B).
	
	
	//----- Calculate the gridpoints, associated Legendre functions and VSH series coefficients, recalculate vector and calculate stdev of vector using the functions in the evolution code -----
	calculate_time_and_rotation_values();	// Needed if bool set_r_max_to_2_R_LC = true.
	calculate_gridpoints();
	calculate_associated_legendre_functions();
	calculate_sqrt_2Lplus1_over_4pi();
	calculate_remapped_chebyshev_polynomials();
	apply_initial_field_values_for_B();
	
	
	//----- Perform the VSH decomposition according to the method chosen -----
	
	B_VSH_coeffs[0] = VSH_decomposition_r( B );
	B_VSH_coeffs[2] = VSH_decomposition_2( B );
	
	switch( test_type ){
		case 1:
			B_VSH_coeffs[1] = VSH_decomposition_1( B );
			for( int f=0; f<3; f++ ){
				B_VSH_coeffs_dr[f] = radial_derivatives_of_VSH_coeff_FD( B_VSH_coeffs, f );
			}
			break;
		case 2:
			B_VSH_coeffs_dr[0] = radial_derivatives_of_VSH_coeff_FD( B_VSH_coeffs, 0 );
			B_VSH_coeffs_dr[2] = radial_derivatives_of_VSH_coeff_FD( B_VSH_coeffs, 2 );
			calculate_B_1_L_and_B_1_L_dr();
			break;
		case 3:
			B_VSH_coeffs_dr[0] = radial_derivatives_of_VSH_coeff_Chebyshev( B_VSH_coeffs, 0 );
			B_VSH_coeffs_dr[2] = radial_derivatives_of_VSH_coeff_Chebyshev( B_VSH_coeffs, 2 );
			calculate_B_1_L_and_B_1_L_dr();
			break;
		default:
			std::cout << "Test type not recognised. Stopping the code.";
			return 1;
	}
	
	
	//----- Calculate divergence and obtain exact values -----
	B_div = divergence_from_VSH_coeffs( B_VSH_coeffs, B_VSH_coeffs_dr );
	
	for( int i=0; i<n_r; i++ ){
		for( int j=0; j<n_t; j++ ){
			B_div_exact[i][j] = B_div_exact_function( r[i], t[j] );
		}
	}
	
	
	//----- Output to CSV -----
	output_file_test_values.open( output_filename_test_values );
	
	output_file_test_values << "i,j,r,t,B_div,B_div_exact\n";
	
	for( int i=0; i<n_r; i++ ){
		for( int j=0; j<n_t; j++ ){
			output_file_test_values << i <<","<< j <<","<< r[i] <<","<< t[j] <<","<< B_div[i][j] <<","<< B_div_exact[i][j] <<"\n";
		}
	}
	
	
	//----- Calculate stdev across all gridpoints -----
	double stdev = stdev_of_array_2d( B_div_exact, B_div );
	std::vector<double> max_dev = max_deviation_of_array_2d( B_div_exact, B_div );
	
	std::cout            << "Stdev of divergence across all gridpoints:\t" << stdev << std::endl;
	output_file_test_log << "Stdev of divergence across all gridpoints:\t" << stdev << "\n";
	
	std::cout            << "Max abs deviation of divergence across all gridpoints:\t" << max_dev[0] <<"\t at ( i=" << max_dev[1] <<" , j=" << max_dev[2] << " )" << std::endl;
	output_file_test_log << "Max abs deviation of divergence across all gridpoints:\t" << max_dev[0] <<"\t at ( i=" << max_dev[1] <<" , j=" << max_dev[2] << " )\n";
	
	
	//----- Code finished -----
	std::cout << "\nCSV file saved:\t" << output_filename_test_values << std::endl;
	std::cout <<   "Log file saved:\t" << output_filename_test_log    << std::endl;
	std::cout << "\nCode finished" << std::endl;
	
	return 0;
	
}




//----- Define functions -----
double B_div_exact_function( double r, double t ){
	return 0;
}