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
#include "../Header_Initial_Conditions_06.h"
#include "../Header_Time_Evolution_33.h"


//----- Declare functions -----
double f_function( double r );
double dfdx_function( double r );



//----- Main code -----
int main(){
	
	
	//----- Extra variables useful only for this testing code -----
	std::string output_filename_test_base = "20240308_Test_Chebyshev_series_arbitrary_interval";	// Don't add ".csv" as the code will append _r, _t, etc to the filename for each csv.
	
	std::string output_filename_test_csv = "../CSV/" + output_filename_test_base + "_function_values.csv";
	std::string output_filename_test_log             = "../Logs/" + output_filename_test_base + "_log.txt";
	
	std::ofstream output_file_test_csv;
	std::ofstream output_file_test_log;
	
	std::vector<double> f_exact     ( n_points_r ); // Exact values of the function.
	std::vector<double> f_series    ( n_points_r ); // Chebyshev series evaluated at the radial gridpoints.
	std::vector<double> dfdx_exact  ( n_points_r ); // Exact values of the first derivative.
	std::vector<double> dfdx_series ( n_points_r ); // First derivative of the function taken by differentiating the Chebyshev series term-by-term.
	
	
	
	//----- Calculate the gridpoints using the function as it appears in the evolution code -----
	
	calculate_gridpoints();
	calculate_remapped_chebyshev_polynomials();
	
	
	//----- Evaluate function at the radial gridpoits -----
	for( int i=0; i<n_points_r; i++ ){
		f_exact   [i] = f_function   ( r[i] );
		dfdx_exact[i] = dfdx_function( r[i] );
	}
	
	
	//----- Calculate Chebyshev series coefficients and output to screen -----
	std::vector<double> a_n = calculate_chebyshev_series_coeffs( f_exact );
	
	output_file_test_log.open( output_filename_test_log );
	
	std::cout            << "n\ta_n" << std::endl;
	output_file_test_log << "n\ta_n\n";
	
	for( int n=0; n<=n_max; n++ ){
		
		std::cout            << n <<"\t"<< a_n[n] << std::endl;
		output_file_test_log << n <<"\t"<< a_n[n] <<"\n";
	}
	
	
	//----- Evaluate Chebyshev series and derivative at all points, calculate stdev and output to CSV -----
	output_file_test_csv.open( output_filename_test_csv );
	output_file_test_csv << "i,r,f_exact,f_series,dfdx_exact,dfdx_series\n";
	
	for( int i=0; i<n_points_r; i++ ){
		
		f_series   [i] = evaluate_chebyshev_series_at_gridpoint( a_n, i );
		dfdx_series[i] = evaluate_chebyshev_series_first_derivative_at_gridpoint( a_n, i );
		
		output_file_test_csv << i <<","<< r[i] <<","<< f_exact[i] <<","<< f_series[i] <<","<< dfdx_exact[i] <<","<< dfdx_series[i] <<"\n";
	}
	
	double f_stdev    = stdev_of_array( f_exact   , f_series    );
	double dfdx_stdev = stdev_of_array( dfdx_exact, dfdx_series );
	
	
	std::cout            << "\nStdev f(x) :\t" << f_stdev    << std::endl;
	output_file_test_log << "\nStdev f(x) :\t" << f_stdev    << "\n";
	std::cout            <<   "Stdev df/dx:\t" << dfdx_stdev << std::endl;
	output_file_test_log <<   "Stdev df/dx:\t" << dfdx_stdev << "\n";
	
	
	
	
	
	
	
	
	
	
	std::cout << "\nCSV file saved:\t" << output_filename_test_csv << std::endl;
	std::cout <<   "Log file saved:\t" << output_filename_test_log << std::endl;
	
		
	
	
	
	
	/*
	
	
	
	
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
			
		
	
	*/
	
	
	
	//----- Code finished -----
	std::cout << "\nCode finished" << std::endl;
	
	return 0;
	
}




//----- Define functions -----
double f_function( double r ){
	return r * r;
}

double dfdx_function( double r ){
	return 2.0 * r;
}