/*
cd OneDrive\PhD\Codes\20240311 Time evolution code\Test codes
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
#include "../Header_Initial_Conditions_08.h"
#include "../Header_Time_Evolution_38.h"


//----- Declare functions -----
double f_function( double r );
double dfdx_function( double r );
double a_n_exact_function( int n );



//----- Main code -----
int main(){
	
	
	//----- Extra variables useful only for this testing code -----
	std::string output_filename_test_base = "20240404_Test_Chebyshev_series_arbitrary_interval_nmax_" + std::to_string( n_max ) + "r_squared";	// Don't add ".csv" as the code will append _r, _t, etc to the filename for each csv.
	
	std::string output_filename_test_csv = "../CSV/"  + output_filename_test_base + "_function_values.csv";
	std::string output_filename_test_log = "../Logs/" + output_filename_test_base + "_log.txt";
	
	std::ofstream output_file_test_csv;
	std::ofstream output_file_test_log;
	
	std::vector<double> f_exact     ( n_r ); // Exact values of the function.
	std::vector<double> f_series    ( n_r ); // Chebyshev series evaluated at the radial gridpoints.
	std::vector<double> dfdx_exact  ( n_r ); // Exact values of the first derivative.
	std::vector<double> dfdx_series ( n_r ); // First derivative of the function taken by differentiating the Chebyshev series term-by-term.
	
	
	
	//----- Calculate the gridpoints using the function as it appears in the evolution code -----
	calculate_time_and_rotation_values();
	calculate_gridpoints();
	calculate_remapped_chebyshev_polynomials();
	
	
	//----- Evaluate function at the radial gridpoits -----
	for( int i=0; i<n_r; i++ ){
		f_exact   [i] = f_function   ( r[i] );
		dfdx_exact[i] = dfdx_function( r[i] );
	}
	
	
	//----- Calculate Chebyshev series coefficients and output to screen -----
	w = precision_csv + 5;
	std::vector<double> a_n = chebyshev_series_coeffs( f_exact );
	
	output_file_test_log.open( output_filename_test_log );
	
	std::cout            <<std::left<<std::setw(w)<< "n" <<std::left<<std::setw(w)<< "numerical" <<std::left<<std::setw(w)<< "exact" <<std::left<<std::setw(w)<< "relative_error" << std::endl;
	output_file_test_log <<std::left<<std::setw(w)<< "n" <<std::left<<std::setw(w)<< "numerical" <<std::left<<std::setw(w)<< "exact" <<std::left<<std::setw(w)<< "relative_error" << "\n";
	
	for( int n=0; n<=n_max; n++ ){
		
		double numerical = a_n[n];
		double exact     = a_n_exact_function( n );
		
		double relative_error = 0;
		if( numerical != 0 ){
			relative_error = 1.0 - exact / numerical;
		}
		
		std::cout            <<std::setprecision(precision_csv)<<std::left<<std::setw(w)<< n <<std::left<<std::setw(w)<< numerical <<std::left<<std::setw(w)<< exact <<std::left<<std::setw(w)<< relative_error << std::endl;
		//output_file_test_log <<std::setprecision(precision_csv)<<std::left<<std::setw(w)<< n <<std::left<<std::setw(w)<< numerical <<std::left<<std::setw(w)<< exact <<std::left<<std::setw(w)<< relative_error << "\n";
		output_file_test_log <<std::setprecision(precision_csv)<< n <<" & "<< numerical <<" & "<< exact <<" & "<< relative_error << "\\\\\n";
	}
	
	
	//----- Evaluate Chebyshev series and derivative at all points, calculate stdev and output to CSV -----
	output_file_test_csv.open( output_filename_test_csv );
	output_file_test_csv << "i,r,f_exact,f_series,dfdx_exact,dfdx_series\n";
	
	for( int i=0; i<n_r; i++ ){
		
		f_series   [i] = evaluate_chebyshev_series_at_gridpoint( a_n, i );
		dfdx_series[i] = evaluate_chebyshev_series_first_derivative_at_gridpoint( a_n, i );
		
		output_file_test_csv << i <<","<< r[i] <<","<< f_exact[i] <<","<< f_series[i] <<","<< dfdx_exact[i] <<","<< dfdx_series[i] <<"\n";
	}
	
	double f_stdev    = stdev_of_array_1d( f_exact   , f_series    );
	double dfdx_stdev = stdev_of_array_1d( dfdx_exact, dfdx_series );
	
	
	std::cout            << "\nStdev f(x) :\t" << f_stdev    << std::endl;
	output_file_test_log << "\nStdev f(x) :\t" << f_stdev    << "\n";
	std::cout            <<   "Stdev df/dx:\t" << dfdx_stdev << std::endl;
	output_file_test_log <<   "Stdev df/dx:\t" << dfdx_stdev << "\n";

	
	
	std::cout << "\nCSV file saved:\t" << output_filename_test_csv << std::endl;
	std::cout <<   "Log file saved:\t" << output_filename_test_log << std::endl;
	
	
	
	//----- Code finished -----
	std::cout << "\nCode finished" << std::endl;
	
	return 0;
	
}




//----- Define functions -----
double f_function( double r ){
	return r * r;
	//return pow( r, -3 );
}

double dfdx_function( double r ){
	return 2.0 * r;
	//return -3.0 * pow( r, -4 );
}

double a_n_exact_function( int n ){
	
	switch( n ){
		case 0:
			return ( 3.0*pow(r_max,2) + 3.0*pow(r_min,2) + 2.0*r_max*r_min ) / 8.0;
		case 1:
			return ( r_max - r_min ) * ( r_max + r_min ) / 2.0;
		case 2:
			return pow( r_max - r_min, 2 ) / 8.0;
		default:
			return 0;
	}
	
	
	return 0;
	
}