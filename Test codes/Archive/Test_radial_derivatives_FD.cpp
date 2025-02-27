/*
cd OneDrive\PhD\Codes\20240311 Time evolution code\Test codes
g++ Test_radial_derivatives_FD.cpp -o Test_radial_derivatives_FD --std=c++11
Test_radial_derivatives_FD

Test the first and second radial derivatives of a function f(r).

REASON ARCHIVED: More accurate finite-difference expressions developed. Test them to different accuracy before adding one to the main code.
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
double f_function    ( double r );
double f_dr_function ( double r );
double f_dr2_function( double r );



//----- Main code -----
int main(){
	
	
	//----- Extra variables useful only for this testing code -----
	std::string output_filename_test_base = "20240501_Test_radial_derivatives_FD_enhanced_endpoints_1_over_r3";	// Don't add ".csv" as the code will append _r, _t, etc to the filename for each csv.
	
	std::string output_filename_test_csv = "../CSV/"  + output_filename_test_base + "_function_values.csv";
	std::string output_filename_test_log = "../Logs/" + output_filename_test_base + "_log.txt";
	
	std::ofstream output_file_test_csv;
	std::ofstream output_file_test_log;
	
	std::vector<double> f          ( n_r );	// Values of the function.
	std::vector<double> f_dr_exact ( n_r );	// First  derivative: exact values.
	std::vector<double> f_dr_FD    ( n_r );	// First  derivative: values from finite differencing.
	std::vector<double> f_dr2_exact( n_r );	// Second derivative: exact values.
	std::vector<double> f_dr2_FD   ( n_r );	// Second derivative: values from finite differencing.
	
	double f_dr_stdev  = 0;
	double f_dr2_stdev = 0;
	
	
	
	//----- Calculate the gridpoints using the function as it appears in the evolution code -----
	calculate_time_and_rotation_values();
	calculate_gridpoints();
	
	
	//----- Build list of values of the function and its exact derivatives -----
	for( int i=0; i<n_r; i++ ){
		f          [i] = f_function    ( r[i] );
		f_dr_exact [i] = f_dr_function ( r[i] );
		f_dr2_exact[i] = f_dr2_function( r[i] );
	}
	
	
	//----- Calculate derivatives by finite differencing -----
	f_dr_FD  = radial_derivatives_1_FD( f );
	f_dr2_FD = radial_derivatives_2_FD( f );
	
	
	//----- Output values to CSV -----
	output_file_test_csv.open( output_filename_test_csv );
	output_file_test_csv << "i,r,f,f_dr_exact,f_dr_FD,abs_rel_err_f_dr,f_dr2_exact,f_dr2_FD,abs_rel_err_f_dr2\n";
	
	for( int i=0; i<n_r; i++ ){
		
		double abs_rel_err_f_dr  = 0;
		double abs_rel_err_f_dr2 = 0;
		if( f_dr_exact[i] != 0 ){
			abs_rel_err_f_dr = abs( 1.0 - f_dr_FD[i] / f_dr_exact[i] );
		}
		if( f_dr2_exact[i] != 0 ){
			abs_rel_err_f_dr2 = abs( 1.0 - f_dr2_FD[i] / f_dr2_exact[i] );
		}
		
		output_file_test_csv << i <<","<< r[i] <<","<< f[i] <<","<< f_dr_exact [i] <<","<< f_dr_FD [i] <<","<< abs_rel_err_f_dr
		                                                    <<","<< f_dr2_exact[i] <<","<< f_dr2_FD[i] <<","<< abs_rel_err_f_dr2 <<"\n";
		
	}
	
	
	//----- Standard deviation of derivatives -----
	
	f_dr_stdev  = stdev_of_array_1d( f_dr_exact , f_dr_FD  );
	f_dr2_stdev = stdev_of_array_1d( f_dr2_exact, f_dr2_FD );
	
	output_file_test_log.open( output_filename_test_log );
	
	std::cout            << "Standard deviation of first  derivative:\t" << f_dr_stdev  << std::endl;
	output_file_test_log << "Standard deviation of first  derivative:\t" << f_dr_stdev  << "\n";
	std::cout            << "Standard deviation of second derivative:\t" << f_dr2_stdev << std::endl;
	output_file_test_log << "Standard deviation of second derivative:\t" << f_dr2_stdev << "\n";
	
	
	//----- Code finished -----
	std::cout << "\nCSV file saved:\t" << output_filename_test_csv << std::endl;
	std::cout <<   "Log file saved:\t" << output_filename_test_log << std::endl;
	
	std::cout << "\nCode finished" << std::endl;
	
	return 0;
	
}




//----- Define functions -----
double f_function    ( double r ){
	//return r * r;
	return pow( r, -3 );
	//return sin( r * r );
}

double f_dr_function ( double r ){
	//return 2.0 * r;
	return -3.0 * pow( r, -4 );
	//return 2.0*r * cos( r*r );
}

double f_dr2_function( double r ){
	//return 2.0;
	return 12.0 * pow( r, -5 );
	//return 2.0 * cos( r*r ) - 4.0 * r*r * sin( r*r );
}