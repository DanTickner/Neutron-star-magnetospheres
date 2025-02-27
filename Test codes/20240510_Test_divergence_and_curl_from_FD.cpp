/*
cd OneDrive\PhD\Codes\20240311 Time evolution code\Test codes
g++ 20240510_Test_divergence_and_curl_from_FD.cpp -o 20240510_Test_divergence_and_curl_from_FD --std=c++11
20240510_Test_divergence_and_curl_from_FD

Test the calculation of the divergence and curl of a vector given its VSH series.
For ease of using the main code, the vector is referred to as B.
Set the known exact values of div(B) and curl(B) in the functions double B_div_exact( double r, double t ) and std::vector<double> B_curl_exact( double r, double t ).
The accuracy of the VSH series itself is already tested in Test_VSH_series_axisymmetric.cpp.

Output a CSV with the calculated and exact divergence and curl of B for each gridpoint.

Calculate the stdev over all gridpoints.
*/

#include <iostream>
#include <iomanip>
#include <vector>
#include <math.h>
#include <fstream>
#include <chrono>

#include "../../Vector_operations.h"
#include "../Header_Initial_Conditions_09.h"
#include "../Header_Time_Evolution_43.h"




//----- Declare functions -----
double B_r_function_test( double r, double t );
double B_t_function_test( double r, double t );
double B_p_function_test( double r, double t );

double B_div_exact_function( double r, double t );

double B_curl_exact_function_r( double r, double t );
double B_curl_exact_function_t( double r, double t );
double B_curl_exact_function_p( double r, double t );




//----- Main code -----
int main(){
	
	
	//----- Extra variables useful only for this testing code -----
	
	int test_type = 1;	// 1 = General; 2 = Divergenceless with finite differencing; 3 = Divergenceless with Chebyshev decomposition;
	
	std::string output_filename_test_base = "20240510_Test_divergence_and_curl_from_FD";
	
	switch( test_type ){
		case 1:
			output_filename_test_base += "_General";
			std::cout << "Test type:\tGeneral" << std::endl;
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
	
	std::vector< std::vector<double> > B_div_exact ( n_r, std::vector<double> ( n_t ) );	// Exact values of div(B).
	std::vector< std::vector<double> > B_div_ARE   ( n_r, std::vector<double> ( n_t ) );	// Absolute relative error.
	std::vector< std::vector< std::vector<double> > > B_curl_exact ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );	// Exact values of curl(B).
	std::vector< std::vector< std::vector<double> > > B_curl_ARE   ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );	// Absolute relative error.
	
	
	//----- Calculate the gridpoints, associated Legendre functions and remapped Chebyshev polynomials using the functions in the evolution code -----
	calculate_gridpoints();
	calculate_functions_of_gridpoints();
	
	
	//----- Apply magnetic field values -----
	for( int i=0; i<n_r; i++ ){
		for( int j=0; j<n_t; j++ ){
			B[0][i][j] = B_r_function_test( r[i], t[j] );
			B[1][i][j] = B_t_function_test( r[i], t[j] );
			B[2][i][j] = B_p_function_test( r[i], t[j] );
		}
	}
	
	
	//----- Calculate radial derivatives and hence divergence and curl of vector -----
	std::vector< std::vector< std::vector<double> > > B_dr ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );
	std::vector< std::vector< std::vector<double> > > B_dt ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );
	std::vector< std::vector< std::vector<double> > > E_dr ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );
	std::vector< std::vector< std::vector<double> > > E_dt ( 3, std::vector< std::vector<double> > ( n_r, std::vector<double> ( n_t ) ) );
	
	for( int f=0; f<3; f++ ){
		B_dr[f] = f_dx_FD( B[f], delta_r );
		B_dt[f] = f_dy_FD( B[f], delta_t );
		E_dr[f] = f_dx_FD( E[f], delta_r );
		E_dt[f] = f_dy_FD( E[f], delta_t );
	}
	
	
	
	//----- Calculate divergence and curl, and obtain exact values and absolute relative error -----
	B_div   = divergence_from_FD_derivatives( B, B_dr, B_dt );
	B_curl  = curl_from_FD_derivatives      ( B, B_dr, B_dt );
	
	std::cout << "code gets here" << std::endl; 
	for( int i=0; i<n_r; i++ ){
		for( int j=0; j<n_t; j++ ){
			
			B_div_exact[i][j] = B_div_exact_function( r[i], t[j] );
			
			B_curl_exact[0][i][j] = B_curl_exact_function_r( r[i], t[j] );
			B_curl_exact[1][i][j] = B_curl_exact_function_t( r[i], t[j] );
			B_curl_exact[2][i][j] = B_curl_exact_function_p( r[i], t[j] );
			
			B_div_ARE[i][j] = ( B_div_exact[i][j] == 0 ) ? ( 0 ) : ( abs( 1.0 - B_div[i][j] / B_div_exact[i][j] ) );
			
			for( int f=0; f<3; f++ ){
				B_curl_ARE[f][i][j] = ( B_curl_exact[f][i][j] == 0 ) ? ( 0 ) : ( abs( 1.0 - B_curl[f][i][j] / B_curl_exact[f][i][j] ) );
			}
			
		}
	}
			
	
	
	
	//----- Output to CSV -----
	output_file_test_values.open( output_filename_test_values );
	
	output_file_test_values << "i,j,r,t,B_div,B_div_exact,B_div_ARE,B_curl_r,B_curl_r_exact,B_curl_r_ARE,B_curl_t,B_curl_t_exact,B_curl_t_ARE,B_curl_p,B_curl_p_exact,B_curl_p_ARE\n";
	
	for( int i=0; i<n_r; i++ ){
		for( int j=0; j<n_t; j++ ){
			
			output_file_test_values << i <<","<< j <<","<< r[i] <<","<< t[j] <<","<< B_div[i][j] <<","<< B_div_exact[i][j] <<","<< B_div_ARE[i][j];
			
			for( int f=0; f<3; f++ ){
				output_file_test_values <<","<< B_curl[f][i][j] <<","<< B_curl_exact[f][i][j] <<","<< B_curl_ARE[f][i][j];
			}
			
			output_file_test_values << "\n";
			
		}
	}
	
	
	//----- Calculate stdev across all gridpoints -----
	double stdev_div    = stdev_of_array_2d( B_div_exact    , B_div     );
	double stdev_curl_r = stdev_of_array_2d( B_curl_exact[0], B_curl[0] );
	double stdev_curl_t = stdev_of_array_2d( B_curl_exact[1], B_curl[1] );
	double stdev_curl_p = stdev_of_array_2d( B_curl_exact[2], B_curl[2] );
	
	std::vector<double> max_dev_div    = max_deviation_of_array_2d( B_div_exact    , B_div     );
	std::vector<double> max_dev_curl_r = max_deviation_of_array_2d( B_curl_exact[0], B_curl[0] );
	std::vector<double> max_dev_curl_t = max_deviation_of_array_2d( B_curl_exact[1], B_curl[1] );
	std::vector<double> max_dev_curl_p = max_deviation_of_array_2d( B_curl_exact[2], B_curl[2] );
	
	std::cout            << "Stdev of divergence across all gridpoints:\t" << stdev_div << std::endl;
	output_file_test_log << "Stdev of divergence across all gridpoints:\t" << stdev_div << "\n";
	
	std::cout            << "Max abs deviation of divergence across all gridpoints:\t" << max_dev_div[0] <<"\t at ( i=" << max_dev_div[1] <<" , j=" << max_dev_div[2] << " )" << std::endl;
	output_file_test_log << "Max abs deviation of divergence across all gridpoints:\t" << max_dev_div[0] <<"\t at ( i=" << max_dev_div[1] <<" , j=" << max_dev_div[2] << " )\n";
	
	std::cout            << "Stdev of curl across all gridpoints:\t" << stdev_curl_r <<"\t"<< stdev_curl_t <<"\t"<< stdev_curl_p << std::endl;
	output_file_test_log << "Stdev of curl across all gridpoints:\t" << stdev_curl_r <<"\t"<< stdev_curl_t <<"\t"<< stdev_curl_p << "\n";
	
	std::cout            << "Max abs deviation of curl across all gridpoints:"<< std::endl;
	std::cout            << "r:\t" << max_dev_curl_r[0] <<"\t at ( i=" << max_dev_curl_r[1] <<" , j=" << max_dev_curl_r[2] << " )" << std::endl;
	std::cout            << "t:\t" << max_dev_curl_t[0] <<"\t at ( i=" << max_dev_curl_t[1] <<" , j=" << max_dev_curl_t[2] << " )" << std::endl;
	std::cout            << "p:\t" << max_dev_curl_p[0] <<"\t at ( i=" << max_dev_curl_p[1] <<" , j=" << max_dev_curl_p[2] << " )" << std::endl;
	output_file_test_log << "Max abs deviation of curl across all gridpoints:\n";
	output_file_test_log << "r:\t" << max_dev_curl_r[0] <<"\t at ( i=" << max_dev_curl_r[1] <<" , j=" << max_dev_curl_r[2] << " )\n";
	output_file_test_log << "t:\t" << max_dev_curl_t[0] <<"\t at ( i=" << max_dev_curl_t[1] <<" , j=" << max_dev_curl_t[2] << " )\n";
	output_file_test_log << "p:\t" << max_dev_curl_p[0] <<"\t at ( i=" << max_dev_curl_p[1] <<" , j=" << max_dev_curl_p[2] << " )\n";


	
	
	
	//----- Code finished -----
	std::cout << "\nCSV file saved:\t" << output_filename_test_values << std::endl;
	std::cout <<   "Log file saved:\t" << output_filename_test_log    << std::endl;
	std::cout << "\nCode finished" << std::endl;
	
	return 0;
	
}




//----- Define functions -----
// These are duplicated for different test functions. Comment-out the blocks of functions that are not being used.

//--- Nonzero-divergence vector (same functional form as static dipole) ---
double B_r_function_test( double r, double t ){
	return cos(t) * pow(r,-3);
}

double B_t_function_test( double r, double t ){
	return 2.0 * sin(t) * pow(r,-3);
}

double B_p_function_test( double r, double t ){
	return sin(2.0*t) * pow(r,-3);
}


double B_div_exact_function( double r, double t ){
	return 3.0 * cos(t) * pow(r,-4);
}


double B_curl_exact_function_r( double r, double t ){
	return ( 2.0 * pow(cos(t),2) - pow(sin(t),2) ) * 2.0 * pow(r,-4);
}

double B_curl_exact_function_t( double r, double t ){
	return 2.0 * sin( 2.0 * t )  * pow(r,-4);
}

double B_curl_exact_function_p( double r, double t ){
	return -3.0 * sin(t) * pow(r,-4);
}


/*
//--- Divergenceless vector (static dipole with A(r,theta) = r^3 sin^2(theta) ) ---
double B_r_function_test( double r, double t ){
	return 2.0 * cos(t) * pow(r,-3);
}

double B_t_function_test( double r, double t ){
	return sin(t) * pow(r,-3);
}

double B_p_function_test( double r, double t ){
	return sin(t) * r*r;
}


double B_div_exact_function( double r, double t ){
	return 0;
}


double B_curl_exact_function_r( double r, double t ){
	return 2.0 * r * cos(t);
}

double B_curl_exact_function_t( double r, double t ){
	return -3.0 * r * sin(t);
}

double B_curl_exact_function_p( double r, double t ){
	return 0;
}
*/