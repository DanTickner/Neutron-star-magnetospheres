/*
cd OneDrive\PhD\Codes\20240311 Time evolution code\Test codes
g++ Test_Finite_Difference_Second_Derivative.cpp -o Test_Finite_Difference_Second_Derivative --std=c++11
Test_Finite_Difference_Second_Derivative

Test the implementation of combinations of finite difference schemes so that all gridpoints have an evaluated first derivative up to required accuracy.
*/

#include <iostream>
#include <iomanip>
#include <vector>
#include <math.h>
#include <fstream>
#include <chrono>

#include "Finite_Difference_Expressions.h"


//----- Declare functions -----
double f_function         ( double x );
double f_dx2_exact_function( double x );



//----- Main code -----
int main(){
	
	
	//----- Extra variables useful only for this testing code -----
	
	double x_0     = 2.0;	// Minimum gridpoint.
	double h       = 0.1;	// Constant grid spacing. Probably doesn't matter but keep positive.
	int    N       = 10;	// Number of gridpoints to use for the entire domain. Keep low so that the endpoints are easily diagnosed.
	int    h_order = 6;		// The code will chose the scheme with accuracy O( h^N ). This is the N.
	
	std::string output_filename_test_base = "Test_finite_difference_second_derivative_20240501";	// Don't add ".csv" as the code will append _r, _t, etc to the filename for each csv.
	
	std::string output_filename_test_csv = "../CSV/"  + output_filename_test_base + "_" + std::to_string( h_order ) + ".csv";
	
	std::ofstream output_file_test_csv;
	
	std::vector<double> x                      ( N );
	std::vector<double> f                      ( N );
	std::vector<double> f_dx2_exact            ( N );
	std::vector<double> f_dx2_FD               ( N );
	std::vector<double> f_dx2_ARE              ( N );	// Absolute relative error of finite-difference derivative at each gridpoint.
	std::vector<double> f_dx2_abs_err          ( N );	// Absolute error.
	std::vector<double> f_dx2_abs_err_expected ( N );	// Expected absolute error ( equal to h^N ).
	
	
	//----- Build vectors of gridpoints and exact values -----
	for( int n=0; n<N; n++ ){
		x[n] = x_0 + n * h;
		f[n] = f_function( x[n] );
		f_dx2_exact[n] = f_dx2_exact_function( x[n] );
	}
	
	
	//----- Calculate derivatives to various accuracies -----
	std::cout << "Calculating second derivatives to order ( h^" << h_order << " )." << std::endl;
	
	switch( h_order ){
		
		case 1:
			f_dx2_FD = f_dx2_FD_function_order_1( f, h );
			break;
		case 2:
			f_dx2_FD = f_dx2_FD_function_order_2( f, h );
			break;
		case 3:
			f_dx2_FD = f_dx2_FD_function_order_3( f, h );
			break;
		case 4:
			f_dx2_FD = f_dx2_FD_function_order_4( f, h );
			break;
		case 5:
			f_dx2_FD = f_dx2_FD_function_order_5( f, h );
			break;
		case 6:
			f_dx2_FD = f_dx2_FD_function_order_6( f, h );
			break;
		default:
			std::cout << "No finite-difference method has been defined for order " << h_order << ". Stopping the code." << std::endl;
			return 1;
	}
	
	for( int n=0; n<N; n++ ){
		f_dx2_ARE             [n] = ( f_dx2_exact[n] != 0 ) ? ( abs( 1.0 - f_dx2_FD[n] / f_dx2_exact[n] ) ) : 0.0;
		f_dx2_abs_err         [n] = abs( f_dx2_FD[n] - f_dx2_exact[n] );
		f_dx2_abs_err_expected[n] = pow( h, h_order );
	}
	
	
	
	//----- Output to CSV -----
	output_file_test_csv.open( output_filename_test_csv );
	output_file_test_csv << "n,f_dx2_exact,f_dx2_FD,ARE,absolute_error,abs_err_expected\n";
	
	for( int n=0; n<N; n++ ){
		output_file_test_csv << n <<","<< f_dx2_exact[n] <<","<< f_dx2_FD[n] <<","<< f_dx2_ARE[n] <<","<< f_dx2_abs_err[n] <<","<< f_dx2_abs_err_expected[n] << "\n";
	}
	
	
	
	//----- Code finished -----
	
	std::cout << "\nCSV file saved:\t" << output_filename_test_csv << std::endl;
	
	std::cout << "\nCode finished" << std::endl;
	
	return 0;
	
}


//----- Define functions -----
double f_function         ( double x ){
	return pow( x, -3.0 );
}


double f_dx2_exact_function( double x ){
	return 12.0 * pow( x, -5.0 );
}