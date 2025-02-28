/*
cd OneDrive\PhD\Codes\20240311 Time evolution code\Test codes
g++ Test_Finite_Difference_2D.cpp -o Test_Finite_Difference_2D --std=c++11
Test_Finite_Difference_2D

Test the implementation of combinations of finite difference schemes of 2D functions so that all gridpoints have an evaluated first derivative up to required accuracy.
*/

#include <iostream>
#include <iomanip>
#include <vector>
#include <math.h>
#include <fstream>
#include <chrono>

#include "Finite_Difference_Expressions.h"


//----- Declare functions -----
double f_function         ( double x, double y );
double f_dx_exact_function( double x, double y );
double f_dy_exact_function( double x, double y );
double stdev_of_array_2d( std::vector< std::vector<double> > exact, std::vector< std::vector<double> > calc );



//----- Main code -----
int main(){
	
	
	//----- Extra variables useful only for this testing code -----
	
	double      x_0     = 2.0;	// Minimum gridpoint.
	double      y_0     = x_0;//3.0;
	double      h_x     = 0.1;	// Constant grid spacing. Probably doesn't matter but keep positive.
	double      h_y     = h_x;//0.2;
	int         N_x     = 15;	// Number of gridpoints to use for the entire domain. Keep low so that the endpoints are easily diagnosed.
	int         N_y     = N_x;//20;	// Number of gridpoints to use for the entire domain. Keep low so that the endpoints are easily diagnosed.
	std::string h_order = "10_with_decimal";//"6_with_8_near_endpoints";		// The code will chose the scheme with accuracy O( h^N ). This is the N. NOT YET IMPLEMENTED.
	
	std::string output_filename_test_base = "20240517_Test_finite_difference_2D";	// Don't add ".csv" as the code will append _r, _t, etc to the filename for each csv.
	
	std::string output_filename_test_csv = "../CSV/"  + output_filename_test_base + "_" + h_order + ".csv";
	
	std::ofstream output_file_test_csv;
	
	std::vector<double> x ( N_x );
	std::vector<double> y ( N_y );
	
	std::vector< std::vector<double> > f          ( N_x, std::vector<double> ( N_y ) );
	std::vector< std::vector<double> > f_dx_exact ( N_x, std::vector<double> ( N_y ) );
	std::vector< std::vector<double> > f_dy_exact ( N_x, std::vector<double> ( N_y ) );
	std::vector< std::vector<double> > f_dx_FD    ( N_x, std::vector<double> ( N_y ) );	// Finite-difference expression.
	std::vector< std::vector<double> > f_dy_FD    ( N_x, std::vector<double> ( N_y ) );	// Finite-difference epxression.
	std::vector< std::vector<double> > f_dx_ARE   ( N_x, std::vector<double> ( N_y ) );	// Absolute relative error of finite-difference derivative at each gridpoint.
	std::vector< std::vector<double> > f_dy_ARE   ( N_x, std::vector<double> ( N_y ) );	// Absolute relative error of finite-difference derivative at each gridpoint.
	
	
	//----- Build vectors of gridpoints and exact values -----
	for( int i=0; i<N_x; i++ ){
		x[i] = x_0 + i * h_x;
	}
	
	for( int j=0; j<N_y; j++ ){
		y[j] = y_0 + j * h_y;
	}
	
	for( int i=0; i<N_x; i++ ){
		for( int j=0; j<N_y; j++ ){
			f         [i][j] = f_function         ( x[i], y[j] );
			f_dx_exact[i][j] = f_dx_exact_function( x[i], y[j] );
			f_dy_exact[i][j] = f_dy_exact_function( x[i], y[j] );
		}
	}
	
	
	//----- Calculate finite-difference expressions -----
	if( h_order == "1"){
		f_dx_FD = f_dx_FD_function_order_1( f, h_x );
		f_dy_FD = f_dy_FD_function_order_1( f, h_y );
	}
	else if( h_order == "6_with_8_near_endpoints" ){
		f_dx_FD = f_dx_FD_function_order_6_with_8_near_endpoints( f, h_x );
		f_dy_FD = f_dy_FD_function_order_6_with_8_near_endpoints( f, h_y );
	}
	else if( h_order == "10" ){
		f_dx_FD = f_dx_FD_function_order_10( f, h_x );
		f_dy_FD = f_dy_FD_function_order_10( f, h_y );
	}
	else if( h_order == "10_with_decimal" ){
		f_dx_FD = f_dx_FD_function_order_10_decimal( f, h_x );
		f_dy_FD = f_dy_FD_function_order_10_decimal( f, h_y );
	}
	else {
		std::cout << "h_order string " << h_order << " not recognised. Stopping the code." << std::endl;
		return 1;
	}
	
	for( int i=0; i<N_x; i++ ){
		for( int j=0; j<N_y; j++ ){
			f_dx_ARE[i][j] = ( f_dx_exact[i][j] != 0 ) ? ( abs( 1.0 - f_dx_FD[i][j] / f_dx_exact[i][j] ) ) : 0.0;
			f_dy_ARE[i][j] = ( f_dy_exact[i][j] != 0 ) ? ( abs( 1.0 - f_dy_FD[i][j] / f_dy_exact[i][j] ) ) : 0.0;
		}
	}
	
	
	//----- Output to CSV -----
	output_file_test_csv.open( output_filename_test_csv );
	output_file_test_csv << "i,j,f_dx_exact,f_dx_FD,f_dx_ARE,f_dy_exact,f_dy_FD,f_dy_ARE\n";
	
	for( int i=0; i<N_x; i++ ){
		for( int j=0; j<N_y; j++ ){
			output_file_test_csv << i <<","<< j <<","<< f_dx_exact[i][j] <<","<< f_dx_FD[i][j] <<","<< f_dx_ARE[i][j]
			                                    <<","<< f_dy_exact[i][j] <<","<< f_dy_FD[i][j] <<","<< f_dy_ARE[i][j] << "\n";
		}
	}
	
	
	
	//----- Code finished -----
	std::cout << "\nStdev of df/dx:\t" << stdev_of_array_2d( f_dx_exact, f_dx_FD ) << std::endl;
	std::cout <<   "Stdev of df/dy:\t" << stdev_of_array_2d( f_dy_exact, f_dy_FD ) << std::endl;
	
	std::cout << "\nCSV file saved:\t" << output_filename_test_csv << std::endl;
	
	std::cout << "\nCode finished" << std::endl;
	
	return 0;
	
}


//----- Define functions -----
double f_function         ( double x, double y ){
	//return pow( x, -3.0 ) * pow( y, 2.5 );
	return pow( x, -3.0 ) * pow( y, -3.0 );
}


double f_dx_exact_function( double x, double y ){
	//return -3.0 * pow( x, -4.0 ) * pow( y, 2.5 );
	return -3.0 * pow( x, -4.0 ) * pow( y, -3.0 );
}

double f_dy_exact_function( double x, double y ){
	//return pow( x, -3.0 ) * 2.5 * pow( y, 1.5 );
	return pow( x, -3.0 ) * -3.0 * pow( y, -4.0 );
}

double stdev_of_array_2d( std::vector< std::vector<double> > exact, std::vector< std::vector<double> > calc ){
	// Calculate standard deviation of an array of calculated results, compared to their known exact values.
	
	double ret = 0;
	
	for( int i=0; i<exact.size(); i++ ){
		for( int j=0; j<exact[0].size(); j++ ){
			ret += pow( exact[i][j] - calc[i][j], 2.0 );
		}
	}
	
	return sqrt( ret / ( (double) exact.size() * exact[0].size() ) );
	
}