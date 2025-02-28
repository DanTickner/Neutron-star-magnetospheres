/*
cd OneDrive\PhD\Codes\20240311 Time evolution code\Test codes
g++ Test_trapezium_rule_1D_Cartesian.cpp -o Test_trapezium_rule_1D_Cartesian --std=c++11
Test_trapezium_rule_1D_Cartesian

Test the 1D trapezium rule for a known integral.
For simplicity, we only use a linear spacing, but the code and functions are written such that any spacing could be used.

*/


#include <iostream>
#include <math.h>
#include <fstream>
#include "Trapezium_rule.h"


//----- Declare functions -----
double f_function( double x );
double integral_indefinite_function( double x );



//----- Main code -----
int main(){
	
	//----- Define variables -----
	std::string output_filename_test_log = "../Logs/20240227_Test_trapezium_rule_1D_Cartesian.txt";
	std::ofstream output_file_test_log;
	
	
	double x_min      = 1.0;
	double x_max      = 2.0;
	int    N_values_x = 10;
	
	std::vector<double> x;
	std::vector<double> f;
	
	
	//----- Divide the interval into subintervals and calculate values of x and f(x) -----
	double delta_x = ( x_max - x_min ) / ( (double) N_values_x - 1.0 );
	
	output_file_test_log.open( output_filename_test_log );
	
	std::cout            << "i\tx[i]\tf[i]" << std::endl;
	output_file_test_log << "i\tx[i]\tf[i]" << "\n";
	
	for( int i=0; i<N_values_x; i++ ){
		x.push_back( x_min + (double) i * delta_x );
		f.push_back( f_function( x[i] ) );
		std::cout            << i <<"\t"<< x[i] <<"\t"<< f[i] << std::endl;
		output_file_test_log << i <<"\t"<< x[i] <<"\t"<< f[i] << "\n";
	}
	
	
	//----- Calculate the integral -----
	//double integral_numerical = integral_trapezium_1D_Cartesian_with_cout( f, x );
	double integral_numerical = integral_trapezium_1D_Cartesian( f, x );
	//double integral_numerical = integral_trapezium_1D_Cartesian_efficient( f, x );
	//double integral_numerical = integral_trapezium_1D_Cartesian_constant_spacing( f, delta_x );
	//double integral_numerical = integral_trapezium_1D_Cartesian_constant_spacing_efficient( f, delta_x );
	
	double integral_exact     = integral_indefinite_function( x_max ) - integral_indefinite_function( x_min );
	
	std::cout            << "\nNumerical:\t" << integral_numerical << std::endl;
	std::cout            <<   "Exact    :\t" << integral_exact     << std::endl;
	
	output_file_test_log << "\nNumerical:\t" << integral_numerical << "\n";
	output_file_test_log <<   "Exact    :\t" << integral_exact     << "\n";
	
	std::cout << "\nLog file saved:\t" << output_filename_test_log << std::endl;
	
	
	//----- Code finished -----
	std::cout << "\nCode finished" << std::endl;
	return 0;
	
}


//----- Define functions -----
double f_function( double x ){
	return x * x;
}

double integral_indefinite_function( double x ){
	return pow( x, 3 ) / 3.0;
}