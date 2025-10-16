/*
cd OneDrive\PhD\Codes\20240311 Time evolution code\Test codes
g++ Test_trapezium_rule_1D_Cartesian_20251004.cpp -o Test_trapezium_rule_1D_Cartesian_20251004 --std=c++11
Test_trapezium_rule_1D_Cartesian_20251004

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
	std::string output_filename_csv = "CSV/20251004_Test_trapezium_rule_1D_Cartesian.csv";
	std::ofstream output_file_csv;
	
	double x_min      = 1.0;
	double x_max      = 2.0;
	
	std::vector<double> N_values_x_vector { 5, 10, 20, 40 };
	
	
	//----- Output coordinates, function values and finite area elements to the CSV -----
	output_file_csv.open( output_filename_csv );
	
	output_file_csv <<  "x_min," << x_min
	                << ",x_max," << x_max << "\n";
	
	output_file_csv << "Nx,integral_numerical,absolute_error,absolute_relative_error\n";
	std::cout << "Nx\tnumerical\tabs_err\tabs_rel_err" << std::endl;
	
	
	//----- Loop through integrals with varying numbers of coordinates -----
	for( int N=0; N<N_values_x_vector.size(); N++ ){
		
		
		
		int    N_values_x = N_values_x_vector[N];
		
		std::vector<double> x;
		std::vector<double> f;
		
		
		//----- Divide the interval into subintervals and calculate values of x and f(x) -----
		double delta_x = ( x_max - x_min ) / ( (double) N_values_x - 1.0 );
		
		for( int i=0; i<N_values_x; i++ ){
			x.push_back( x_min + (double) i * delta_x );
			f.push_back( f_function( x[i] ) );
		}
		
		
		//----- Calculate the integral -----
		//double integral_numerical = integral_trapezium_1D_Cartesian_with_cout( f, x );
		double integral_numerical = integral_trapezium_1D_Cartesian( f, x );
		//double integral_numerical = integral_trapezium_1D_Cartesian_efficient( f, x );
		//double integral_numerical = integral_trapezium_1D_Cartesian_constant_spacing( f, delta_x );
		//double integral_numerical = integral_trapezium_1D_Cartesian_constant_spacing_efficient( f, delta_x );
		
		double integral_exact     = integral_indefinite_function( x_max ) - integral_indefinite_function( x_min );
		
		double abs_err = abs( integral_numerical - integral_exact );
		double abs_rel_err = abs( abs_err / integral_exact );
		
		std::cout       << N_values_x <<"\t"<< integral_numerical <<"\t"<< abs_err <<"\t"<< abs_rel_err << std::endl;
		output_file_csv << N_values_x <<"," << integral_numerical <<"," << abs_err <<"," << abs_rel_err << "\n";
		
	}
	
	
	
	
	//----- Code finished -----
	std::cout << "\nLog file saved:\t" << output_filename_csv << std::endl;
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