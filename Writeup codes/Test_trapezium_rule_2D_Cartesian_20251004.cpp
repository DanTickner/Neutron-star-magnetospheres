/*
cd OneDrive\PhD\20250915 Updated codes
g++ Test_trapezium_rule_2D_Cartesian_20251004.cpp -o Test_trapezium_rule_2D_Cartesian_20251004 --std=c++11
Test_trapezium_rule_2D_Cartesian_20251004

Test the 2D trapezium rule for a known integral in Cartesian coordinates.
For simplicity, we only use a linear spacing, but the code and functions are written such that any spacing could be used.

*/


#include <iostream>
#include <math.h>
#include <fstream>
#include "Trapezium_rule.h"


//----- Declare functions -----
double f_function( double u1, double u2 );
double integral_indefinite_function( double u1, double u2 );



//----- Main code -----
int main(){
	
	//----- Define variables -----
	std::string output_filename_csv = "CSV/20251004_Test_trapezium_rule_2D_Cartesian.csv";
	std::ofstream output_file_csv;
	
	double u1_min      = 1.0;
	double u1_max      = 2.0;
	double u2_min      = 3.5;
	double u2_max      = 5.0;
	
	
	std::vector<double> N_values_u1_vector { 5, 10, 20, 40 };
	std::vector<double> N_values_u2_vector { 5, 10, 20, 40 };
	
	
	
	
	//----- Output coordinates, function values and finite area elements to the CSV -----
	output_file_csv.open( output_filename_csv );
	
	output_file_csv <<  "u1_min," << u1_min
	                << ",u1_max," << u1_max
					<< ",u2_min," << u2_min
					<< ",u2_max," << u2_max << "\n";
	
	output_file_csv << "N1,N2,integral_numerical,absolute_error,absolute_relative_error\n";
	std::cout << "N1\tN2\tnumerical\tabs_err\tabs_rel_err" << std::endl;
	
	
	
	
	//----- Loop through integrals with varying numbers of coordinates -----
	for( int N1=0; N1<N_values_u1_vector.size(); N1++ ){
		for( int N2=0; N2<N_values_u2_vector.size(); N2++ ){
			
			
			int    N_values_u1 = N_values_u1_vector[N1];
			int    N_values_u2 = N_values_u2_vector[N2];
			
			std::vector<double> u1 ( N_values_u1 );
			std::vector<double> u2 ( N_values_u2 );
			std::vector< std::vector<double> > f  ( N_values_u1, std::vector<double> ( N_values_u2 ) );
			std::vector< std::vector<double> > dA ( N_values_u1, std::vector<double> ( N_values_u2 ) );
			
			
			//----- Split the intervals and calculate the coordinates -----
			double delta_u1 = ( u1_max - u1_min ) / ( (double) N_values_u1 - 1.0 );
			double delta_u2 = ( u2_max - u2_min ) / ( (double) N_values_u2 - 1.0 );
			
			for( int i=0; i<N_values_u1; i++ ){
				u1[i] = u1_min + (double) i * delta_u1;
			}
			for( int j=0; j<N_values_u2; j++ ){
				u2[j] = u2_min + (double) j * delta_u2;
			}
			
			
			//----- Calculate function values -----
			for( int i=0; i<N_values_u1; i++ ){
				for( int j=0; j<N_values_u2; j++ ){
					f[i][j] = f_function( u1[i], u2[j] );
				}
			}
			
			
			//----- Calculate finite area element values -----
			// dA[i][j] depends on i+1 and j+1 so must be done after all function values calculated.
			for( int i=0; i<N_values_u1; i++ ){
				for( int j=0; j<N_values_u2; j++ ){
					dA[i][j] = ( u1[i+1] - u1[i] ) * ( u2[j+1] - u2[j] );		// This will need to change for different coordinate systems.
				}
			}
			
			
			//----- Calculate the integral -----
			//double integral_numerical = integral_trapezium_2D_with_cout( f, dA );
			double integral_numerical = integral_trapezium_2D( f, dA );
			double integral_exact     = integral_indefinite_function( u1_max, u2_max ) + integral_indefinite_function( u1_min, u2_min ) - integral_indefinite_function( u1_min, u2_max ) - integral_indefinite_function( u1_max, u2_min );
			
			double abs_err = abs( integral_numerical - integral_exact );
			double abs_rel_err = abs( abs_err / integral_exact );
			
			std::cout       << N_values_u1 <<"\t"<< N_values_u2 <<"\t"<< integral_numerical <<"\t"<< abs_err <<"\t"<< abs_rel_err << std::endl;
			output_file_csv << N_values_u1 <<"," << N_values_u2 <<"," << integral_numerical <<"," << abs_err <<"," << abs_rel_err << "\n";
			
		}
		std::cout << std::endl;
	}
	
	
	
	
	//----- Code finished -----
	std::cout << "\nLog file saved:\t" << output_filename_csv << std::endl;
	std::cout << "\nCode finished" << std::endl;
	
	return 0;
	
}


//----- Define functions -----
double f_function( double u1, double u2 ){
	return u1 * u1 + u2 * u2;
}

double integral_indefinite_function( double u1, double u2 ){
	return (1.0/3.0) * ( pow(u1,3) * u2 + u1 * pow(u2,3) );
}