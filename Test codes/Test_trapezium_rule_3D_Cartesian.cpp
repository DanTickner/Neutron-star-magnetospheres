/*
cd OneDrive\PhD\Codes\20230810 Time evolution with updated equations\Test codes
g++ Test_trapezium_rule_3D_Cartesian.cpp -o Test_trapezium_rule_3D_Cartesian --std=c++11
Test_trapezium_rule_3D_Cartesian

Test the 3D trapezium rule for a known integral in Cartesian coordinates.
For simplicity, we only use a linear spacing, but the code and functions are written such that any spacing could be used.

*/


#include <iostream>
#include <math.h>
#include <fstream>
#include "Trapezium_rule.h"


//----- Declare functions -----
double f_function( double u1, double u2, double u3 );
double integral_indefinite_function( double u1, double u2, double u3 );



//----- Main code -----
int main(){
	
	//----- Define variables -----
	std::string output_filename_test_log = "../Logs/20240227_Test_trapezium_rule_3D_Cartesian.txt";
	std::ofstream output_file_test_log;
	
	double u1_min      = 1.0;
	double u1_max      = 2.0;
	double u2_min      = 3.5;
	double u2_max      = 5.0;
	double u3_min      = 2.5;
	double u3_max      = 6.0;
	int    N_values_u1 = 10;
	int    N_values_u2 = 20;
	int    N_values_u3 = 30;
	
	std::vector<double> u1 ( N_values_u1 );
	std::vector<double> u2 ( N_values_u2 );
	std::vector<double> u3 ( N_values_u3 );
	std::vector< std::vector< std::vector< double> > > f  ( N_values_u1, std::vector< std::vector<double> > ( N_values_u2, std::vector<double> ( N_values_u3 ) ) );
	std::vector< std::vector< std::vector< double> > > dV ( N_values_u1, std::vector< std::vector<double> > ( N_values_u2, std::vector<double> ( N_values_u3 ) ) );
	
	
	//----- Split the intervals and calculate the coordinates -----
	double delta_u1 = ( u1_max - u1_min ) / ( (double) N_values_u1 - 1.0 );
	double delta_u2 = ( u2_max - u2_min ) / ( (double) N_values_u2 - 1.0 );
	double delta_u3 = ( u3_max - u3_min ) / ( (double) N_values_u3 - 1.0 );
	
	for( int i=0; i<N_values_u1; i++ ){
		u1[i] = u1_min + (double) i * delta_u1;
	}
	for( int j=0; j<N_values_u2; j++ ){
		u2[j] = u2_min + (double) j * delta_u2;
	}
	for( int k=0; k<N_values_u3; k++ ){
		u3[k] = u3_min + (double) k * delta_u3;
	}
	
	
	//----- Calculate function values -----
	for( int i=0; i<N_values_u1; i++ ){
		for( int j=0; j<N_values_u2; j++ ){
			for( int k=0; k<N_values_u3; k++ ){
				f[i][j][k] = f_function( u1[i], u2[j], u3[k] );
			}
		}
	}
	
	
	//----- Calculate finite volume element values -----
	// dV[i][j][k] depends on i+1, j+1 and k+1 so must be done after all function values calculated.
	for( int i=0; i<N_values_u1; i++ ){
		for( int j=0; j<N_values_u2; j++ ){
			for( int k=0; k<N_values_u3; k++ ){
				dV[i][j][k] = ( u1[i+1] - u1[i] ) * ( u2[j+1] - u2[j] ) * ( u3[k+1] - u3[k] );		// This will need to change for different coordinate systems.
			}
		}
	}
	
	
	//----- Output coordinates, function values and finite volume elements to the screen -----
	output_file_test_log.open( output_filename_test_log );
	
	std::cout            << "i\tj\tk\tu1[i]\tu2[j]\tu3[k]\tf[i][j][k]\tdV[i][j][k]\tDelta x*Delta y*Delta z" << std::endl;
	output_file_test_log << "i\tj\tk\tu1[i]\tu2[j]\tu3[k]\tf[i][j][k]\tdV[i][j][k]\tDelta x*Delta y*Delta z" << "\n";
	
	for( int i=0; i<N_values_u1; i++ ){
		for( int j=0; j<N_values_u2; j++ ){
			for( int k=0; k<N_values_u3; k++ ){
				std::cout            << i <<"\t"<< j <<"\t"<< k <<"\t"<< u1[i] <<"\t"<< u2[j] <<"\t"<< u3[k] <<"\t"<< f[i][j][k] <<"\t"<< dV[i][j][k] <<"\t"<< delta_u1 * delta_u2 * delta_u3 << std::endl;
				output_file_test_log << i <<"\t"<< j <<"\t"<< k <<"\t"<< u1[i] <<"\t"<< u2[j] <<"\t"<< u3[k] <<"\t"<< f[i][j][k] <<"\t"<< dV[i][j][k] <<"\t"<< delta_u1 * delta_u2 * delta_u3 << "\n";
			}
		}
	}
	
	
	//----- Calculate the integral -----
	double integral_numerical =   integral_trapezium_3D( f, dV );
	double integral_exact     =   integral_indefinite_function( u1_max, u2_max, u3_max ) - integral_indefinite_function( u1_max, u2_max, u3_min ) - integral_indefinite_function( u1_max, u2_min, u3_max ) + integral_indefinite_function( u1_max, u2_min, u3_min )
	                            - integral_indefinite_function( u1_min, u2_max, u3_max ) + integral_indefinite_function( u1_min, u2_max, u3_min ) + integral_indefinite_function( u1_min, u2_min, u3_max ) - integral_indefinite_function( u1_min, u2_min, u3_min );
	
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
double f_function( double u1, double u2, double u3 ){
	return u1 * u1 + u2 * u2 + u3 * u3;
}

double integral_indefinite_function( double u1, double u2, double u3 ){
	return (1.0/3.0) * ( pow(u1,3) * u2 * u3 + u1 * pow(u2,3) * u3 + u1 * u2 * pow(u3,3) );
}