/*
cd OneDrive\PhD\Codes\20240225 Multidimensional trapezium rule redux
g++ Test_2D_polar_20240225.cpp -o Test_2D_polar_20240225 --std=c++11
Test_2D_polar_20240225

Test the 2D trapezium rule for a known integral in circular polar coordinates.
For simplicity, we only use a linear spacing, but the code and functions are written such that any spacing could be used.

*/


#include <iostream>
#include <math.h>
#include "Trapezium_rule.h"


//----- Declare functions -----
double f_function( double u1, double u2 );
double integral_indefinite_function( double u1, double u2 );



//----- Main code -----
int main(){
	
	//----- Define variables -----
	const double pi    = acos( -1 );
	double u1_min      = pi / 3.0;
	double u1_max      = pi / 2.0;
	double u2_min      = pi / 4.0;
	double u2_max      = pi;
	int    N_values_u1 = 10;
	int    N_values_u2 = 20;
	
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
			dA[i][j] = 0.5 * ( pow(u1[i+1],2) - pow(u1[i],2) ) * ( u2[j+1] - u2[j] );		// This will need to change for different coordinate systems.
		}
	}
	
	
	//----- Output coordinates, function values and finite area elements to the screen -----
	std::cout << "i\tj\tu1[i]\tu2[j]\tf[i][j]\tdA[i][j]" << std::endl;
	
	for( int i=0; i<N_values_u1; i++ ){
		for( int j=0; j<N_values_u2; j++ ){
			std::cout << i <<"\t"<< j <<"\t"<< u1[i] <<"\t"<< u2[j] <<"\t"<< f[i][j] <<"\t"<< dA[i][j] << std::endl;
		}
	}
	
	
	//----- Calculate the integral -----
	//double integral_numerical = integral_trapezium_2D_with_cout( f, dA );
	double integral_numerical = integral_trapezium_2D( f, dA );
	double integral_exact     = integral_indefinite_function( u1_max, u2_max ) + integral_indefinite_function( u1_min, u2_min ) - integral_indefinite_function( u1_min, u2_max ) - integral_indefinite_function( u1_max, u2_min );
	
	std::cout << "\nNumerical:\t" << integral_numerical << std::endl;
	std::cout <<   "Exact    :\t" << integral_exact     << std::endl;
	
	
	//----- Code finished -----
	std::cout << "\nCode finished" << std::endl;
	return 0;
	
}


//----- Define functions -----
double f_function( double u1, double u2 ){
	return u1 * cos( u2 );
}

double integral_indefinite_function( double u1, double u2 ){
	return (1.0/3.0) * pow(u1,3) * sin( u2 );
}