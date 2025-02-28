/*

cd OneDrive\PhD\Codes\VSH series
g++ AVSH_series_from_CSV_01.cpp -o AVSH_series_from_CSV_01 --std=c++11
AVSH_series_from_CSV_01

Read a vector function 
A(r,theta) = A_r(r,theta) e_r + A_theta(r,theta) e_theta + A_phi(r,theta) e_phi
from a CSV file of discrete datapoints and calculate its AVSH (axisymmetric vector spherical harmonic) series coefficients at a given value of r.

CSV reading code lines taken from ../Helloworld codes/Read_CSV.cpp.

Based on VSH_series_01.cpp.

*/


#include <iostream>
#include <iomanip>
#include <math.h>
#include <vector>
#include <chrono>
#include <fstream>
#include <sstream>
#include <string>

#include "../Generate_Associated_Legendre_Functions.h"	// Includes "../Vector_Operations.h", which we also need here.
#include "../Numerical_Integration.h"

const double pi = acos(-1);




int main(){
	
	//----- Define variables -----
	int ell_max = 11;
	int w = 18;							// Fixed width for screen output.
	
	double test_r     = 1.0;//0.321;				// Value of x at which to test the series.
	double test_theta = 0.432;
	
	std::string                             filename   = "CSV/20230925.csv";
	std::vector< std::vector<std::string> > content;
	std::vector<             std::string  > row;
	std::string                             line, word;
	
	
	std::vector<double> theta_values;
	
	std::vector<double> integrand_values_r;
	std::vector<double> integrand_values_1;
	std::vector<double> integrand_values_2;
	
	std::vector<double> coeffs_r ( ell_max+1 );
	std::vector<double> coeffs_1 ( ell_max+1 );
	std::vector<double> coeffs_2 ( ell_max+1 );
	
	std::vector<double> A_r_values;
	std::vector<double> A_theta_values;
	std::vector<double> A_phi_values;
	
	

	
	//----- Start timing -----
	auto time_start = std::chrono::high_resolution_clock::now();	
	
	
	
	
	//----- Read the CSV file -----
	std::fstream input_file ( filename, std::ios::in );
	
	if( input_file.is_open() ){
		
		while( std::getline( input_file, line ) ){
			
			row.clear();
			
			std::stringstream str( line );
			
			while( std::getline( str, word, ',' ) ){	// MUST use single quotes ',' and not double quotes ","
				row.push_back( word );
			}
			
			content.push_back( row  );
		}
	}
	
	else std::cout << "Could not open the file " << filename << std::endl;
	
	
	
	
	//----- Build lists of coordinates and precalculate vector values -----
	for( int i=1; i<content.size(); i++ ){
		theta_values  .push_back( std::stod( content[i][0] ) );
		A_r_values    .push_back( std::stod( content[i][1] ) );
		A_theta_values.push_back( std::stod( content[i][2] ) );
		A_phi_values  .push_back( std::stod( content[i][3] ) );
	}
	
	
	
	
	//----- Calculate VSH series coefficients -----
	
	integrand_values_r = std::vector<double> ( theta_values.size() );
	integrand_values_1 = std::vector<double> ( theta_values.size() );
	integrand_values_2 = std::vector<double> ( theta_values.size() );
	
	std::cout << "ell" << "\t"<< "0" <<"\t|\t"<< std::left<<std::setw(w)<< "A^{r,ell}_m" <<std::left<<std::setw(w)<< "A^{(1),ell}_m" <<std::left<<std::setw(w)<< "A^{(2),ell}_m" << std::endl;
	
	for( int ell=0; ell<=ell_max; ell++ ){
		
		//--- Build lists of integrand values ---
		for( int i=0; i<theta_values.size(); i++ ){
				
			double theta = theta_values[i];
			
			double P_ell_0 = p_ell  ( cos(theta), ell );
			double P_ell_1 = p_ell_1( cos(theta), ell );
			
			integrand_values_r[i] = A_r_values    [i] * P_ell_0 * sin( theta );
			integrand_values_1[i] = A_theta_values[i] * P_ell_1 * sin( theta );
			integrand_values_2[i] = A_phi_values  [i] * P_ell_1 * sin( theta );
			
		}
		
		
		//--- Integrate and put values into array ---
		coeffs_r[ell] = trapezium( integrand_values_r, theta_values ) * sqrt( (2.0*ell+1.0)*pi );
		coeffs_1[ell] = trapezium( integrand_values_1, theta_values ) * sqrt( (2.0*ell+1.0)*pi ) / ( (double)ell*(ell+1.0) );
		coeffs_2[ell] = trapezium( integrand_values_2, theta_values ) * sqrt( (2.0*ell+1.0)*pi ) / ( (double)ell*(ell+1.0) );
		
		
		if( ell == 0 ){
			coeffs_1[0] = 0;
			coeffs_2[0] = 0;
		}
		
		std::cout << ell << "\t"<< 0 <<"\t|\t"<< std::left<<std::setw(w)<< coeffs_r[ell] <<std::left<<std::setw(w)<< coeffs_1[ell] <<std::left<<std::setw(w)<< coeffs_2[ell] << std::endl;
	}
	
	std::cout << integrand_values_r.size() << std::endl;
	
	
	
	
	//----- Evaluate VSH series at the chosen gridpoint -----
	
	std::vector< std::complex<double> > A_series { 0, 0, 0 };
	
	for( int ell=0; ell<=ell_max; ell++ ){
		
		double sqrt_factor = sqrt( (2.0*ell+1.0) / ( 4.0 * pi ) );
		
		double P_ell_0 = p_ell  ( cos(test_theta), ell );
		double P_ell_1 = p_ell_1( cos(test_theta), ell );
		
		A_series[0] += sqrt_factor * P_ell_0 * coeffs_r[ell];
		A_series[1] += sqrt_factor * P_ell_1 * coeffs_1[ell];
		A_series[2] += sqrt_factor * P_ell_1 * coeffs_2[ell];
		
	}
	
	std::cout << std::endl;
	cout_vector_setw( A_series, w, "A_series" );
	
	
	
	
	//----- Output execution time and finish code -----
	auto   time_stop   = std::chrono::high_resolution_clock::now();
	double exec_time   = std::chrono::duration_cast<std::chrono::nanoseconds>( time_stop - time_start ).count() * 1e-9;
	int    exec_time_h = exec_time / 3600;
	int    exec_time_m = exec_time / 60 - exec_time_h * 60;
	int    exec_time_s = exec_time - exec_time_h * 3600 - exec_time_m * 60;
	std::cout << "\nExecution time (h:mm:ss):\t" << exec_time_h;
	if( exec_time_m < 10 ){
		std::cout << ":0" << exec_time_m;
	} else {
		std::cout << ":"  << exec_time_m;
	}
	if( exec_time_s < 10 ){
		std::cout << ":0" << exec_time_s << std::endl;
	} else {
		std::cout << ":"  << exec_time_s << std::endl;
	}
	
	
	std::cout << "done" << std::endl;
	return 0;
}