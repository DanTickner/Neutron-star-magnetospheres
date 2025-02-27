/*
cd OneDrive\PhD\Codes\20230810 Time evolution with updated equations\Test codes
g++ Test_coordinates.cpp -o Test_coordinates --std=c++11
Test_coordinates

Test the creation of lists of coordinates by the evolution code, to ensure that all endpoints are considered, and the correct number of datapoints is indeed created.
Option to override parameters set in Header_Initial_Conditions_xx.h for smaller numbers of coordinates.
Ouput a CSV file which can easily be checked.
*/

#include <iostream>
#include <iomanip>
#include <vector>
#include <math.h>
#include <fstream>
#include <chrono>

#include "../../Vector_operations.h"
#include "../Header_Initial_Conditions_04.h"
#include "../Header_Vector_B_Rotating_Dipole_Ramp.h"
#include "../Header_Time_Evolution_31.h"



int main(){
	
	
	//----- Extra variables useful only for this testing code -----
	std::string output_filename_test_base = "20240223_Test_Coordinates";
	
	std::string output_filename_test_r   = "../CSV/"  + output_filename_test_base + "_r.csv";
	std::string output_filename_test_t   = "../CSV/"  + output_filename_test_base + "_t.csv";
	std::string output_filename_test_log = "../Logs/" + output_filename_test_base + ".txt";
	
	std::ofstream output_file_test_r;
	std::ofstream output_file_test_t;
	std::ofstream output_file_test_log;
	
	
	double mean_abs_dr_by_di  = 0;
	double abs_d2r_by_di2     = 0;
	double max_abs_d2r_by_di2 = 0;
	
	double mean_abs_dt_by_dj  = 0;
	double abs_d2t_by_dj2     = 0;
	double max_abs_d2t_by_dj2 = 0;
	
	double expected_abs_dr_by_di = ( r_max - r_min ) / ( (double) n_points_r - 1.0 );
	double expected_abs_dt_by_dj = pi / ( (double) n_points_t - 1.0 );
	
	
	
	//----- Calculate the gridpoints using the function as it appears in the evolution code -----
	
	calculate_gridpoints();
	
	std::vector<double> dr_by_di  ( r.size() );
	std::vector<double> d2r_by_di2( r.size() );
	
	std::vector<double> dt_by_dj  ( t.size() );
	std::vector<double> d2t_by_dj2( t.size() );
	
	//----- Output to CSV -----
	
	//--- Radial coordinate ---
	output_file_test_r.open( output_filename_test_r );
	output_file_test_r << "i,r_i,(dr_by_di)_i,(d2r_by_di2)_i\n";
	
	for( int i=0; i<r.size(); i++ ){
		
		if( i >= 1 ){
			dr_by_di[i] = r[i] - r[i-1];
		}
		
		if( i >= 2 ){
			d2r_by_di2[i] = dr_by_di[i] - dr_by_di[i-1];
		}
		
		mean_abs_dr_by_di += std::abs( dr_by_di[i] );
		
		abs_d2r_by_di2 = std::abs( d2r_by_di2[i] );
		if( abs_d2r_by_di2 >= max_abs_d2r_by_di2 ){
			max_abs_d2r_by_di2 = abs_d2r_by_di2;
		}
		
		output_file_test_r << i <<","<< std::setprecision(15) << r[i] <<","<< dr_by_di[i] <<","<< d2r_by_di2[i] <<"\n";
	}
	
	std::cout << "\nCSV file saved:\t" << output_filename_test_r << std::endl;
	
	
	mean_abs_dr_by_di /= (double) r.size() - 1.0;	// The i=0 term is ignored, so there is one less datapoint than radial coordinates.
	
	std::cout << "Mean value of |dr_by_di|  :\t" << std::setprecision(15) << mean_abs_dr_by_di     << std::endl;
	std::cout << "Exp. value of |dr_by_di|  :\t" << std::setprecision(15) << expected_abs_dr_by_di << std::endl;
	std::cout << "Max  value of |d2r_by_di2|:\t" << std::setprecision(15) << max_abs_d2r_by_di2    << std::endl;
	
	output_file_test_log.open( output_filename_test_log );
	output_file_test_log << "Mean value of |dr_by_di|  :\t" << std::setprecision(15) << mean_abs_dr_by_di     << "\n";
	output_file_test_log << "Exp. value of |dr_by_di|  :\t" << std::setprecision(15) << expected_abs_dr_by_di << "\n";
	output_file_test_log << "Max  value of |d2r_by_di2|:\t" << std::setprecision(15) << max_abs_d2r_by_di2    << "\n";
	
	
	//--- Polar coordinate ---
	output_file_test_t.open( output_filename_test_t );
	output_file_test_t << "j,t_j,(dt_by_dj)_j,(d2t_by_dj2)_j\n";
	
	for( int j=0; j<t.size(); j++ ){
		
		if( j >= 1 ){
			dt_by_dj[j] = t[j] - t[j-1];
		}
		
		if( j >= 2 ){
			d2t_by_dj2[j] = dt_by_dj[j] - dt_by_dj[j-1];
		}
		
		mean_abs_dt_by_dj += std::abs( dt_by_dj[j] );
		
		abs_d2t_by_dj2 = std::abs( d2t_by_dj2[j] );
		if( abs_d2t_by_dj2 >= max_abs_d2t_by_dj2 ){
			max_abs_d2t_by_dj2 = abs_d2t_by_dj2;
		}
		
		output_file_test_t << j <<","<< std::setprecision(15) << t[j] <<"\n";
	}
	
	std::cout << "\nCSV file saved:\t" << output_filename_test_t << std::endl;
	
	mean_abs_dt_by_dj /= (double) t.size() - 1.0;	// The i=0 term is ignored, so there is one less datapoint than radial coordinates.
	
	std::cout << "\nMean value of |dt_by_dj|  :\t" << std::setprecision(15) << mean_abs_dt_by_dj     << std::endl;
	std::cout <<   "Exp. value of |dt_by_dj|  :\t" << std::setprecision(15) << expected_abs_dt_by_dj << std::endl;
	std::cout <<   "Max  value of |d2t_by_dj2|:\t" << std::setprecision(15) << max_abs_d2t_by_dj2    << std::endl;
	
	output_file_test_log << "\nMean value of |dt_by_dj|  :\t" << std::setprecision(15) << mean_abs_dt_by_dj     << "\n";
	output_file_test_log <<   "Exp. value of |dt_by_dj|  :\t" << std::setprecision(15) << expected_abs_dt_by_dj << "\n";
	output_file_test_log <<   "Max  value of |d2t_by_dj2|:\t" << std::setprecision(15) << max_abs_d2t_by_dj2    << "\n";
	
	std::cout << "\nLog file saved:\t" << output_filename_test_log << std::endl;
	
	
	
	
	//----- Code finished -----
	std::cout << "\nCode finished" << std::endl;
	
	return 0;
	
}