/*
cd OneDrive\PhD\Codes\20240311 Time evolution code\Test codes
g++ Test_Electric_field_rampup.cpp -o Test_Electric_field_rampup --std=c++11
Test_Electric_field_rampup

Test whether the final result from ramping up the electric field matches the expected electric field if rotation was switched on instantaneously.

*/

#include <iostream>
#include <iomanip>
#include <vector>
#include <math.h>
#include <fstream>
#include <chrono>

#include "../../Vector_operations.h"
#include "../Header_Initial_Conditions_08.h"
#include "../Header_Time_Evolution_38.h"



int main(){
	
	
	//----- Extra variables useful only for this testing code -----
	std::string output_filename_test_base = "20240403_Test_Eletric_field_rampup";
	
	std::string output_filename_test_csv_one_point  = "../CSV/"  + output_filename_test_base + "_one_point.csv";
	std::string output_filename_test_csv_all_points = "../CSV/"  + output_filename_test_base + "_all_points.csv";
	std::string output_filename_test_log            = "../Logs/" + output_filename_test_base + ".txt";
	
	std::ofstream output_file_test_csv_one_point;
	std::ofstream output_file_test_csv_all_points;
	std::ofstream output_file_test_log;
	
	std::vector< std::vector< std::vector<double> > > E_final_expected ( 3, std::vector< std::vector<double> >( n_r, std::vector<double> ( n_t ) ) );
	
	
	//----- Calculate the gridpoints using the function as it appears in the evolution code -----
	calculate_time_and_rotation_values();	// Needed in case bool set_r_max_to_2_R_LC = true.
	calculate_gridpoints();
	
	
	//----- Ramp up the electric field -----
	output_file_test_csv_one_point.open( output_filename_test_csv_one_point );
	output_file_test_csv_one_point << "T_index,E_r,E_t,E_p\n";
	
	while( T_index < n_T ){
		
		ramp_up_electric_fields_for_rotation();
		
		output_file_test_csv_one_point <<std::setprecision(precision_csv)<< T_index <<","<< E[0][cout_i][cout_j] <<","<< E[1][cout_i][cout_j] <<","<< E[2][cout_i][cout_j] <<","<< pow(r[cout_i],-2) <<","<< r[cout_i]*pow(r[cout_i],-3) << "\n";
		
		T_index ++;
		
	}
	
	
	//----- Compare final values to expected values -----
	output_file_test_csv_all_points.open( output_filename_test_csv_all_points );
	output_file_test_csv_all_points << "T_index,E_r,E_t,E_p,E_r_expected,E_t_expected,E_p_expected\n";
	
	for( int i=0; i<n_r; i++ ){
		for( int j=0; j<n_t; j++ ){
			
			double factor = Omega.back() * r[i] * st[j];
			
			E_final_expected[0][i][j] = factor *  B_t_function( r[i], t[j] );
			E_final_expected[1][i][j] = factor * -B_r_function( r[i], t[j] );
			E_final_expected[2][i][j] = 0;
			
			output_file_test_csv_all_points << E[0][i][j] <<","<< E[1][i][j] <<","<< E[2][i][j] <<","<<
			                                   E_final_expected[0][i][j] <<","<< E_final_expected[1][i][j] <<","<< E_final_expected[2][i][j];
			
		}
	}
	
	double stdev_E_r_final = stdev_of_array_2d( E[0], E_final_expected[0] );
	double stdev_E_t_final = stdev_of_array_2d( E[1], E_final_expected[1] );
	double stdev_E_p_final = stdev_of_array_2d( E[2], E_final_expected[2] );
	
	output_file_test_log.open( output_filename_test_log );
	
	std::cout << "Stdev of final electric field in radial direction:\t" << stdev_E_r_final << std::endl;
	std::cout << "Stdev of final electric field in theta  direction:\t" << stdev_E_t_final << std::endl;
	std::cout << "Stdev of final electric field in phi    direction:\t" << stdev_E_p_final << std::endl;
	
	output_file_test_log << "Stdev of final electric field in radial direction:\t" << stdev_E_r_final << "\n";
	output_file_test_log << "Stdev of final electric field in theta  direction:\t" << stdev_E_t_final << "\n";
	output_file_test_log << "Stdev of final electric field in phi    direction:\t" << stdev_E_p_final << "\n";
	
	
	//----- Code finished -----
	std::cout << "\nCSV file saved:\t" << output_filename_test_csv_one_point  << std::endl;
	std::cout <<   "CSV file saved:\t" << output_filename_test_csv_all_points << std::endl;
	std::cout <<   "Log file saved:\t" << output_filename_test_log            << std::endl;
	
	std::cout << "\nCode finished" << std::endl;
	
	return 0;
	
}