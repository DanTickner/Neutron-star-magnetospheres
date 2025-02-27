/*
cd OneDrive\PhD\Codes\20240311 Time evolution code\Test codes
g++ 20241015_Test_twist_smoothing_function.cpp -o 20241015_Test_Twist_Smoothing_Function --std=c++11
20241015_Test_twist_smoothing_function

Output the smoothed twist over the coordinates to a CSV, for plotting with a second code and testing that the smoothing works as expected.

*/

#include <iostream>
#include <iomanip>
#include <vector>
#include <math.h>
#include <fstream>
#include <chrono>

#include "../../Vector_operations.h"
#include "../Header_Initial_Conditions_11.h"
#include "../Header_Time_Evolution_45.h"




//----- Main code -----
int main(){
	
	
	//----- Calculate the gridpoints and build the array of smoothed values -----
	calculate_gridpoints();
	calculate_functions_of_gridpoints();
	setup_twist_final_Omega();
	
	
	//----- Output to CSV -----
	std::string output_filename = "../CSV/20250217_Test_twist_smoothing_function_values.csv";
	
	std::ofstream output_file;
	
	output_file.open( output_filename );
	
	output_file << "r_1," << twist_r_min << ",r_2," << twist_r_max << ",t_1," << twist_t_min << ",t_2," << twist_t_max <<"\n";
	
	output_file << "i,j,r,t,in_twisted_region,twist,twist_over_unsmoothed\n";
	
	for( int i=0; i<n_r; i++ ){
		for( int j=0; j<n_t; j++ ){
			
			bool in_twisted_region = ( r[i] >= twist_r_min ) && ( r[i] <= twist_r_max ) && ( t[j] >= twist_t_min ) && ( t[j] <= twist_t_max );
			
			output_file << i <<","<< j <<","<< r[i] <<","<< t[j] <<","<< in_twisted_region <<","<< twist_final_Omega_values[i][j] <<","<< twist_final_Omega_values[i][j] / twist_final_Omega << "\n";
			
		}
	}
	
	
	
	//----- Code finished -----
	std::cout << "\nCSV file saved:\t" << output_filename << std::endl;
	std::cout << "\nCode finished" << std::endl;
	
	return 0;
	
}