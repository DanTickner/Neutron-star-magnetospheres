/*
cd OneDrive\PhD\Codes\20240311 Time evolution code\Test codes
g++ Test_time_and_rotation.cpp -o Test_time_and_rotation --std=c++11
Test_time_and_rotation

Test that the values of time and angular velocity agree with those expected. Test that rotation ramps up between the specified points and reaches the intended maximum rate.
*/

#include <iostream>
#include <iomanip>
#include <vector>
#include <math.h>
#include <fstream>
#include <chrono>

#include "../Header_Initial_Conditions_08.h"
#include "../Header_Time_Evolution_38.h"



int main(){
	
	
	//----- Extra variables useful only for this testing code -----
	std::string output_filename_test_base = "20240404_Test_time_and_rotation_new";
	
	std::string output_filename_test_csv = "../CSV/"  + output_filename_test_base + ".csv";
	
	std::ofstream output_file_test_csv;

	
	
	
	//----- Calculate the time, angular velocity and star rotation angle values using the function as it appears in the evolution code -----
	calculate_time_and_rotation_values();
	
	
	//----- Output to CSV -----
	output_file_test_csv.open( output_filename_test_csv );
	output_file_test_csv << "T_index,T,T_SI,Omega,Omega_SI,P,P_SI,R_LC,R_LC_SI,star_rotation_angle\n";
	
	for( int T_index=0; T_index<T.size(); T_index++ ){
		output_file_test_csv <<std::setprecision(precision_csv)<< T_index <<","<< T    [T_index] <<","<< T_SI    [T_index]
		                     <<","<< Omega[T_index] <<","<< Omega_SI[T_index]
							 <<","<< P    [T_index] <<","<< P_SI    [T_index]
							 <<","<< R_LC [T_index] <<","<< R_LC_SI [T_index]
							 <<","<< star_rotation_angle[T_index] << std::endl;
	}
	
	std::cout << "CSV file saved:\t" << output_filename_test_csv << std::endl;
	
	
	
	//----- Code finished -----
	std::cout << "\nCode finished" << std::endl;
	
	return 0;
	
}