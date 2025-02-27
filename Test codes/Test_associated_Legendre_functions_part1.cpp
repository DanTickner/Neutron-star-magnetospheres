/*
cd OneDrive\PhD\Codes\20240311 Time evolution code\Test codes
g++ Test_associated_Legendre_functions_part1.cpp -o Test_associated_Legendre_functions_part1 --std=c++11
Test_associated_Legendre_functions_part1

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
#include "../Header_Initial_Conditions_07.h"
#include "../Header_Time_Evolution_37.h"



int main(){
	
	
	//----- Extra variables useful only for this testing code -----
	std::string output_filename_test = "../CSV/20240329_Test_Associated_Legendre_Functions.csv";
	
	std::ofstream output_file_test;
	
	
	
	//----- Calculate the gridpoints and associated Legendre functions using the function as it appears in the evolution code -----
	
	calculate_gridpoints();
	
	calculate_associated_legendre_functions();
	
	
	//----- Output to CSV -----
	output_file_test.open( output_filename_test );
	output_file_test << "j,t_j,ell,P_ell_0,P_ell_1\n";
	
	for( int j=0; j<t.size(); j++ ){
		
		for( int ell=0; ell<P0[0].size(); ell++ ){
			
			output_file_test << j <<","<<std::setprecision(15)<< t[j] <<","<< ell <<","<<std::setprecision(15)<< P0[j][ell] <<","<<std::setprecision(15)<< P1[j][ell] <<"\n";
			
		}
		
	}
	
	std::cout << "CSV file saved:\t" << output_filename_test << std::endl;
	
	
	//----- Code finished -----
	std::cout << "\nCode finished" << std::endl;
	
	return 0;
	
}
	