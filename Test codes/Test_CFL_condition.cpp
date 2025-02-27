/*
cd OneDrive\PhD\Codes\20240311 Time evolution code\Test codes
g++ Test_CFL_condition.cpp -o Test_CFL_condition --std=c++11
Test_CFL_condition

Calculate the separation between all adjacent gridpoints.
Determine the minimum separation between them, which will be the maximum timestep by the CFL condition.
Verify the expressions in Cartesian coordinates, spherical polar coordinates and spherical polar coordinates with constant grid spacing by changing int test_type.
*/

#include <iostream>
#include <iomanip>
#include <vector>
#include <math.h>
#include <fstream>
#include <chrono>

#include "../Header_Initial_Conditions_07.h"
#include "../Header_Time_Evolution_36.h"


//----- Declare functions -----
double gridpoint_separation_spherical( double r1, double t1, double r2, double t2 );
double gridpoint_separation_cartesian( double x1, double z1, double x2, double z2 );


//----- Main code -----
int main(){
	
	
	//----- Extra variables useful only for this testing code -----
	
	std::string output_filename_test_base = "20240313_Test_CFL_condition";
	
	int test_type = 4;	// 1 = Cartesian, 2 = spherical polar, 3 = spherical polar for constant grid spacing, 4 = spherical polar constant approximation
	
	switch( test_type ){
		case 1:
			std::cout << "Testing separation expression for Cartesian coordinates." << std::endl;
			output_filename_test_base += "_Cartesian";
			break;
		case 2:
			std::cout << "Testing separation expression for spherical polar coordinates." << std::endl;
			output_filename_test_base += "_Spherical";
			break;
		case 3:
			std::cout << "Testing separation expression for spherical polar coordinates with constant grid spacing." << std::endl;
			output_filename_test_base += "_expressions_for_constant_spacing";
			break;
		case 4:
			std::cout << "Testing separation expression for spherical polar coordinates with constant grid spacing (approximation)." << std::endl;
			output_filename_test_base += "_expressions_for_constant_spacing_approximation";
			break;
		default:
			std::cout << "test_type = " << test_type << " not recognised. Please enter 1 (Cartesian), 2 (spherical polar), 3 (spherical polar for constant grid spacing) or 4 (spherical polar for constant grid spacing approximation). Stopping the code." << std::endl;
			return 1;
	}
	
	std::string output_filename_test_csv = "../CSV/"  + output_filename_test_base + ".csv";
	std::string output_filename_test_log = "../Logs/" + output_filename_test_base + ".txt";
	
	std::ofstream output_file_test_csv;
	std::ofstream output_file_test_log;
	
	output_file_test_log.open( output_filename_test_log );
	switch( test_type ){
		case 1:
			output_file_test_log << "Testing separation expression for Cartesian coordinates." << "\n";
			break;
		case 2:
			output_file_test_log << "Testing separation expression for spherical polar coordinates." << "\n";
			break;
		case 3:
			output_file_test_log << "Testing separation expression for spherical polar coordinates with constant grid spacing." << "\n";
			break;
		case 4:
			output_file_test_log << "Testing separation expression for spherical polar coordinates with constant grid spacing (approximation)." << "\n";
			break;
		default:
			return 1;
	}
	
	std::vector<double> gridpoint_with_smallest_separation_for_constant_r { 999, 999, 999 };
	std::vector<double> gridpoint_with_smallest_separation_for_constant_t { 999, 999, 999 };
	
	double separation_for_constant_r = 0;
	double separation_for_constant_t = 0;
	
	
	
	
	//----- Calculate the gridpoints using the function as it appears in the evolution code -----
	calculate_gridpoints();
	
	
	//----- Output gridpoint spacing for reference -----
	std::cout << "\ndelta_r:\t" <<std::setprecision(csv_precision)<< delta_r << std::endl;
	std::cout <<   "delta_t:\t" <<std::setprecision(csv_precision)<< delta_t << std::endl;
	
	output_file_test_log << "\ndelta_r:\t" <<std::setprecision(csv_precision)<< delta_r << "\n";
	output_file_test_log <<   "delta_t:\t" <<std::setprecision(csv_precision)<< delta_t << "\n";
	
	
	//----- Output to CSV -----
	output_file_test_csv.open( output_filename_test_csv );
	output_file_test_csv << "i,j,r,t,x,z";
	output_file_test_csv << ",separation_for_constant_r";
	output_file_test_csv << "\n";
	
	for( int i=0; i<n_points_r; i++ ){
		for( int j=0; j<n_points_t; j++ ){
			
			//--- Separation with constant r requires us to not be on the final theta-gridpoint ---
			if( j == n_points_t - 1 ){
				separation_for_constant_r = 0;
			}
			else{
				
				switch( test_type ){
					case 1:
						separation_for_constant_r = gridpoint_separation_cartesian( x[i][j], z[i][j], x[i][j+1], z[i][j+1] );
						break;
					case 2:
						separation_for_constant_r = gridpoint_separation_spherical( r[i], t[j], r[i], t[j+1] );
						break;
					case 3:
						separation_for_constant_r = r[i] * sqrt(2.0) * sqrt( 1.0 - cos( delta_t ) );
						break;
					case 4:
						separation_for_constant_r = r[i] * delta_t;
						break;
					default:
						return 1;
				}
				
				
				if( separation_for_constant_r <= gridpoint_with_smallest_separation_for_constant_r[0] ){
					gridpoint_with_smallest_separation_for_constant_r[0] = separation_for_constant_r;
					gridpoint_with_smallest_separation_for_constant_r[1] = r[i];
					gridpoint_with_smallest_separation_for_constant_r[2] = t[j];
				}
				
			}
			
			
			//--- Separation with constant theta requires us to not be on the final r-gridpoint ---
			if( i == n_points_r - 1 ){
				separation_for_constant_t = 0;
			}
			else{
				
				switch( test_type ){
					case 1:
						separation_for_constant_t = gridpoint_separation_cartesian( x[i][j], z[i][j], x[i+1][j], z[i+1][j] );
						break;
					case 2:
						separation_for_constant_t = gridpoint_separation_spherical( r[i], t[j], r[i+1], t[j] );
						break;
					case 3:
						separation_for_constant_t = delta_r;
						break;
					case 4:
						separation_for_constant_t = delta_r;
						break;
					default:
						return 1;
				}
				
				if( separation_for_constant_t <= gridpoint_with_smallest_separation_for_constant_t[0] ){
					gridpoint_with_smallest_separation_for_constant_t[0] = separation_for_constant_t;
					gridpoint_with_smallest_separation_for_constant_t[1] = r[i];
					gridpoint_with_smallest_separation_for_constant_t[2] = t[j];
				}
				
			}
			
			output_file_test_csv << i <<","<< j <<","<<std::setprecision(csv_precision)<< r[i] <<","<< t[j] <<","<< x[i][j] <<","<< z[i][j] <<","<< separation_for_constant_r <<","<< separation_for_constant_t << "\n";
			
		}
	}
	
	std::cout << "\nCSV file saved:\t" << output_filename_test_csv << std::endl;
	
	std::cout << "\nSmallest separation with constant r:\t" <<std::setprecision(csv_precision)<< gridpoint_with_smallest_separation_for_constant_r[0] << "\tat (r,theta) = ( " << gridpoint_with_smallest_separation_for_constant_r[1] <<" , "<< gridpoint_with_smallest_separation_for_constant_r[2] << " )." << std::endl;
	std::cout <<   "Smallest separation with constant t:\t" <<std::setprecision(csv_precision)<< gridpoint_with_smallest_separation_for_constant_t[0] << "\tat (r,theta) = ( " << gridpoint_with_smallest_separation_for_constant_t[1] <<" , "<< gridpoint_with_smallest_separation_for_constant_t[2] << " )." << std::endl;
	
	output_file_test_log << "\nSmallest separation with constant r:\t" <<std::setprecision(csv_precision)<< gridpoint_with_smallest_separation_for_constant_r[0] << "\tat (r,theta) = ( " << gridpoint_with_smallest_separation_for_constant_r[1] <<" , "<< gridpoint_with_smallest_separation_for_constant_r[2] << " )." << "\n";
	output_file_test_log <<   "Smallest separation with constant t:\t" <<std::setprecision(csv_precision)<< gridpoint_with_smallest_separation_for_constant_t[0] << "\tat (r,theta) = ( " << gridpoint_with_smallest_separation_for_constant_t[1] <<" , "<< gridpoint_with_smallest_separation_for_constant_t[2] << " )." << "\n";
	
	std::cout << "\nLog file saved:\t" << output_filename_test_log << std::endl;
	
	
	
	
	//----- Code finished -----
	std::cout << "\nCode finished" << std::endl;
	
	return 0;
	
}




//----- Define functions -----
double gridpoint_separation_spherical( double r1, double t1, double r2, double t2 ){
	return sqrt( r2*r2 + r1*r1 - 2.0*r2*r1*( sin(t2)*sin(t1) + cos(t2)*cos(t1) ) );
}

double gridpoint_separation_cartesian( double x1, double z1, double x2, double z2 ){
	return sqrt( pow( x2-x1, 2 ) + pow( z2-z1, 2 ) );
}