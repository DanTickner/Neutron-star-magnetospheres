/*
cd OneDrive\PhD\Codes\20230810 Time evolution with updated equations/Test codes
g++ Test_VSH_series_axisymmetric_zero_divergence_FD.cpp -o Test_VSH_series_axisymmetric_zero_divergence_FD --std=c++11
Test_VSH_series_axisymmetric_zero_divergence_FD

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
#include "../Header_B_static_dipole_A_equals_minus2Psi.h"
#include "../Header_Time_Evolution_31.h"




//----- Main code -----
int main(){
	
	
	//----- Extra variables useful only for this testing code -----
	//n_points_r = 10;
	//n_points_t = 5;
	
	//std::string output_filename_test_base = "20240227_Test_VSH_series_axisymmetric_zero_divergence_FD";
	std::string output_filename_test_base = "20240227_Test_VSH_series_axisymmetric_zero_divergence_Chebyshev_fakedata";
	
	std::string output_filename_test_coeffs         = "../CSV/"  + output_filename_test_base + "_coeffs.csv";
	std::string output_filename_test_coeffs_stdev   = "../CSV/"  + output_filename_test_base + "_coeffs_stdev.csv";
	std::string output_filename_test_values         = "../CSV/"  + output_filename_test_base + "_values.csv";
	std::string output_filename_test_values_stdev_r = "../CSV/"  + output_filename_test_base + "_values_stdev_r.csv";
	std::string output_filename_test_values_stdev_t = "../CSV/"  + output_filename_test_base + "_values_stdev_t.csv";
	std::string output_filename_test_log            = "../Logs/" + output_filename_test_base + ".txt";
	
	std::ofstream output_file_test_coeffs;
	std::ofstream output_file_test_coeffs_stdev;
	std::ofstream output_file_test_values;
	std::ofstream output_file_test_values_stdev_r;
	std::ofstream output_file_test_values_stdev_t;
	std::ofstream output_file_test_log;
	
	std::vector< std::vector< std::vector<double> > > VSH_coeffs_exact      ( 3, std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( ell_max+1  ) ) );	// Known VSH coefficients, to be input manually.
	std::vector< std::vector< std::vector<double> > > VSH_coeffs_difference ( 3, std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( ell_max+1  ) ) );	// Difference between calculated and exact VSH coefficients.
	std::vector< std::vector< std::vector<double> > > difference            ( 3, std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) ) );	// Difference between recalculated and original vector.
	
	std::vector< std::vector< double> > stdev_coeffs   ( 3, std::vector<double> ( ell_max+1 ) );
	std::vector< std::vector< double> > stdev_values_r ( 3, std::vector<double> ( n_points_r ) );	// Stdev of recalculated vector as a function of r.
	std::vector< std::vector< double> > stdev_values_t ( 3, std::vector<double> ( n_points_t ) );	// Stdev of recalculated vector as a function of theta.
	
	bool output_recalculated_values_to_csv = false;	// Option to user. Filesizes can become extremely large when using large numbers of gridpoints.
	
	
	//----- Calculate the gridpoints, associated Legendre functions and VSH series coefficients, recalculate vector and calculate stdev of vector using the functions in the evolution code -----
	calculate_gridpoints();
	calculate_associated_legendre_functions();
	calculate_sqrt_2Lplus1_over_4pi();
	apply_initial_field_values_for_B();
	VSH_decomposition_for_B();
	
	
	
	
	//----- Calculate B_1_L -----
	for( int ell=1; ell<=ell_max; ell++ ){
		
		//--- Calculate first d/dr B^{r,ell} ---
		// Right derivatives at inner boundary.
		B_VSH_coeffs_dr [0][0][ell] = ( B_VSH_coeffs[0][1][ell] - B_VSH_coeffs[0][0][ell] ) / delta_r;
		
		// Symmetric derivatives for all intermediate points.
		for( int i=1; i<n_points_r-1; i++ ){
			B_VSH_coeffs_dr [0][i][ell] = ( B_VSH_coeffs[0][i+1][ell] - B_VSH_coeffs[0][i-1][ell] ) / ( 2.0*delta_r );
		}
		
		// Left derivatives at outer boundary.
		B_VSH_coeffs_dr[0] .back()[ell] = ( B_VSH_coeffs[0].back()[ell] - B_VSH_coeffs[0][n_points_r-2][ell] ) / delta_r;
		
		//--- Calculate B^{(1),ell} ---
		for( int i=0; i<n_points_r; i++ ){
			B_VSH_coeffs   [1][i][ell] = ( r[i] * B_VSH_coeffs_dr [0][i][ell] + 2.0 * B_VSH_coeffs   [0][i][ell] ) / ( (double) ell*(ell+1.0) );
			B_VSH_coeffs[1][i][ell] *= 0.875;	// Fake data until I can get the Chebyshev method to work.
		}
		
	}
	
	
	
	
	//----- Recalculate vector and calculate stdev of vector using the functions in the evolution code -----
	recalculate_B_from_VSH_coeffs();
	calculate_stdev_of_VSH_decomposition();
	
	std::cout << "\nstdev of r-component across all coordinates:\t" << B_stdev[0] << std::endl;
	std::cout <<   "stdev of t-component across all coordinates:\t" << B_stdev[1] << std::endl;
	std::cout <<   "stdev of p-component across all coordinates:\t" << B_stdev[2] << std::endl;
	
	output_file_test_log.open( output_filename_test_log );
	output_file_test_coeffs << "k,i,ell,VSH_coeff_calc,VSH_coeff_exact,difference\n";
	output_file_test_log << "stdev of r-component across all coordinates:\t" << B_stdev[0] << "\n";
	output_file_test_log << "stdev of t-component across all coordinates:\t" << B_stdev[1] << "\n";
	output_file_test_log << "stdev of p-component across all coordinates:\t" << B_stdev[2] << "\n";
	
	std::cout << "Log file saved:\t" << output_filename_test_log << std::endl;
	
	
	//----- Calculate difference between recalculated and exact vector -----
	for( int i=0; i<n_points_r; i++ ){
		for( int j=0; j<n_points_t; j++ ){
			for( int k=0; k<3; k++ ){
				difference[k][i][j] = B[k][i][j] - B_from_VSH_coeffs[k][i][j];
			}
		}
	}
	
	
	//----- Build array of known exact VSH series coefficients -----
	// These values must be put in manually, and will change depending on the vector function being tested!
	
	for( int i=0; i<n_points_r; i++ ){
		VSH_coeffs_exact[0][i][1] = 4.0 * sqrt( pi / 3.0 ) * pow( r[i], -3 );
		VSH_coeffs_exact[1][i][1] = -0.5 * VSH_coeffs_exact[0][i][1];
	}
	
	
	//----- Output VSH series coefficients to CSV file -----
	output_file_test_coeffs.open( output_filename_test_coeffs );
	output_file_test_coeffs << "k,i,ell,VSH_coeff_calc,VSH_coeff_exact,difference\n";
	
	for( int k=0; k<3; k++ ){
		for( int i=0; i<n_points_r; i++ ){
			for( int ell=0; ell<=ell_max; ell++ ){
				
				VSH_coeffs_difference[k][i][ell] = B_VSH_coeffs[k][i][ell] - VSH_coeffs_exact[k][i][ell];
				
				output_file_test_coeffs << k <<","<< i <<","<< ell <<","<<std::setprecision(15)<< B_VSH_coeffs[k][i][ell] <<","<<std::setprecision(15)<< VSH_coeffs_exact[k][i][ell] <<","<<std::setprecision(15)<< VSH_coeffs_difference[k][i][ell] << "\n";
				
			}
		}
	}
	
	std::cout << "CSV file saved:\t" << output_filename_test_coeffs << std::endl;
	
	
	//----- Stdev for each value of ell for the r, (1) and (2) coefficients -----
	output_file_test_coeffs_stdev.open( output_filename_test_coeffs_stdev );
	output_file_test_coeffs_stdev << "k,ell,stdev\n";
	
	for( int k=0; k<3; k++ ){
		for( int ell=0; ell<=ell_max; ell++ ){
			
			for( int i=0; i<n_points_r; i++ ){
				stdev_coeffs[k][ell] += pow( VSH_coeffs_difference[k][i][ell], 2 );
			}
			
			stdev_coeffs[k][ell] = sqrt( stdev_coeffs[k][ell] / ( (double) n_points_r ) );
			
			output_file_test_coeffs_stdev << k <<","<< ell <<","<< stdev_coeffs[k][ell] <<"\n";
			
		}
	}
	
	std::cout << "CSV file saved:\t" << output_filename_test_coeffs_stdev << std::endl;
	
	
	
	
	//----- Output vector values to CSV file -----
	
	// Failsafe to avoid large file creations. Comment-out if output wanted.
	if( n_points_total > 3e5 ){
		output_recalculated_values_to_csv = false;
	}
	
	if( output_recalculated_values_to_csv ){
		output_file_test_values.open( output_filename_test_values );
		output_file_test_values << "i,j,exact_r,exact_t,exact_p,recalculated_r,recalculated_t,recalculated_p,difference_r,difference_t,difference_p\n";
		
		for( int i=0; i<n_points_r; i++ ){
			for( int j=0; j<n_points_t; j++ ){
				output_file_test_values << i <<","<< j
										<<","<< B                [0][i][j] <<","<< B                [1][i][j] <<","<< B                [2][i][j]
										<<","<< B_from_VSH_coeffs[0][i][j] <<","<< B_from_VSH_coeffs[1][i][j] <<","<< B_from_VSH_coeffs[2][i][j]
										<<","<< difference       [0][i][j] <<","<< difference       [1][i][j] <<","<< difference       [2][i][j] <<"\n";
			}
		}
		
		std::cout << "CSV file saved:\t" << output_filename_test_values << std::endl;
	} else{
		std::cout << "***** "<< output_filename_test_values << " NOT CREATED *****" << std::endl;
	}
	
	
	//----- Stdev for each vector component as a function of r -----
	output_file_test_values_stdev_r.open( output_filename_test_values_stdev_r );
	output_file_test_values_stdev_r << "i,stdev_r,stdev_t,stdev_p\n";
	
	for( int i=0; i<n_points_r; i++ ){
		
		for( int k=0; k<3; k++ ){
			
			for( int j=0; j<n_points_t; j++ ){
				stdev_values_r[k][i] += pow( difference[k][i][j], 2 );
			}
			
			stdev_values_r[k][i] = sqrt( stdev_values_r[k][i] / ( (double) n_points_t ) );
			
		}
			
		output_file_test_values_stdev_r << i <<","<< stdev_values_r[0][i] <<","<< stdev_values_r[1][i] <<","<< stdev_values_r[2][i] << "\n";
			
	}
	
	std::cout << "CSV file saved:\t" << output_filename_test_values_stdev_r << std::endl;
	
	
	//----- Stdev for each vector component as a function of theta -----
	output_file_test_values_stdev_t.open( output_filename_test_values_stdev_t );
	output_file_test_values_stdev_t << "j,stdev_r,stdev_t,stdev_p\n";
	
	for( int j=0; j<n_points_t; j++ ){
		
		for( int k=0; k<3; k++ ){
			
			for( int i=0; i<n_points_r; i++ ){
				stdev_values_t[k][j] += pow( difference[k][i][j], 2 );
			}
			
			stdev_values_t[k][j] = sqrt( stdev_values_t[k][j] / ( (double) n_points_t ) );
			
		}
			
		output_file_test_values_stdev_t << j <<","<< stdev_values_t[0][j] <<","<< stdev_values_t[1][j] <<","<< stdev_values_t[2][j] << "\n";
			
	}
	
	std::cout << "CSV file saved:\t" << output_filename_test_values_stdev_t << std::endl;
	
	
	
	
	//----- Code finished -----
	std::cout << "\nCode finished" << std::endl;
	
	return 0;
	
}