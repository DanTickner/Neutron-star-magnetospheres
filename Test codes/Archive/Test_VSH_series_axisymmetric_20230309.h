/*
cd OneDrive\PhD\Codes\20230810 Time evolution with updated equations\Test codes
g++ Test_VSH_series_axisymmetric.cpp -o Test_VSH_series_axisymmetric --std=c++11
Test_VSH_series_axisymmetric

Test the calculation of the VSH series of an axisymmetric spherical polar vector A.
Three possible expressions to test.
1) The "original" method where we calculate the (1)-coefficients as normal. You must use this option if the vector has nonzero divergence.
2) The divergenceless method where we calculate d/dr A^{r,\ell} by finite differencing.
3) The divergenceless method where we calculate d/dr A^{r,\ell} by a Chebyshev series.
This is chosen by int_test_type at the beginning of the code.

Ouput CSV files with vector values, the VSH series coefficients as a function of r, and standard deviations as a function of each coordinate.
*/

#include <iostream>
#include <iomanip>
#include <vector>
#include <math.h>
#include <fstream>
#include <chrono>

#include "../../Vector_operations.h"
#include "../Header_Initial_Conditions_06.h"
#include "../Header_Time_Evolution_35.h"




//----- Main code -----
int main(){
	
	
	//----- Extra variables useful only for this testing code -----
	
	int test_type = 3;	// 1 = Original; 2 = Divergenceless with finite differencing; 3 = Divergenceless with Chebyshev decomposition; 
	
	std::string output_filename_test_base = "20240309_Test_VSH_series_axisymmetric";
	
	switch( test_type ){
		case 1:
			output_filename_test_base += "_Original";
			std::cout << "Test type:\tOriginal" << std::endl;
			break;
		case 2:
			output_filename_test_base += "_Divergenceless_FD";
			std::cout << "Test type:\tDivergenceless with finite differencing" << std::endl;
			break;
		case 3:
			output_filename_test_base += "_Divergenceless_Chebyshev_nmax_" + std::to_string( n_max );
			std::cout << "Test type:\tDivergenceless with Chebyshev series" << std::endl;
			break;
		default:
			std::cout << "Test type not recognised. Stopping the code.";
			return 1;
	}
	
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
	
	std::vector< std::vector< std::vector<double> > > VSH_coeffs_exact      ( 3, std::vector< std::vector<double> > ( ell_max+1, std::vector<double> ( n_points_r ) ) );	// Known VSH coefficients, to be input manually.
	
	std::vector< std::vector< double> > stdev_coeffs   ( 3, std::vector<double> ( ell_max+1 ) );
	std::vector< std::vector< double> > stdev_values_r ( 3, std::vector<double> ( n_points_r ) );	// Stdev of recalculated vector as a function of r.
	std::vector< std::vector< double> > stdev_values_t ( 3, std::vector<double> ( n_points_t ) );	// Stdev of recalculated vector as a function of theta.
	
	bool output_recalculated_values_to_csv = false;	// Option to user. Filesizes can become extremely large when using large numbers of gridpoints.
	int n_points_total_limit_for_vector_output_to_csv = 3e5;	// Failsafe to avoid accidentally making huge files when the above is left as true.
	
	
	//----- Calculate the gridpoints, associated Legendre functions and VSH series coefficients, recalculate vector and calculate stdev of vector using the functions in the evolution code -----
	calculate_gridpoints();
	calculate_associated_legendre_functions();
	calculate_sqrt_2Lplus1_over_4pi();
	calculate_remapped_chebyshev_polynomials();
	apply_initial_field_values_for_B();
	
	
	//----- Perform the VSH decomposition according to the method chosen -----
	switch( test_type ){
		case 1:
			VSH_decomposition( B, B_VSH_coeffs );
			break;
		case 2:
			VSH_decomposition_divergenceless( B, B_VSH_coeffs );
			calculate_radial_derivatives_of_VSH_coeffs_FD();
			calculate_B_1_L_and_B_1_L_dr();
			break;
		case 3:
			VSH_decomposition_divergenceless( B, B_VSH_coeffs );
			calculate_radial_derivatives_of_VSH_coeffs_Chebyshev();
			calculate_B_1_L_and_B_1_L_dr();
			break;
		default:
			std::cout << "Test type not recognised. Stopping the code.";
			return 1;
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
	
	
	//----- Build array of known exact VSH series coefficients -----
	// These values must be put in manually, and will change depending on the vector function being tested!
	
	for( int i=0; i<n_points_r; i++ ){
		VSH_coeffs_exact[0][1][i] = 4.0 * sqrt( pi / 3.0 ) * pow( r[i], -3 );
		VSH_coeffs_exact[1][1][i] = -0.5 * VSH_coeffs_exact[0][1][i];
	}
	
	//----- Output VSH series coefficients to CSV file -----
	output_file_test_coeffs.open( output_filename_test_coeffs );
	output_file_test_coeffs << "k,i,ell,VSH_coeff_calc,VSH_coeff_exact\n";
	
	for( int k=0; k<3; k++ ){
		for( int i=0; i<n_points_r; i++ ){
			for( int ell=0; ell<=ell_max; ell++ ){
				output_file_test_coeffs << k <<","<< i <<","<< ell <<","<<std::setprecision(15)<< B_VSH_coeffs[k][ell][i] <<","<<std::setprecision(15)<< VSH_coeffs_exact[k][ell][i] << "\n";
			}
		}
	}
	
	std::cout << "CSV file saved:\t" << output_filename_test_coeffs << std::endl;
	
	
	//----- Stdev for each value of ell for the r, (1) and (2) coefficients -----
	output_file_test_coeffs_stdev.open( output_filename_test_coeffs_stdev );
	output_file_test_coeffs_stdev << "k,ell,stdev\n";
	
	for( int k=0; k<3; k++ ){
		for( int ell=0; ell<=ell_max; ell++ ){
			
			stdev_coeffs[k][ell] = stdev_of_array( B_VSH_coeffs[k][ell], VSH_coeffs_exact[k][ell] );
			
			output_file_test_coeffs_stdev << k <<","<< ell <<","<< stdev_coeffs[k][ell] <<"\n";
			
		}
	}
	
	std::cout << "CSV file saved:\t" << output_filename_test_coeffs_stdev << std::endl;
	
	
	
	
	//----- Output vector values to CSV file -----
	
	// Failsafe to avoid large file creations. Comment-out if output wanted.
	if( n_points_total > n_points_total_limit_for_vector_output_to_csv ){
		output_recalculated_values_to_csv = false;
	}
	
	if( output_recalculated_values_to_csv ){
		output_file_test_values.open( output_filename_test_values );
		output_file_test_values << "i,j,exact_r,exact_t,exact_p,recalculated_r,recalculated_t,recalculated_p\n";
		
		for( int i=0; i<n_points_r; i++ ){
			for( int j=0; j<n_points_t; j++ ){
				output_file_test_values << i <<","<< j
										<<","<< B                [0][i][j] <<","<< B                [1][i][j] <<","<< B                [2][i][j]
										<<","<< B_from_VSH_coeffs[0][i][j] <<","<< B_from_VSH_coeffs[1][i][j] <<","<< B_from_VSH_coeffs[2][i][j] <<"\n";
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
				stdev_values_r[k][i] += pow( B[k][i][j] - B_from_VSH_coeffs[k][i][j], 2 );
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
				stdev_values_t[k][j] += pow( B[k][i][j] - B_from_VSH_coeffs[k][i][j], 2 );
			}
			
			stdev_values_t[k][j] = sqrt( stdev_values_t[k][j] / ( (double) n_points_r ) );
			
		}
			
		output_file_test_values_stdev_t << j <<","<< stdev_values_t[0][j] <<","<< stdev_values_t[1][j] <<","<< stdev_values_t[2][j] << "\n";
			
	}
	
	std::cout << "CSV file saved:\t" << output_filename_test_values_stdev_t << std::endl;
	
	
	
	
	//----- Code finished -----
	std::cout << "\nCode finished" << std::endl;
	
	return 0;
	
}