/*
cd OneDrive\PhD\Codes\20230810 Time evolution with updated equations/Test codes
g++ Test_VSH_series_axisymmetric_zero_divergence.cpp -o Test_VSH_series_axisymmetric_zero_divergence --std=c++11
Test_VSH_series_axisymmetric_zero_divergence

Test the creation of lists of coordinates by the evolution code, to ensure that all endpoints are considered, and the correct number of datapoints is indeed created.
Option to override parameters set in Header_Initial_Conditions_xx.h for smaller numbers of coordinates.
Ouput a CSV file which can easily be checked.

REASON ARCHIVED: Want to use as many functions from the original header file as possible. Replaced E_recalculated by a new dedicated function in the header file, and also use the original stdev for vector values instead of defining it again.
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


//----- Declare functions -----
double E_r_function( double r, double t );
double E_t_function( double r, double t );
double E_p_function( double r, double t );



//----- Main code
int main(){
	
	
	//----- Extra variables useful only for this testing code -----
	//n_points_r = 10;
	//n_points_t = 5;
	
	std::string output_filename_test_base = "../CSV/20240226_Test_VSH_series_axisymmetric_nonzero_divergence";	// Don't add ".csv" as the code will append _r, _t, etc to the filename for each csv.
	
	std::string output_filename_test_coeffs         = output_filename_test_base + "_coeffs.csv";
	std::string output_filename_test_coeffs_stdev   = output_filename_test_base + "_coeffs_stdev.csv";
	std::string output_filename_test_values         = output_filename_test_base + "_values.csv";
	std::string output_filename_test_values_stdev_r = output_filename_test_base + "_values_stdev_r.csv";
	std::string output_filename_test_values_stdev_t = output_filename_test_base + "_values_stdev_t.csv";
	
	std::ofstream output_file_test_coeffs;
	std::ofstream output_file_test_coeffs_stdev;
	std::ofstream output_file_test_values;
	std::ofstream output_file_test_values_stdev_r;
	std::ofstream output_file_test_values_stdev_t;
	
	std::vector< std::vector< std::vector<double> > > E_VSH_coeffs_exact      ( 3, std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( ell_max+1  ) ) );	// Known VSH coefficients, to be input manually.
	std::vector< std::vector< std::vector<double> > > E_VSH_coeffs_difference ( 3, std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( ell_max+1  ) ) );	// Difference between calculated and exact VSH coefficients.
	std::vector< std::vector< std::vector<double> > > E_recalculated          ( 3, std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) ) );	// Original vector recalculated from the VSH coefficients.
	std::vector< std::vector< std::vector<double> > > E_difference            ( 3, std::vector< std::vector<double> > ( n_points_r, std::vector<double> ( n_points_t ) ) );	// Difference between recalculated and original vector.
	
	std::vector< std::vector< double> > stdev_coeffs   ( 3, std::vector<double> ( ell_max+1 ) );
	std::vector< std::vector< double> > stdev_values_r ( 3, std::vector<double> ( n_points_r ) );	// Stdev of recalculated vector as a function of r.
	std::vector< std::vector< double> > stdev_values_t ( 3, std::vector<double> ( n_points_t ) );	// Stdev of recalculated vector as a function of theta.
	
	std::vector<double> stdev_values_tot ( 3 ); 	// Stdev of recalculated vector across all coordinate values.
	
	bool output_recalculated_values_to_csv = false;	// Option to user. Filesizes can become extremely large when using large numbers of gridpoints.
	
	
	//----- Calculate the gridpoints and associated Legendre functions using the function as it appears in the evolution code -----
	calculate_gridpoints();
	calculate_associated_legendre_functions();
	calculate_sqrt_2Lplus1_over_4pi();
	
	
	
	//----- Build array of exact vector values -----
	for( int i=0; i<n_points_r; i++ ){
		for( int j=0; j<n_points_t; j++ ){
			E[0][i][j] = E_r_function( r[i], t[j] );
			E[1][i][j] = E_t_function( r[i], t[j] );
			E[2][i][j] = E_p_function( r[i], t[j] );
		}
	}
	
	
	//----- Build array of known exact VSH series coefficients -----
	// These values must be put in manually, and will change depending on the vector function being tested!
	
	for( int i=0; i<n_points_r; i++ ){
		E_VSH_coeffs_exact[0][i][1] = 4.0 * sqrt( pi / 3.0 ) * pow( r[i], -3 );
		E_VSH_coeffs_exact[1][i][1] = -0.5 * E_VSH_coeffs_exact[0][i][1];
	}
	
	
	//----- Calculate VSH series coefficients -----
	VSH_decomposition_for_E();
	
	
	//----- Output VSH series coefficients to CSV file -----
	output_file_test_coeffs.open( output_filename_test_coeffs );
	output_file_test_coeffs << "k,i,ell,VSH_coeff_calc,VSH_coeff_exact,difference\n";
	
	for( int k=0; k<3; k++ ){
		for( int i=0; i<n_points_r; i++ ){
			for( int ell=0; ell<=ell_max; ell++ ){
				
				E_VSH_coeffs_difference[k][i][ell] = E_VSH_coeffs[k][i][ell] - E_VSH_coeffs_exact[k][i][ell];
				
				output_file_test_coeffs << k <<","<< i <<","<< ell <<","<<std::setprecision(15)<< E_VSH_coeffs[k][i][ell] <<","<<std::setprecision(15)<< E_VSH_coeffs_exact[k][i][ell] <<","<<std::setprecision(15)<< E_VSH_coeffs_difference[k][i][ell] << "\n";
				
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
				stdev_coeffs[k][ell] += pow( E_VSH_coeffs_difference[k][i][ell], 2 );
			}
			
			stdev_coeffs[k][ell] = sqrt( stdev_coeffs[k][ell] / ( (double) n_points_r ) );
			
			output_file_test_coeffs_stdev << k <<","<< ell <<","<< stdev_coeffs[k][ell] <<"\n";
			
		}
	}
	
	std::cout << "CSV file saved:\t" << output_filename_test_coeffs_stdev << std::endl;
	
	
	//----- Recalculate vector values from VSH series -----
	for( int i=0; i<n_points_r; i++ ){
		for( int j=0; j<n_points_t; j++ ){
			
			for( int ell=0; ell<=ell_max; ell++ ){
				E_recalculated[0][i][j] += E_VSH_coeffs[0][i][ell] * sqrt_2Lplus1_over_4pi[ell] * P0[j][ell];
				E_recalculated[1][i][j] += E_VSH_coeffs[1][i][ell] * sqrt_2Lplus1_over_4pi[ell] * P1[j][ell];
				E_recalculated[2][i][j] += E_VSH_coeffs[2][i][ell] * sqrt_2Lplus1_over_4pi[ell] * P1[j][ell];
			}
			
			for( int k=0; k<3; k++ ){
				E_difference[k][i][j] = E[k][i][j] - E_recalculated[k][i][j];
			}
		}
	}
	
	
	//----- Output vector values to CSV file -----
	
	// Failsafe to avoid large file creations. Comment-out if output wanted.
	if( n_points_r * n_points_t > 3e5 ){
		output_recalculated_values_to_csv = false;
	}
	
	if( output_recalculated_values_to_csv ){
		output_file_test_values.open( output_filename_test_values );
		output_file_test_values << "i,j,exact_r,exact_t,exact_p,recalculated_r,recalculated_t,recalculated_p,difference_r,difference_t,difference_p\n";
		
		for( int i=0; i<n_points_r; i++ ){
				for( int j=0; j<n_points_t; j++ ){
					output_file_test_values << i <<","<< j
											<<","<< E             [0][i][j] <<","<< E             [1][i][j] <<","<< E             [2][i][j]
											<<","<< E_recalculated[0][i][j] <<","<< E_recalculated[1][i][j] <<","<< E_recalculated[2][i][j]
											<<","<< E_difference  [0][i][j] <<","<< E_difference  [1][i][j] <<","<< E_difference  [2][i][j] <<"\n";
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
				stdev_values_r[k][i] += pow( E_difference[k][i][j], 2 );
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
				stdev_values_t[k][j] += pow( E_difference[k][i][j], 2 );
			}
			
			stdev_values_t[k][j] = sqrt( stdev_values_t[k][j] / ( (double) n_points_t ) );
			
		}
			
		output_file_test_values_stdev_t << j <<","<< stdev_values_t[0][j] <<","<< stdev_values_t[1][j] <<","<< stdev_values_t[2][j] << "\n";
			
	}
	
	std::cout << "CSV file saved:\t" << output_filename_test_values_stdev_t << std::endl;
	
	
	//----- Stdev for each vector across both coordinates -----
	for( int k=0; k<3; k++ ){
		
		for( int i=0; i<n_points_r; i++ ){
			for( int j=0; j<n_points_t; j++ ){
				stdev_values_tot[k] += pow( E_difference[k][i][j], 2 );
			}
		}
		
		stdev_values_tot[k] = sqrt( stdev_values_tot[k] / ( (double) n_points_r * n_points_t ) );
		
	}
		
	std::cout << "\nstdev of r-component across all coordinates:\t" << stdev_values_tot[0] << std::endl;
	std::cout <<   "stdev of t-component across all coordinates:\t" << stdev_values_tot[1] << std::endl;
	std::cout <<   "stdev of p-component across all coordinates:\t" << stdev_values_tot[2] << std::endl;
			
	
	
	//----- Code finished -----
	std::cout << "\nCode finished" << std::endl;
	
	return 0;
	
}





//----- Define functions -----
double E_r_function( double r, double t ){
	return 2.0 * cos( t ) * pow( r, -3 );
}

double E_t_function( double r, double t ){
	return sin( t ) * pow( r, -3 );
}

double E_p_function( double r, double t ){
	return 0;
}