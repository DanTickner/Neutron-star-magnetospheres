/*
cd OneDrive\PhD\Codes\20230810 Time evolution with updated equations
g++ Time_Evolution_01.cpp -o Time_Evolution_01
Time_Evolution_01
Evolve the magnetic and electric fields over time, by finite differencing the radial direction and using a VSH expansion in the theta direction.
	 

Based on 20230707 Time Evolution with radial finite differencing/Time_Evolution_05.cpp
*/

#include <iostream>
#include <iomanip>
#include <vector>
#include <math.h>
#include <fstream>
#include <chrono>

#include "Headers/Initial_conditions_01.h"
#include "../Vector_operations.h"
#include "Headers/Vector_B_dipole.h"
#include "Headers/Vector_E_rotating_dipole.h"


//----- Define variables -----
int n_points_r  = 200;
int n_points_t  = 100;
int n_timesteps = 300;

double delta_T = 0.01;

int w = 14;			// fixed width for text output


std::string output_filename = "CSV/20230816_apply_FF_condition_delete.csv";

int csv_write_freq_r = n_points_r/50;	// Write to CSV after this many steps in the radial direction.
int csv_write_freq_t = n_points_t/50;	// Write to CSV after this many steps in the polar  direction.
int csv_write_freq_T = 1;				// Write to CSV after this many timesteps.

int cout_i = 4;
int cout_j = 7;
int cout_freq_T = 1;

#include "Headers/Header_Time_Evolution_01.h"



int main(){
	
	//----- Fail if cout gridpoints exceed number of gridpoints, to prevent code failure due to list overflow with no accompanying message -----
	if( cout_i > n_points_r ){
		std::cout << "Index for screen output in radial direction exceeds number of gridpoints (" << cout_i << " vs " << n_points_r <<"). Stopping the code." << std::endl;
		return 0;
	}
	if( cout_j > n_points_t ){
		std::cout << "Index for screen output in colatitude direction exceeds number of gridpoints (" << cout_j << " vs " << n_points_t <<"). Stopping the code." << std::endl;
		return 0;
	}
	
	
	//----- Allow user to determine whether CFL condition verified and inform of total length of simulation -----
	double csv_filesize_factor = 3e-4;	// Adjust to give more accurate esimates.
	double csv_filesize_guess  = csv_filesize_factor * n_timesteps * n_points_r * n_points_t / ( (double) csv_write_freq_T * csv_write_freq_r * csv_write_freq_t );
	
	std::cout << "Parameter values (given in code units, then in SI units)" << std::endl;
	std::cout << "Omega           :\t" <<std::left<<std::setw(16)<< Omega << std::endl;
	std::cout << "R_LC            :\t" <<std::left<<std::setw(16)<< R_LC << std::endl;
	std::cout << "r_min, r_max    :\t" <<std::left<<std::setw(16)<< r_min <<std::left<<std::setw(16)<< r_max << std::endl;
	std::cout << "Timestep        :\t" <<std::left<<std::setw(16)<< delta_T 	<<std::left<<std::setw(16)<< delta_T * T_factor << std::endl;
	std::cout << "CFL max timestep:\t" <<std::left<<std::setw(16)<< delta_r <<std::left<<std::setw(16)<< delta_r * T_factor << std::endl;
	std::cout << "Length of simulation:\t" <<std::left<<std::setw(16)<< n_timesteps * delta_T <<  n_timesteps * delta_T * T_factor <<"\tsteps: " << n_timesteps << std::endl;
	std::cout << "Est. CSV size (MB):\t" <<  csv_filesize_guess << std::endl;
	
	
	
	
	//----- Build lists of gridpoints and precalculate area element for integral -----
	calculate_gridpoints();
	
	std::cout << n_max <<"\t"<< ell_max <<"\t" <<cout_i <<"\t"<< cout_j <<"\t"<< r[cout_i] <<"\t"<< t[cout_j] <<"\t"<< Omega <<"\n"<< std::endl;
	
	
	
	
	//----- Precalculate basis functions and AVSH series coeffs of cross products of VSHs -----
	calculate_associated_legendre_functions();
	calculate_AVSH_series_coeffs_of_VSH_cross_products();
	
	
	
	
	//----- Prepare output CSV and print headers to screen -----
	output_file.open( output_filename );
	output_headers_to_csv();
	output_headers_to_screen();
	
	//----- Build initial arrays of VSH coefficient values. Start timing now. Output initial values. -----
	time_start_seconds = std::chrono::high_resolution_clock::now().time_since_epoch().count() * 1e-9;
	
	apply_initial_VSH_coeff_values();
	calculate_radial_derivatives_of_VSH_coeffs_inner_boundary();
	calculate_radial_derivatives_of_VSH_coeffs();
	calculate_B_1_L_and_B_1_L_dr_inner_boundary( false );
	calculate_B_1_L_and_B_1_L_dr( false );
	calculate_alpha_and_beta();
	calculate_Cartesian_vector_components();
	output_to_screen( 0 );
	output_to_csv( 0 );
	
	//----- Evolve to future timesteps -----
	
	
	for( int T_index=1; T_index<=n_timesteps; T_index++ ){
	
		T += delta_T;
		
		reset_arrays();
		
		for( int i=1; i<n_points_r; i++ ){
			// Not evolving the i=0 coefficients is equivalent to applying Petri's inner boundary conditions (Petri 2012 S3.3 and Book 20, pp.16-17).
		
			//--- Precalculate integrals of alpha and beta times associated Legendre polmynomials ---
			calculate_alpha_and_beta_integrals( i );
			
			//--- Calculate time derivatives ---
			calculate_ExB_VSH_coeffs( i );
			calculate_time_derivatives_of_VSH_coeffs( i );
			integrate_VSH_coeffs_wrt_time( i );
			
			//--- Apply force-free conditions ---	// This makes the code take ~16 times as long, due to many more integrals needing to be performed.
			/*
			calculate_B_tilde_integrals( i );
			calculate_E_tilde_coeffs( i );
			calculate_E_corrected_integrals( i );
			calculate_E_corrected_coeffs( i );
			*/
		}
		
		update_E_within_light_cylinder_if_greater_than_Bsq_and_redo_VSH_decomposition();
		
		
		
		//--- Calculate 2D arrays of alpha(r,theta) and beta(r,theta) ---
		calculate_alpha_and_beta();
		
		//--- Calculate B^{(1),L} and output ---
		calculate_radial_derivatives_of_VSH_coeffs();
		calculate_B_1_L_and_B_1_L_dr();
		calculate_Cartesian_vector_components();
		output_to_csv( T_index );
		
		output_to_screen( T_index );
	}
	
	
	
	
	
	/* THIS WAS ALREADY A MULTILINE COMMENT 20230810
	
	//----- Evolve to future timesteps -----
	
	for( int T_index=0; T_index<=n_timesteps; T_index++ ){
		
		calculate_CVSH_coeffs();
		
		
		//--- Calculate spatial derivatives and perform Euler integration ---
		
		for( int i=0; i<n_points_r; i++ ){
			for( int j=0; j<n_points_t; j++ ){
				
				calculate_spatial_derivatives( i, j );
				calculate_time_derivatives( i, j );
				output_to_csv( i, j, T_index );
								
				//--- Perform Euler integration and overwrite vector values ---
				// Recall that integration gives the value at the NEXT time step, so this should be done after CSV output.
				
				B_r[i][j] += B_r_dT * delta_T;
				B_t[i][j] += B_t_dT * delta_T;
				B_p[i][j] += B_p_dT * delta_T;
				E_r[i][j] += E_r_dT * delta_T;
				E_t[i][j] += E_t_dT * delta_T;
				E_p[i][j] += E_p_dT * delta_T;
				
			}
		}
		
		
		//--- Enforce boundary conditions and force-free conditions ---
		enforce_boundary_condition_outer();
		enforce_force_free_condition_E_less_than_B();
		
		
		//--- Output to screen ---
		output_to_screen( T_index );
		
		T += delta_T;
	}
	*/
	
	//----- Output execution time and finish code -----
	output_execution_time();

	std::cout << "\nCSV file saved:\t" << output_filename << std::endl;
	return 0;
}