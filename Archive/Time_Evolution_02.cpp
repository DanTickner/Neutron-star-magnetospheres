/*
cd OneDrive\PhD\Codes\20230810 Time evolution with updated equations
g++ Time_Evolution_02.cpp -o Time_Evolution_02
Time_Evolution_02
Evolve the magnetic and electric fields over time, by finite differencing the radial direction and using a VSH expansion in the theta direction.
	 

V02: Don't try and massage the VSH coefficients. Just evaluate the field at every step, do the evolution, apply the FF conditions, and then redo the decomposition.


The process at each timestep is as follows.
1. Apply the second force-free condition.
2. Perform a VSH decomposition of the magnetic and electric fields. (Do after 1 because 1 potentially changes the fields.)
3. Recalculate the radial derivatives of the VSH coefficients of the fields.
4. Use the VSH coefficients to calculate the spatial derivatives of the fields and hence the time derivatives of the fields.
5. Output to screen and CSV. Do this before integration because integration moves the values to those they would take at the next timestep.
6. Perform the integration.

The integrated fields after the very last timestep are never output.
The code is intended to have very many timesteps, so this is not a problem. Just add an extra timestep if you desire the fields at a certain time.
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
int n_points_r  = 1000;
int n_points_t  = 500;
int n_timesteps = 150;

double delta_T = 0.001;

int w = 14;			// fixed width for text output


std::string output_filename = "CSV/20230818_inside_LC.csv";

int csv_write_freq_r = n_points_r/50;	// Write to CSV after this many steps in the radial direction.
int csv_write_freq_t = n_points_t/50;	// Write to CSV after this many steps in the polar  direction.
int csv_write_freq_T = 10;				// Write to CSV after this many timesteps.

int cout_i = 4;
int cout_j = 7;
int cout_freq_T = 10;

#include "Headers/Header_Time_Evolution_03.h"



int main(){
	
	//----- Failure conditions -----
	// Fail if cout gridpoints exceed number of gridpoints, to prevent code failure due to list overflow with no accompanying message.
	// Fail if CFL condition violated at any gridpoint. Code would still run but results not trustworthy.
	if( cout_i > n_points_r ){
		std::cout << "Index for screen output in radial direction exceeds number of gridpoints (" << cout_i << " vs " << n_points_r <<"). Stopping the code." << std::endl;
		return 1;
	}
	if( cout_j > n_points_t ){
		std::cout << "Index for screen output in colatitude direction exceeds number of gridpoints (" << cout_j << " vs " << n_points_t <<"). Stopping the code." << std::endl;
		return 1;
	}
	if( delta_T > delta_r ){
		std::cout << "delta_T = " << delta_T <<" > delta_r = " << delta_r <<". This (albeit crude) application of the CFL condition has been violated. Stopping the code." << std::endl;
		return 1;
	}
	
	
	//----- Build lists of gridpoints and precalculated basis functions -----
	calculate_gridpoints();
	calculate_associated_legendre_functions();
	calculate_sqrt_2Lplus1_over_4pi();
	calculate_n_gridpoints_inside_light_cylinder();
	
	
	
	//----- Output useful information to the screen -----
	double csv_filesize_factor = 1.9e-4;	// Adjust to give more accurate esimates.
	double csv_filesize_guess  = csv_filesize_factor * n_timesteps * n_points_r * n_points_t / ( (double) csv_write_freq_T * csv_write_freq_r * csv_write_freq_t );
	
	std::cout << "Parameter values (given in code units, then in SI units)" << std::endl;
	std::cout << "Omega               :\t" <<std::left<<std::setw(16)<< Omega   << std::endl;
	std::cout << "R_LC                :\t" <<std::left<<std::setw(16)<< R_LC    << std::endl;
	std::cout << "r_min, r_max        :\t" <<std::left<<std::setw(16)<< r_min   <<std::left<<std::setw(16)<< r_max              << std::endl;
	std::cout << "Timestep            :\t" <<std::left<<std::setw(16)<< delta_T <<std::left<<std::setw(16)<< delta_T * T_factor << std::endl;
	std::cout << "CFL max timestep    :\t" <<std::left<<std::setw(16)<< delta_r <<std::left<<std::setw(16)<< delta_r * T_factor << std::endl;
	std::cout << "Length of simulation:\t" <<std::left<<std::setw(16)<< n_timesteps * delta_T <<  n_timesteps * delta_T * T_factor <<"\tsteps: " << n_timesteps << std::endl;
	std::cout << "Est. CSV size (MB)  :\t" <<  csv_filesize_guess << std::endl;
	std::cout << "n_max = " << n_max <<"\tell_max = " << ell_max  << std::endl;
	std::cout << "Printing values at ( r[" << cout_i <<"], theta[" << cout_j <<"] ) = ( " << r[cout_i] <<", " << t[cout_j] << " )." << std::endl;
	
	std::cout << "\nTo cut evolution short, press CTRL+C at any time. This is perfectly safe and the CSV file will still be saved up to that point.\n" << std::endl;
	
	
	
	
	//----- Prepare output CSV and print headers to screen -----
	output_file.open( output_filename );
	output_headers_to_csv();
	output_headers_to_screen();
	
	
	
	//----- Step at T=0 is subtly different -----
	
	time_start_seconds = std::chrono::high_resolution_clock::now().time_since_epoch().count() * 1e-9;
	
	apply_initial_field_values_and_VSH_decompositions();
	
	apply_second_force_free_condition_within_light_cylinder();
	
	VSH_decomposition_for_E();	// B cannot be altered by application of the first FF condition, so no need to recalculate its VSH series in the first step.
	calculate_force_free_conditions();
	
	calculate_radial_derivatives_of_VSH_coeffs_inner_boundary();
	calculate_radial_derivatives_of_VSH_coeffs();
	calculate_B_1_L_and_B_1_L_dr_inner_boundary();
	calculate_B_1_L_and_B_1_L_dr_initial();
	
	calculate_time_derivatives();
	
	calculate_Cartesian_vector_components();
	evaluate_radial_derivatives();
	
	output_to_screen( 0 );
	output_to_csv   ( 0 );
	
	integrate_vectors_wrt_time();
	apply_outer_boundary_conditions();
	
	
	
	
	//----- Evolve to future timesteps -----
	for( int T_index=1; T_index<=n_timesteps; T_index++ ){
	
		T += delta_T;
		
		apply_second_force_free_condition_within_light_cylinder();
		
		VSH_decomposition_for_B();
		VSH_decomposition_for_E();
		calculate_force_free_conditions();
		
		calculate_radial_derivatives_of_VSH_coeffs();
		calculate_B_1_L_and_B_1_L_dr();
		
		calculate_time_derivatives();
		
		calculate_Cartesian_vector_components();
		evaluate_radial_derivatives();
		
		output_to_screen( T_index );
		output_to_csv   ( T_index );
		//if( T_index >= 80 ){ output_to_csv   ( T_index ); }
		
		integrate_vectors_wrt_time();
		apply_outer_boundary_conditions();
		
	}
	
	
	
	
	//----- Output execution time and finish code -----
	output_execution_time();
	
	std::cout << "\nCSV file saved:\t" << output_filename << std::endl;
	return 0;
}