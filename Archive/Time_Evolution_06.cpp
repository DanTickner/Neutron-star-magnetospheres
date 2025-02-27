/*
cd OneDrive\PhD\Codes\20230810 Time evolution with updated equations
g++ Time_Evolution_06.cpp -o Time_Evolution_06
Time_Evolution_06
Evolve the magnetic and electric fields over time, by finite differencing the radial direction and using a VSH expansion in the theta direction.
	 

V06: Support function renaming
     void apply_second_force_free_condition_within_light_cylinder() -> void apply_second_force_free_condition_within_force_free_region()
	 Remove option for Chebyshev decomposition for now, as it seems very inaccurate.
	 Replace crude CSV filesize estimate by a function which more carefully counts the number rows, and can handle choices of min/max values of r,T to output.


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
//#include "Headers/Vector_E_rotating_dipole.h"
#include "Headers/Vector_E_zero.h"
//#include "Headers/Vector_B_Michel_monopole.h"
//#include "Headers/Vector_E_rotating_Michel_monopole.h"


//----- Define variables -----
// Define these before including the main header file so that they are known to the functions within it.
int n_points_r  = 1000;
int n_points_t  = 500;
int n_timesteps = 600;

double delta_T = 0.015;

int w = 14;			// fixed width for text output


// Choose output filename, including parent file structure. Don't include ".csv" extension because this will be modified to add _profiles and _history.
std::string output_filename = "CSV/20230921_zero_longer";

int csv_profiles_write_freq_r = n_points_r/50;	// Write to profiles CSV after this many steps in the radial direction.
int csv_profiles_write_freq_t = n_points_t/50;	// Write to profiles CSV after this many steps in the polar  direction.
int csv_profiles_write_freq_T = 10;				// Write to profiles CSV after this many timesteps.
int csv_history_write_freq_T  = 1;				// Write to history  CSV after this many timesteps.

int csv_profiles_write_i_min = 0;				// Only write to profiles CSV if i in this range. Useful for troubleshooting problematic radial intervals.
int csv_profiles_write_i_max = n_points_r;		// Only write to profiles CSV if i in this range. Useful for troubleshooting problematic radial intervals.
int csv_profiles_write_T_min = 0;				// Only write to profiles CSV if T_index in this range. Useful for troubleshooting problematic time intervals.
int csv_profiles_write_T_max = n_timesteps;		// Only write to profiles CSV if T_index in this range. Useful for troubleshooting problematic time intervals.

int cout_i = 999;
int cout_j = 250;
int cout_freq_T = 1;

#include "Headers/Header_Time_Evolution_07.h"



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
	calculate_n_gridpoints_within_force_free_region();
	
	
	
	//----- Output useful information to the screen -----
	// No need to estimate the size of the history CSV since it's far smaller than the profiles CSV.
	// For the future, perhaps split into multiple profiles CSVs. Then, add them all together to get an estimated total space required. A nice feature would be to fail if that much space isn't available, or if it exceeds some specified threshold.
	double T_max = n_timesteps * delta_T;
	
	std::cout << "Parameter values (given in code units, then in SI units)" << std::endl;
	std::cout << "Omega               :\t" <<std::left<<std::setw(16)<< Omega   << std::endl;
	std::cout << "R_LC                :\t" <<std::left<<std::setw(16)<< R_LC    << std::endl;
	std::cout << "r_min, r_max        :\t" <<std::left<<std::setw(16)<< r_min   <<std::left<<std::setw(16)<< r_max              << std::endl;
	std::cout << "Timestep            :\t" <<std::left<<std::setw(16)<< delta_T <<std::left<<std::setw(16)<< delta_T * T_factor << std::endl;
	std::cout << "CFL max timestep    :\t" <<std::left<<std::setw(16)<< delta_r <<std::left<<std::setw(16)<< delta_r * T_factor << std::endl;
	std::cout << "Length of simulation:\t" <<std::left<<std::setw(16)<< T_max << T_max * T_factor <<"\tRotation periods: "<< T_max / ( 2.0/Omega ) <<"\tsteps: " << n_timesteps << std::endl;
	std::cout << "Est. CSV size (MB)  :\t" << estimate_csv_profiles_filesize() << std::endl;
	std::cout << "n_max = " << n_max <<"\tell_max = " << ell_max  << std::endl;
	std::cout << "Printing values at ( r[" << cout_i <<"], theta[" << cout_j <<"] ) = ( " << r[cout_i] <<", " << t[cout_j] << " ) = ( " << r[cout_i] <<", " << t[cout_j]/pi << " pi )." << std::endl;
	
	std::cout << "\nTo cut evolution short, press CTRL+C at any time. This is perfectly safe and the CSV file will still be saved up to that point.\n" << std::endl;
	
	
	
	
	//----- Prepare output CSV and print headers to screen -----
	output_file_profiles.open( output_filename + "_1_profiles.csv" );
	output_file_history .open( output_filename + "_2_history.csv"  );
	output_headers_to_csv_profiles();
	output_headers_to_csv_history ();
	output_headers_to_screen();
	
	
	
	//----- Step at T=0 is subtly different -----
	
	time_start_seconds = std::chrono::high_resolution_clock::now().time_since_epoch().count() * 1e-9;
	
	apply_initial_field_value_and_VSH_decomposition_for_B();
	apply_initial_field_value_and_VSH_decomposition_for_E();
	
	apply_second_force_free_condition_within_force_free_region();
	
	VSH_decomposition_for_E();	// B cannot be altered by application of the first FF condition, so no need to recalculate its VSH series in the first step.
	calculate_stdev();
	calculate_force_free_conditions();
	
	calculate_radial_derivatives_of_VSH_coeffs();
	calculate_B_1_L_and_B_1_L_dr_initial();
	
	calculate_time_derivatives();
	
	calculate_Cartesian_vector_components();
	evaluate_radial_derivatives();
	
	calculate_magnetic_energy();
	
	output_to_screen      ( 0 );
	output_to_csv_profiles( 0 );
	output_to_csv_history ( 0 );
	
	integrate_vectors_wrt_time();
	apply_inner_boundary_conditions();
	apply_outer_boundary_conditions();
	
	
	
	
	//----- Evolve to future timesteps -----
	for( int T_index=1; T_index<=n_timesteps; T_index++ ){
	
		T += delta_T;
		
		apply_second_force_free_condition_within_force_free_region();
		
		VSH_decomposition_for_B();
		VSH_decomposition_for_E();
		calculate_stdev();
		calculate_force_free_conditions();
		
		calculate_radial_derivatives_of_VSH_coeffs();
		calculate_B_1_L_and_B_1_L_dr();
		
		calculate_time_derivatives();
		
		calculate_Cartesian_vector_components();
		evaluate_radial_derivatives();
		
		calculate_magnetic_energy();
		
		output_to_screen      ( T_index );
		output_to_csv_profiles( T_index );
		output_to_csv_history ( T_index );
		
		integrate_vectors_wrt_time();
		apply_inner_boundary_conditions();
		apply_outer_boundary_conditions();
		
	}
	
	
	
	
	//----- Output execution time and finish code -----
	output_execution_time();
	
	std::cout << "\nCSV file saved:\t" << output_filename << "_1_profiles.csv" << std::endl;
	std::cout <<   "CSV file saved:\t" << output_filename << "_2_history.csv"  << std::endl;
	
	return 0;
}