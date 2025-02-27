/*
cd OneDrive\PhD\Codes\20240311 Time evolution code
g++ Time_Evolution_23.cpp -o Time_Evolution_23 --std=c++11
Time_Evolution_23
Evolve the magnetic and electric fields over time, by finite differencing the radial direction and using a VSH expansion in the theta direction.
	 
V23: Move all steps to a single loop, with an if-statement choosing the Adams-Bashforth function to use.

*/

#include <iostream>
#include <iomanip>
#include <vector>
#include <math.h>
#include <fstream>
#include <chrono>

#include "Header_Initial_Conditions_07.h"
#include "Header_Time_Evolution_36.h"



int main(){
	
	//----- Failure conditions -----
	// Fail if cout gridpoints exceed number of gridpoints, to prevent code failure due to list overflow with no accompanying message.
	// Fail if ell_max exceeds n_points_t. I imagine this will cause issues with integration, but the code actually fails first on double estimate_csv_profiles_filesize().
	// Give warning if CFL condition violated at any gridpoint. Code will still run but results not trustworthy.
	if( cout_i > n_points_r ){
		std::cout << "Index for screen output in radial direction exceeds number of gridpoints (" << cout_i << " vs " << n_points_r <<"). Stopping the code." << std::endl;
		return 1;
	}
	if( cout_j > n_points_t ){
		std::cout << "Index for screen output in colatitude direction exceeds number of gridpoints (" << cout_j << " vs " << n_points_t <<"). Stopping the code." << std::endl;
		return 1;
	}
	if( ell_max >= n_points_t ){
		std::cout << "Maximum value of ell for VSH decomposition exceeds number of gridpoints (" << ell_max << " vs " << n_points_t <<"). Stopping the code." << std::endl;
		return 1;
	}
	
	calculate_CFL_max_timestep();
	
	if( delta_T > delta_T_CFL ){
		std::cout << "delta_T = " << delta_T <<" > delta_T_CFL = " << delta_T_CFL <<". The CFL condition has been violated. Recommend decreasing timestep or increasing gridpoint size (decreasing number of gridpoints).\n" << std::endl;
	}
	
	
	//----- Build lists of gridpoints and precalculated basis functions, and calculate initial magnetic field -----
	calculate_gridpoints();
	calculate_associated_legendre_functions();
	calculate_sqrt_2Lplus1_over_4pi();
	calculate_time_and_rotation_values();
	calculate_friction_for_sponge_layer();
	apply_initial_field_values_for_B();
	
	
	
	
	//----- Prepare output CSV and output to screen and to log file -----
	output_file_profiles  .open( output_filename_profiles   );
	output_file_history   .open( output_filename_history    );
	output_file_VSH_coeffs.open( output_filename_VSH_coeffs );
	output_file_BCs       .open( output_filename_BCs        );
	output_file_log       .open( output_filename_log        );
	
	output_headers_to_csv_profiles  ();
	output_headers_to_csv_history   ();
	output_headers_to_csv_VSH_coeffs();
	output_headers_to_csv_BCs       ();
	
	output_parameters_to_screen     ();
	output_parameters_to_log_file   ();
	output_headers_to_screen        ();
	
	
	
	
	//----- Perform the time evolution -----
	
	while( T_index <= n_T ){
		
		
		//--- Integrate the fields ---
		switch( T_index ){
			case 0:
				// Don't integrate.
				break;
			case 1:
				integrate_vectors_wrt_time_Adams_Bashforth_Order_1();
				break;
			case 2:
				integrate_vectors_wrt_time_Adams_Bashforth_Order_2();
				break;
			default:
				integrate_vectors_wrt_time_Adams_Bashforth_Order_3();
		}
		
		
		//--- Apply the remaining functions ---
		apply_inner_boundary_conditions();
		apply_outer_boundary_conditions();
		
		ramp_up_electric_fields_for_rotation();
		
		calculate_dot_products                              ();
		apply_force_free_conditions_within_force_free_region();
		
		VSH_decomposition( E, E_VSH_coeffs );
		
		VSH_decomposition_divergenceless( B, B_VSH_coeffs );
		calculate_radial_derivatives_of_VSH_coeffs_FD();
		calculate_B_1_L_and_B_1_L_dr();
		
		calculate_stdev_of_VSH_decomposition();
		
		calculate_radial_derivatives_of_vector_components();
		calculate_time_derivatives_of_vector_components  ();
		
		calculate_Cartesian_vector_components();
		
		calculate_magnetic_energy();
		count_nans               ();
		
		
		//--- Output ---
		output_to_screen        ();
		output_to_csv_profiles  ();
		output_to_csv_history   ();
		output_to_csv_VSH_coeffs();
		output_to_csv_BCs       ();
		
		// If need profile at a specific timestep, use the template below. For multiple specific timesteps, just copy-and-paste and change the value of T_index.
		//if( T_index == 123 ){ output_to_csv_profiles( true ); }
		
		
		//--- Increment the timestep ---
		T_index ++;
		
	}
	
	
	
	
	//----- Output execution time and finish code -----
	output_execution_time();
	
	return 0;
}