/*
cd OneDrive\PhD\Codes\20240311 Time evolution code
g++ Time_Evolution_28.cpp -o Time_Evolution --std=c++11
Time_Evolution
Evolve the magnetic and electric fields over time, by finite differencing the radial direction and using a VSH expansion in the theta direction.
	 
V28: Replace VSH decomposition by finite differencing in the theta direction.
     Remove all functionality for Chebyshev and VSH series.

*/

#include <iostream>
#include <iomanip>
#include <vector>
#include <math.h>
#include <fstream>
#include <string>
#include <sstream>
#include <chrono>

#include "Header_Initial_Conditions_10.h"
#include "Header_Time_Evolution_43.h"



int main(){
	
	//--- Check that no CSV output frequencies were accidentally set to zero ---
	if( csv_profiles_write_freq_r   == 0 ){ csv_profiles_write_freq_r   = 1; }
	if( csv_profiles_write_freq_t   == 0 ){ csv_profiles_write_freq_t   = 1; }
	if( csv_profiles_write_freq_T   == 0 ){ csv_history_write_freq_T    = 1; }
	if( csv_BCs_write_freq_T        == 0 ){ csv_BCs_write_freq_T        = 1; }
	if( csv_fields_write_freq_T     == 0 ){ csv_fields_write_freq_T     = 1; }
	if( csv_profiles_write_i_min    == 0 ){ csv_profiles_write_i_min    = 1; }
	if( csv_profiles_write_i_max    == 0 ){ csv_profiles_write_i_max    = 1; }
	if( csv_profiles_write_T_min    == 0 ){ csv_profiles_write_T_min    = 1; }
	if( csv_profiles_write_T_max    == 0 ){ csv_profiles_write_T_max    = 1; }
	
	
	//----- Prepare output CSV and output to screen and to log file -----
	if( output_csv_fields ){
		output_file_fields.open( output_filename_fields );
		output_headers_to_csv_fields();
	}
	if( output_csv_profiles ){
		output_file_profiles.open( output_filename_profiles );
		output_headers_to_csv_profiles();
	}
	if( output_csv_history ){
		output_file_history.open( output_filename_history );
		output_headers_to_csv_history();
	}
	if( output_csv_BCs ){
		output_file_BCs.open( output_filename_BCs );
		output_headers_to_csv_BCs();
	}
	
	output_file_log.open( output_filename_log );
	
	
	//----- Set time values, gridpoints and initial magnetic field -----
	if( load_fields_from_previous_run ){
		T_index = T_index_initial + 1;	// The entire evolution and post-analysis (BC application, etc) was performed at this step in the previous run, so we can go to the next.
		std::cout       << "Loading configuration from previous run\nFile:\t" << input_file_previous_run_base <<"\nT_index:\t"<< T_index_initial << std::endl;
		output_file_log << "Loading configuration from previous run\nFile:\t" << input_file_previous_run_base <<"\nT_index:\t"<< T_index_initial << "\n";
		set_time_and_rotation_values_from_CSV();
		set_gridpoints_from_CSV();
		apply_initial_field_values_from_CSV();
	}
	else {
		T_index_initial = 0;
		T_index         = 0;
		calculate_time_and_rotation_values();
		calculate_gridpoints();
		apply_initial_field_values_for_B();
		output_file_gridpoints .open( output_filename_gridpoints  );
		output_file_time_values.open( output_filename_time_values );
		output_to_csv_gridpoints ();
		output_to_csv_time_values();
	}
	
	
	//----- Failure conditions -----
	// Fail if cout gridpoints exceed number of gridpoints, to prevent code failure due to list overflow with no accompanying message.
	// Fail if ell_max exceeds n_t. I imagine this will cause issues with integration, but the code actually fails first on double estimate_csv_profiles_filesize().
	// Give warning if CFL condition violated at any gridpoint. Code will still run but results not trustworthy.
	if( cout_i > n_r ){
		std::cout       << "Index for screen output in radial direction exceeds number of gridpoints (" << cout_i << " vs " << n_r <<"). Stopping the code." << std::endl;
		output_file_log << "Index for screen output in radial direction exceeds number of gridpoints (" << cout_i << " vs " << n_r <<"). Stopping the code.\n";
		return 1;
	}
	if( cout_j > n_t ){
		std::cout       << "Index for screen output in colatitude direction exceeds number of gridpoints (" << cout_j << " vs " << n_t <<"). Stopping the code." << std::endl;
		output_file_log << "Index for screen output in colatitude direction exceeds number of gridpoints (" << cout_j << " vs " << n_t <<"). Stopping the code.\n";
		return 1;
	}
	
	calculate_CFL_max_timestep();
	
	if( delta_T > delta_T_CFL ){
		std::cout       << "delta_T = " << delta_T <<" > delta_T_CFL = " << delta_T_CFL <<". The CFL condition has been violated. Recommend decreasing timestep or increasing gridpoint size (decreasing number of gridpoints).\n" << std::endl;
		output_file_log << "delta_T = " << delta_T <<" > delta_T_CFL = " << delta_T_CFL <<". The CFL condition has been violated. Recommend decreasing timestep or increasing gridpoint size (decreasing number of gridpoints).\n";
	}
	
	
	//----- Build lists of gridpoints and precalculated basis functions, and calculate initial magnetic field -----
	calculate_functions_of_gridpoints();
	calculate_friction_for_sponge_layer();
	
	
	//----- Output chosen parameters and the headers row for the in-evolution monitoring -----
	output_parameters( std::cout );
	output_parameters( output_file_log );
	output_headers_to_screen_and_log_file();
	
	
	//----- Perform the time evolution -----
	
	while( T_index < T_index_max ){
		
		
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
			case 3:
				integrate_vectors_wrt_time_Adams_Bashforth_Order_3();
				break;
			default:
				integrate_vectors_wrt_time_Adams_Bashforth_Order_4();
		}
		
		
		//--- Apply the remaining functions ---
		apply_inner_boundary_conditions();
		apply_outer_boundary_conditions();
		
		ramp_up_electric_fields_for_rotation();
		
		calculate_dot_products                              ();
		apply_force_free_conditions_within_force_free_region();
		
		for( int f=0; f<3; f++ ){
			B_dr[f] = f_dx_FD( B[f], delta_r );
			B_dt[f] = f_dy_FD( B[f], delta_t );
			E_dr[f] = f_dx_FD( E[f], delta_r );
			E_dt[f] = f_dy_FD( E[f], delta_t );
		}
		
		B_div   = divergence_from_FD_derivatives( B, B_dr, B_dt );
		E_div   = divergence_from_FD_derivatives( E, E_dr, E_dt );
		B_curl  = curl_from_FD_derivatives      ( B, B_dr, B_dt );
		E_curl  = curl_from_FD_derivatives      ( E, E_dr, E_dt );
		
		calculate_current_density();
		calculate_time_derivatives_of_vector_components();
		
		calculate_Cartesian_vector_components();
		
		calculate_electric_and_magnetic_energy();
		calculate_normalised_B_div();
		count_nans                ();
		
		
		//--- Output ---
		if( output_csv_fields ){
			output_to_csv_fields();
		}
		if( output_csv_profiles ){
			output_to_csv_profiles();
		}
		if( output_csv_history ){
			output_to_csv_history();
		}
		if( output_csv_BCs ){
			output_to_csv_BCs();
		}
		
		output_to_screen_and_log_file();
		
		// If need profile at a specific timestep, use the template below. For multiple specific timesteps, just copy-and-paste and change the value of T_index.
		//if( T_index == 123 ){ output_to_csv_profiles( true ); }
		
		
		//--- Stop evolution if the code has blown up ---
		/*
		// Some C++ compilers require std::isnan and std::isinf, whereas some require isnan and isinf. Comment-out as appropriate.
		//if( ( std::isnan( total_magnetic_energy ) ) or ( std::isinf( total_magnetic_energy ) ) ){
		if( ( isnan( total_magnetic_energy ) ) or ( isinf( total_magnetic_energy ) ) ){
			std::cout       << "Magnetic energy is no longer finite. The evolution has probably blown up. Stopping the code." << std::endl;
			output_file_log << "Magnetic energy is no longer finite. The evolution has probably blown up. Stopping the code.\n";
			return 1;
		}
		*/
		
		//--- Increment the timestep ---
		T_index ++;
		
	}
	
	
	
	
	//----- Output execution time and finish code -----
	output_execution_time();
	
	return 0;
}