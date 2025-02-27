/*
cd OneDrive\PhD\Codes\20240311 Time evolution code
g++ Time_Evolution_27.cpp -o Time_Evolution --std=c++11
Time_Evolution
Evolve the magnetic and electric fields over time, by finite differencing the radial direction and using a VSH expansion in the theta direction.
	 
V27: Add functionality to continue evolution from fields saved at a previous timestep. To do this, we need to set the gridpoints to be those used in the previous evolution. This also means redeclaring all our array variables to have the right dimensions n_r,n_t.
     Replace n_T by T_index_max.
	 Change void output_headers_to_screen_and_log_file() to void output_headers_to_screen_and_log_file() and void output_to_screen_and_log_file() to void output_to_screen_and_log_file().
	 Output the chosen-gridpoint values to the log file within these functions as well as the console. This is useful when doing back-to-back runs with similar parameters, so you can quickly see whether there is a big effect before waiting for the run to finish.

*/

#include <iostream>
#include <iomanip>
#include <vector>
#include <math.h>
#include <fstream>
#include <string>
#include <sstream>
#include <chrono>

#include "Header_Initial_Conditions_09.h"
#include "Header_Time_Evolution_42.h"



int main(){
	
	//--- Check that no CSV output frequencies were accidentally set to zero ---
	if( csv_profiles_write_freq_r   == 0 ){ csv_profiles_write_freq_r   = 1; }
	if( csv_profiles_write_freq_t   == 0 ){ csv_profiles_write_freq_t   = 1; }
	if( csv_profiles_write_freq_T   == 0 ){ csv_history_write_freq_T    = 1; }
	if( csv_VSH_coeffs_write_freq_T == 0 ){ csv_VSH_coeffs_write_freq_T = 1; }
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
	if( output_csv_VSH_coeffs ){
		output_file_VSH_coeffs.open( output_filename_VSH_coeffs );
		output_headers_to_csv_VSH_coeffs();
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
	if( ell_max >= n_t ){
		std::cout       << "Maximum value of ell for VSH decomposition exceeds number of gridpoints (" << ell_max << " vs " << n_t <<"). Stopping the code." << std::endl;
		output_file_log << "Maximum value of ell for VSH decomposition exceeds number of gridpoints (" << ell_max << " vs " << n_t <<"). Stopping the code.\n";
		return 1;
	}
	
	calculate_CFL_max_timestep();
	
	if( delta_T > delta_T_CFL ){
		std::cout       << "delta_T = " << delta_T <<" > delta_T_CFL = " << delta_T_CFL <<". The CFL condition has been violated. Recommend decreasing timestep or increasing gridpoint size (decreasing number of gridpoints).\n" << std::endl;
		output_file_log << "delta_T = " << delta_T <<" > delta_T_CFL = " << delta_T_CFL <<". The CFL condition has been violated. Recommend decreasing timestep or increasing gridpoint size (decreasing number of gridpoints).\n";
	}
	
	
	//----- Build lists of gridpoints and precalculated basis functions, and calculate initial magnetic field -----
	calculate_functions_of_gridpoints();
	calculate_associated_legendre_functions();
	calculate_sqrt_2Lplus1_over_4pi();
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
			default:
				integrate_vectors_wrt_time_Adams_Bashforth_Order_3();
		}
		
		
		//--- Apply the remaining functions ---
		ramp_up_electric_fields_for_rotation();
		
		apply_inner_boundary_conditions();
		apply_outer_boundary_conditions();
		calculate_dot_products();	 //Redo the calculation because the above functions will change E and B.
		
		check_force_free_conditions();
		calculate_n_gridpoints_within_FF_region();
		apply_force_free_conditions_within_force_free_region();
		calculate_dot_products();	 //Redo the calculation because the above functions will change E and B.
		
		B_VSH_coeffs[0] = VSH_decomposition_r( B );
		B_VSH_coeffs[2] = VSH_decomposition_2( B );
		E_VSH_coeffs[0] = VSH_decomposition_r( E );
		E_VSH_coeffs[1] = VSH_decomposition_1( E );
		E_VSH_coeffs[2] = VSH_decomposition_2( E );
		
		for( int ell=0; ell<=ell_max; ell++ ){
			B_VSH_coeffs_dr[0][ell] = radial_derivatives_1_FD( B_VSH_coeffs[0][ell] );
			B_VSH_coeffs_dr[2][ell] = radial_derivatives_1_FD( B_VSH_coeffs[2][ell] );
			E_VSH_coeffs_dr[0][ell] = radial_derivatives_1_FD( E_VSH_coeffs[0][ell] );
			E_VSH_coeffs_dr[1][ell] = radial_derivatives_1_FD( E_VSH_coeffs[1][ell] );
			E_VSH_coeffs_dr[2][ell] = radial_derivatives_1_FD( E_VSH_coeffs[2][ell] );
		}
		
		//- Calculate VSH series of B by normal method (comment-out as appropriate) -
		B_VSH_coeffs[1] = VSH_decomposition_1( B );
		for( int ell=0; ell<=ell_max; ell++ ){
			B_VSH_coeffs_dr[1][ell] = radial_derivatives_1_FD( B_VSH_coeffs[1][ell] );
		}
		
		/*
		//- Calculate VSH series of B by divergenceless method (comment-out as appropriate )-
		for( int ell=0; ell<=ell_max; ell++ ){
			B_VSH_coeffs_dr2[0][ell] = radial_derivatives_2_FD( B_VSH_coeffs[0][ell] );
		}
		calculate_B_1_ell_and_derivative_divergenceless();
		*/
		
		//--- Remaining operations at each step ---
		calculate_stdev_of_VSH_decomposition();
		
		calculate_radial_derivatives_of_vector_components();
		
		B_div  = divergence_from_VSH_coeffs( B_VSH_coeffs, B_VSH_coeffs_dr );
		E_div  = divergence_from_VSH_coeffs( E_VSH_coeffs, E_VSH_coeffs_dr );
		B_curl = curl_from_VSH_coeffs      ( B_VSH_coeffs, B_VSH_coeffs_dr );
		E_curl = curl_from_VSH_coeffs      ( E_VSH_coeffs, E_VSH_coeffs_dr );
		
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
		if( output_csv_VSH_coeffs ){
			output_to_csv_VSH_coeffs();
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