/*
cd OneDrive\PhD\Codes\20240311 Time evolution code
g++ Time_Evolution_26.cpp -o Time_Evolution --std=c++11
Time_Evolution
Evolve the magnetic and electric fields over time, by finite differencing the radial direction and using a VSH expansion in the theta direction.
	 
V26: Split the VSH decomposition and radial derivative calculation into component-wise functions. See OneNote 20240329Fri.
	 Change while( T_index =< n_T ) to while( T_index < n_T ) because the former gave one too many timesteps.
	 Move output file creation before generation of gridpoints etc and before failure conditions so that everything can be recorded.
	 Combine void output_parameters_to_screen() and void output_parameters_to_log_file() into a single file with an ofstream argument, saving code repetition and guaranteeing consistency.
	 void calculate_normalised_B_div() to track the size of div(B) across the domain over time.
	 Add CSV _5_fields containing the magnetic and electric field values at all gridpoints, to be used when reading-in values for rerunning simulations with e.g. higher output frequency. 

*/

#include <iostream>
#include <iomanip>
#include <vector>
#include <math.h>
#include <fstream>
#include <string>
#include <sstream>
#include <chrono>

#include "Header_Initial_Conditions_08.h"
#include "Header_Time_Evolution_40.h"



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
	if( output_csv_fields ){
		output_file_fields.open( output_filename_fields );
		output_headers_to_csv_fields();
	}
	
	output_file_log.open( output_filename_log );
	
	
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
	calculate_time_and_rotation_values();
	calculate_gridpoints();
	calculate_associated_legendre_functions();
	calculate_sqrt_2Lplus1_over_4pi();
	calculate_friction_for_sponge_layer();
	apply_initial_field_values_for_B();
	
	output_file_gridpoints.open( output_filename_gridpoints );
	output_to_csv_gridpoints();
	std::cout << r[3] << std::endl;
	
	
	//----- Output chosen parameters and the headers row for the in-evolution monitoring -----
	output_parameters( std::cout );
	output_parameters( output_file_log );
	output_headers_to_screen();
	
	
	//----- Perform the time evolution -----
	
	while( T_index < n_T ){
		
		
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
		
		/*
		//- Calculate VSH series of B by normal method (comment-out as appropriate) -
		B_VSH_coeffs[1] = VSH_decomposition_1( B );
		for( int ell=0; ell<=ell_max; ell++ ){
			B_VSH_coeffs_dr[1][ell] = radial_derivatives_1_FD( B_VSH_coeffs[1][ell] );
		}
		*/
		
		//- Calculate VSH series of B by divergenceless method (comment-out as appropriate )-
		for( int ell=0; ell<=ell_max; ell++ ){
			B_VSH_coeffs_dr2[0][ell] = radial_derivatives_2_FD( B_VSH_coeffs[0][ell] );
		}
		calculate_B_1_ell_and_derivative_divergenceless();
		
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
		if( output_csv_fields ){
			output_to_csv_fields();
		}
		
		output_to_screen();
		
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