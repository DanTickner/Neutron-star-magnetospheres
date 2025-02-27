/*
cd OneDrive\PhD\Codes\20230810 Time evolution with updated equations
g++ Time_Evolution_15_traditional_VSH_for_B.cpp -o Time_Evolution_15_traditional_VSH_for_B --std=c++11
Time_Evolution_15_traditional_VSH_for_B
Evolve the magnetic and electric fields over time, by finite differencing the radial direction and using a VSH expansion in the theta direction.
	 

V15: Calculate B using the "normal VSH decomposition" instead of taking advantage of its divergencelessness. This gives more immediate accuracy; we may find issues at later times if div(B) grows without limit.
     Only compatible with Header_Time_Evolution_18_traditional_VSH_for_B.h.

*/

#include <iostream>
#include <iomanip>
#include <vector>
#include <math.h>
#include <fstream>
#include <chrono>

#include "Headers/Initial_conditions_03.h"
#include "../Vector_operations.h"
#include "Headers/Vector_B_rotating_dipole_ramp.h"
#include "Headers/Header_Time_Evolution_18_traditional_VSH_for_B.h"



int main(){
	
	//----- Failure conditions -----
	// Fail if cout gridpoints exceed number of gridpoints, to prevent code failure due to list overflow with no accompanying message.
	// Give warning if CFL condition violated at any gridpoint. Code will still run but results not trustworthy.
	if( cout_i > n_points_r ){
		std::cout << "Index for screen output in radial direction exceeds number of gridpoints (" << cout_i << " vs " << n_points_r <<"). Stopping the code." << std::endl;
		return 1;
	}
	if( cout_j > n_points_t ){
		std::cout << "Index for screen output in colatitude direction exceeds number of gridpoints (" << cout_j << " vs " << n_points_t <<"). Stopping the code." << std::endl;
		return 1;
	}
	if( delta_T > delta_r ){
		std::cout << "delta_T = " << delta_T <<" > delta_r = " << delta_r <<". This (albeit crude) application of the CFL condition has been violated. Recommend changing timestep or gridpoint size.\n" << std::endl;
	}
	
	
	//----- Build lists of gridpoints and precalculated basis functions -----
	calculate_gridpoints();
	calculate_associated_legendre_functions();
	calculate_sqrt_2Lplus1_over_4pi();
	
	
	
	//----- Output useful information to the screen -----
	// No need to estimate the size of the history CSV since it's far smaller than the profiles CSV.
	// For the future, perhaps split into multiple profiles CSVs. Then, add them all together to get an estimated total space required. A nice feature would be to fail if that much space isn't available, or if it exceeds some specified threshold.
	// NOTE: Number of rotation periods needs updating to include the ramp. Calculate on paper by integrating dOmega_by_dT.
	// Or add a variable rotation_angle which just gets Omega * delta_T added to it at each timestep, and can be recorded in the _history.csv file.
	
	std::cout << "Time start:\t" << time_start_string << std::endl;
	
	std::cout << "Parameter values (given in code units, then in SI units)" << std::endl;
	std::cout << "Omega_max              :\t" <<std::left<<std::setw(16)<< Omega_max<< std::endl;
	std::cout << "R_LC_max               :\t" <<std::left<<std::setw(16)<< 1.0 / Omega_max    << std::endl;
	std::cout << "r_min, r_max           :\t" <<std::left<<std::setw(16)<< r_min   <<std::left<<std::setw(16)<< r_max              << std::endl;
	std::cout << "n_points_r, delta_r    :\t" <<std::left<<std::setw(16)<< n_points_r <<std::left<<std::setw(16)<< delta_r << std::endl;
	std::cout << "n_points_t, delta_t    :\t" <<std::left<<std::setw(16)<< n_points_t <<std::left<<std::setw(16)<< delta_t << std::endl;
	std::cout << "Timestep               :\t" <<std::left<<std::setw(16)<< delta_T <<std::left<<std::setw(16)<< delta_T * T_factor << std::endl;
	std::cout << "CFL max timestep       :\t" <<std::left<<std::setw(16)<< delta_r <<std::left<<std::setw(16)<< delta_r * T_factor << std::endl;
	std::cout << "Length of simulation   :\t" <<std::left<<std::setw(16)<< T_max << T_max * T_factor <<"\tsteps: " << n_timesteps << std::endl;
	std::cout << "Est. CSV size (MB)     :\t" << estimate_csv_profiles_filesize() << std::endl;
	std::cout << "ell_max_rotation_off/on:\t" << ell_max_rotation_off <<"\t"<< ell_max_rotation_on << std::endl;
	std::cout << "Printing values at ( r[" << cout_i <<"], theta[" << cout_j <<"] ) = ( " << r[cout_i] <<", " << t[cout_j] << " ) = ( " << r[cout_i] <<", " << t[cout_j]/pi << " pi )." << std::endl;
	std::cout << "dOmega_by_dT_index     :\t" << dOmega_by_dT_index <<"\tdOmega_by_dT:\t"<< dOmega_by_dT << std::endl;
	std::cout << "Ramp start             :\t" << T_index_Omega_ramp_start <<"\t"<< T_index_Omega_ramp_start * delta_T << std::endl;
	std::cout << "Ramp stop              :\t" << T_index_Omega_ramp_stop  <<"\t"<< T_index_Omega_ramp_stop  * delta_T << std::endl;
	
	std::cout << "\nuse_outer_sponge_layer:\t" << use_outer_sponge_layer << std::endl;
	std::cout << "sigma_0, gamma, beta    :\t" << friction_sigma_0 <<"\t"<< friction_gamma <<"\t"<< friction_beta << std::endl;
	
	
	//----- Prepare output CSV and print headers to screen -----
	output_file_profiles  .open( "CSV/"  + output_filename + "_1_profiles.csv"    );
	output_file_history   .open( "CSV/"  + output_filename + "_2_history.csv"     );
	output_file_VSH_coeffs.open( "CSV/"  + output_filename + "_3_VSH_coeffs.csv"  );
	output_file_log       .open( "Logs/" + output_filename + ".txt"               );
	
	std::cout << "\nCSV file saved:\tCSV/"  << output_filename << "_1_profiles.csv"   << std::endl;
	std::cout <<   "CSV file saved:\tCSV/"  << output_filename << "_2_history.csv"    << std::endl;
	std::cout <<   "CSV file saved:\tCSV/"  << output_filename << "_3_VSH_coeffs.csv" << std::endl;
	std::cout <<   "Log file saved:\tLogs/" << output_filename << ".txt"              << std::endl;
	
	std::cout << "\nTo cut evolution short, press CTRL+C at any time. This is perfectly safe and the CSV file will still be saved up to that point.\n" << std::endl;
	
	output_headers_to_csv_profiles  ();
	output_headers_to_csv_history   ();
	output_headers_to_csv_VSH_coeffs();
	output_headers_to_screen();
	output_parameters_to_log_file();
	
	
	//----- Step at T=0 is subtly different -----
	determine_whether_rotation_on_and_set_ell_max_for_this_timestep( 0 );
	
	apply_initial_field_value_and_VSH_decomposition_for_B();

	calculate_stdev( 0 );
	calculate_force_free_conditions( 0 );
	
	calculate_radial_derivatives_of_VSH_coeffs();
	
	calculate_time_derivatives( 0, 0 );
	
	calculate_Cartesian_vector_components( 0 );
	evaluate_radial_derivatives();
	
	calculate_magnetic_energy( 0 );
	
	count_nans              ( 0 );
	output_to_screen        ( 0, 0 );
	output_to_csv_profiles  ( 0, 0 );
	output_to_csv_history   ( 0 );
	output_to_csv_VSH_coeffs( 0 );
	
	integrate_vectors_wrt_time_Adams_Bashforth_Order_1();
	apply_inner_boundary_conditions( 0 );
	apply_outer_boundary_conditions( 0 );
	
	
	
	
	//----- Evolve to future timesteps -----
	for( int T_index=1; T_index<=n_timesteps; T_index++ ){
		
		determine_whether_rotation_on_and_set_ell_max_for_this_timestep( T_index );
		
		int n = 1;	// Order of Adams-Bashforth method. Not really supposed to be a variable you can change, but used during code development.
	
		T += delta_T;
		if( ( T_index >= T_index_Omega_ramp_start ) and ( T_index <= T_index_Omega_ramp_stop ) ){
			Omega += dOmega_by_dT_index;
		}
		star_rotation_angle += Omega * delta_T * T_factor;
		
		apply_force_free_conditions_within_force_free_region( n );
		
		VSH_decomposition_for_E( n );
		VSH_decomposition_for_B( n );
		calculate_radial_derivatives_of_VSH_coeffs();
		
		calculate_stdev( n );
		calculate_force_free_conditions( n );
		
		calculate_time_derivatives( T_index, n );
		
		calculate_Cartesian_vector_components( n );
		evaluate_radial_derivatives();
		
		calculate_magnetic_energy( n );
		
		count_nans( n );
		output_to_screen        ( T_index, n );
		output_to_csv_profiles  ( T_index, n );
		output_to_csv_history   ( T_index );
		output_to_csv_VSH_coeffs( T_index );
		
		integrate_vectors_wrt_time_Adams_Bashforth_Order_2();
		apply_inner_boundary_conditions( n );
		apply_outer_boundary_conditions( n );
		
	}
	
	
	
	
	//----- Output execution time and finish code -----
	output_execution_time();
	
	return 0;
}