/*
cd OneDrive\PhD\Codes\20230810 Time evolution with updated equations
g++ Time_Evolution_18.cpp -o Time_Evolution_18 --std=c++11
Time_Evolution_18
Evolve the magnetic and electric fields over time, by finite differencing the radial direction and using a VSH expansion in the theta direction.

Do not use with earlier version than Header_Time_Evolution_24.h.
	 
V18: Using E_ramp_dT is inaccurate. Instead, add to the electric field at each timestep after integration has been performed, with void ramp_up_electric_fields_for_rotation(). See OneNote notes 20231020Fri.

*/

#include <iostream>
#include <iomanip>
#include <vector>
#include <math.h>
#include <fstream>
#include <chrono>

#include "Headers/Initial_conditions_04.h"
#include "../Vector_operations.h"
#include "Headers/Vector_B_rotating_dipole_ramp.h"
#include "Headers/Header_Time_Evolution_29.h"



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
	calculate_friction_for_sponge_layer();
	
	
	
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
	std::cout << "ell_max                :\t" << ell_max << std::endl;
	std::cout << "Printing values at ( r[" << cout_i <<"], theta[" << cout_j <<"] ) = ( " << r[cout_i] <<", " << t[cout_j] << " ) = ( " << r[cout_i] <<", " << t[cout_j]/pi << " pi )." << std::endl;
	std::cout << "dOmega_by_dT_index     :\t" << dOmega_by_dT_index <<"\tdOmega_by_dT:\t"<< dOmega_by_dT << std::endl;
	std::cout << "Ramp start             :\t" << T_index_Omega_ramp_start <<"\t"<< T_index_Omega_ramp_start * delta_T << std::endl;
	std::cout << "Ramp stop              :\t" << T_index_Omega_ramp_stop  <<"\t"<< T_index_Omega_ramp_stop  * delta_T << std::endl;
	
	std::cout << "\nuse_outer_sponge_layer:\t" << use_outer_sponge_layer << std::endl;
	std::cout << "sigma_0, gamma, beta    :\t" << friction_sigma_0 <<"\t"<< friction_gamma <<"\t"<< friction_beta << std::endl;
	
	
	//----- Prepare output CSV and print headers to screen -----
	output_file_profiles  .open( "CSV/"  + output_filename + "_1_profiles.csv"   );
	output_file_history   .open( "CSV/"  + output_filename + "_2_history.csv"    );
	output_file_VSH_coeffs.open( "CSV/"  + output_filename + "_3_VSH_coeffs.csv" );
	output_file_BCs       .open( "CSV/"  + output_filename + "_4_BCs.csv"        );
	output_file_log       .open( "Logs/" + output_filename + ".txt"              );
	
	std::cout << "\nCSV file saved:\tCSV/"  << output_filename << "_1_profiles.csv"   << std::endl;
	std::cout <<   "CSV file saved:\tCSV/"  << output_filename << "_2_history.csv"    << std::endl;
	std::cout <<   "CSV file saved:\tCSV/"  << output_filename << "_3_VSH_coeffs.csv" << std::endl;
	std::cout <<   "CSV file saved:\tCSV/"  << output_filename << "_4_BCs.csv"        << std::endl;
	std::cout <<   "Log file saved:\tLogs/" << output_filename << ".txt"              << std::endl;
	
	std::cout << "\nTo cut evolution short, press CTRL+C at any time. This is perfectly safe and the CSV file will still be saved up to that point.\n" << std::endl;
	
	output_headers_to_csv_profiles  ();
	output_headers_to_csv_history   ();
	output_headers_to_csv_VSH_coeffs();
	output_headers_to_csv_BCs       ();
	output_headers_to_screen        ();
	output_parameters_to_log_file   ();
	
	
	
	
	//----- T_index = 0: No integration but apply BCs and FF conditions -----
	
	apply_initial_field_values_for_B();
	
	apply_inner_boundary_conditions( true );
	apply_outer_boundary_conditions( true );
	
	apply_force_free_conditions_within_force_free_region();
	
	VSH_decomposition_for_B();
	VSH_decomposition_for_E();
	
	calculate_radial_derivatives_of_VSH_coeffs       ();
	calculate_B_1_L_and_B_1_L_dr                     ();
	calculate_radial_derivatives_of_vector_components();
	calculate_dot_products                           ();
	calculate_time_derivatives_of_vector_components  ();
	
	calculate_Cartesian_vector_components();
	
	calculate_stdev          ();
	calculate_magnetic_energy();
	count_nans               ();
	
	output_to_screen        ();
	output_to_csv_profiles  ();
	output_to_csv_history   ();
	output_to_csv_VSH_coeffs();
	output_to_csv_BCs       ();
	
	
	
	//----- Step at T_index = 1: First-order Adams-Bashforth technique -----
	T_index = 1;
	T += delta_T;
	if( ( T_index >= T_index_Omega_ramp_start ) and ( T_index <= T_index_Omega_ramp_stop ) ){
		Omega += dOmega_by_dT_index;
		star_rotation_angle += Omega * delta_T * T_factor;
	}
	
	integrate_vectors_wrt_time_Adams_Bashforth_Order_1();
	
	apply_inner_boundary_conditions( true );
	apply_outer_boundary_conditions( true );
	
	ramp_up_electric_fields_for_rotation();
	
	calculate_dot_products                              ();
	apply_force_free_conditions_within_force_free_region();
	
	VSH_decomposition_for_B();
	VSH_decomposition_for_E();
	
	calculate_radial_derivatives_of_VSH_coeffs       ();
	calculate_B_1_L_and_B_1_L_dr                     ();
	calculate_radial_derivatives_of_vector_components();
	calculate_dot_products                           ();
	calculate_time_derivatives_of_vector_components  ();
	
	calculate_Cartesian_vector_components();
	
	calculate_stdev          ();
	calculate_magnetic_energy();
	count_nans               ();
	
	output_to_screen        ();
	output_to_csv_profiles  ();
	output_to_csv_history   ();
	output_to_csv_VSH_coeffs();
	output_to_csv_BCs       ();
	
	
	
	
	//----- Step at T_index = 2: Second-order Adams-Bashforth technique -----
	T_index = 2;
	T += delta_T;
	if( ( T_index >= T_index_Omega_ramp_start ) and ( T_index <= T_index_Omega_ramp_stop ) ){
		Omega += dOmega_by_dT_index;
		star_rotation_angle += Omega * delta_T * T_factor;
	}
	
	integrate_vectors_wrt_time_Adams_Bashforth_Order_2();
	
	apply_inner_boundary_conditions( true );
	apply_outer_boundary_conditions( true );
	
	ramp_up_electric_fields_for_rotation();
	
	calculate_dot_products                              ();
	apply_force_free_conditions_within_force_free_region();
	
	VSH_decomposition_for_B();
	VSH_decomposition_for_E();
	
	calculate_radial_derivatives_of_VSH_coeffs       ();
	calculate_B_1_L_and_B_1_L_dr                     ();
	calculate_radial_derivatives_of_vector_components();
	calculate_dot_products                           ();
	calculate_time_derivatives_of_vector_components  ();
	
	calculate_Cartesian_vector_components();
	
	calculate_stdev          ();
	calculate_magnetic_energy();
	count_nans               ();
	
	output_to_screen        ();
	output_to_csv_profiles  ();
	output_to_csv_history   ();
	output_to_csv_VSH_coeffs();
	output_to_csv_BCs       ();
	
	
	
	
	//----- T_index >= 3: Second-order Adams-Bashforth technique -----
	while( T_index <= n_timesteps ){
		
		T_index += 1;
		T       += delta_T;
		if( ( T_index >= T_index_Omega_ramp_start ) and ( T_index <= T_index_Omega_ramp_stop ) ){
			Omega += dOmega_by_dT_index;
			star_rotation_angle += Omega * delta_T * T_factor;
		}
		
		integrate_vectors_wrt_time_Adams_Bashforth_Order_3();
	
		apply_inner_boundary_conditions( true );
		apply_outer_boundary_conditions( true );
		
		ramp_up_electric_fields_for_rotation();
		
		calculate_dot_products                              ();
		apply_force_free_conditions_within_force_free_region();
		
		VSH_decomposition_for_B();
		VSH_decomposition_for_E();
		
		calculate_radial_derivatives_of_VSH_coeffs       ();
		calculate_B_1_L_and_B_1_L_dr                     ();
		calculate_radial_derivatives_of_vector_components();
		
		calculate_time_derivatives_of_vector_components();
		
		calculate_Cartesian_vector_components();
		
		calculate_stdev          ();
		calculate_magnetic_energy();
		count_nans               ();
		
		output_to_screen        ();
		output_to_csv_profiles  ();
		output_to_csv_history   ();
		output_to_csv_VSH_coeffs();
		output_to_csv_BCs       ();
		
	}
	
	
	
	
	//----- Output execution time and finish code -----
	output_execution_time();
	
	return 0;
}