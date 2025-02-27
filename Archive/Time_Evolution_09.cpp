/*
cd OneDrive\PhD\Codes\20230810 Time evolution with updated equations
g++ Time_Evolution_09.cpp -o Time_Evolution_09 --std=c++11
Time_Evolution_09
Evolve the magnetic and electric fields over time, by finite differencing the radial direction and using a VSH expansion in the theta direction.
	 

V09: Now fully committed to starting from E=0 and ramping up linearly.
     Expression for dE/dT during the ramping-up phase of rotation depends on the chosen magnetic field, so has been relegated to the header file.
     Since we always start from E=0, we no longer need a header file for the electric field.
	 Generalise inner boundary conditions for any NS surface magnetic field, by using double B_r_function() etc from header file.
	 Relegate ALL user-chosen variables to the initial_conditions_xx.h header file, for consistency. Having some here may have saved time moving between files, but the choice of what to place where was always arbitrary and would be confusing to other users.
	 Add third CSV file to store calculated VSH series coeffs. May be useful to see their evolution over time, and analyse whether a certain ell_max is really necessary.
	 Add function to count the number of nans in the B and E components. Later, can add if-statement to stop code if this exceeds zero or a set threshold.
	 void calculate_n_gridpoints_within_force_free_region() has been moved inside void apply_force_free_conditions_within_force_free_region(), so don't call it at the start of the code.
	 Added log file output.

*/

#include <iostream>
#include <iomanip>
#include <vector>
#include <math.h>
#include <fstream>
#include <chrono>

#include "Headers/Initial_conditions_02.h"
#include "../Vector_operations.h"
#include "Headers/Vector_B_rotating_dipole_ramp.h"
#include "Headers/Header_Time_Evolution_12.h"



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
	std::cout << "Omega_max           :\t" <<std::left<<std::setw(16)<< Omega_max<< std::endl;
	std::cout << "R_LC_max            :\t" <<std::left<<std::setw(16)<< 1.0 / Omega_max    << std::endl;
	std::cout << "r_min, r_max        :\t" <<std::left<<std::setw(16)<< r_min   <<std::left<<std::setw(16)<< r_max              << std::endl;
	std::cout << "Timestep            :\t" <<std::left<<std::setw(16)<< delta_T <<std::left<<std::setw(16)<< delta_T * T_factor << std::endl;
	std::cout << "CFL max timestep    :\t" <<std::left<<std::setw(16)<< delta_r <<std::left<<std::setw(16)<< delta_r * T_factor << std::endl;
	std::cout << "Length of simulation:\t" <<std::left<<std::setw(16)<< T_max << T_max * T_factor <<"\tsteps: " << n_timesteps << std::endl;
	std::cout << "Est. CSV size (MB)  :\t" << estimate_csv_profiles_filesize() << std::endl;
	std::cout << "ell_max             :\t" << ell_max << std::endl;
	std::cout << "Printing values at ( r[" << cout_i <<"], theta[" << cout_j <<"] ) = ( " << r[cout_i] <<", " << t[cout_j] << " ) = ( " << r[cout_i] <<", " << t[cout_j]/pi << " pi )." << std::endl;
	std::cout << "dOmega_by_dT_index  :\t" << dOmega_by_dT_index <<"\tdOmega_by_dT:\t"<< dOmega_by_dT << std::endl;
	std::cout << "Ramp start          :\t" << T_index_Omega_ramp_start <<"\t"<< T_index_Omega_ramp_start * delta_T << std::endl;
	std::cout << "Ramp stop           :\t" << T_index_Omega_ramp_stop  <<"\t"<< T_index_Omega_ramp_stop  * delta_T << std::endl;
	
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
	
	apply_initial_field_value_and_VSH_decomposition_for_B();
	
	VSH_decomposition( B_r, B_t, B_p, B_r_L, B_1_L, B_2_L );	// 20230922 Trying during troubleshooting why B_1_L coeffs inaccurate.
	
	apply_force_free_conditions_within_force_free_region();

	calculate_stdev();
	calculate_force_free_conditions();
	
	calculate_radial_derivatives_of_VSH_coeffs();
	calculate_B_1_L_and_B_1_L_dr_initial();
	
	calculate_time_derivatives( 0 );
	
	calculate_Cartesian_vector_components();
	evaluate_radial_derivatives();
	
	calculate_magnetic_energy();
	
	count_nans();
	output_to_screen        ( 0 );
	output_to_csv_profiles  ( 0 );
	output_to_csv_history   ( 0 );
	output_to_csv_VSH_coeffs( 0 );
	
	integrate_vectors_wrt_time();
	apply_inner_boundary_conditions( 0 );
	apply_outer_boundary_conditions();
	
	
	
	
	//----- Evolve to future timesteps -----
	for( int T_index=1; T_index<=n_timesteps; T_index++ ){
	
		T += delta_T;
		if( ( T_index >= T_index_Omega_ramp_start ) and ( T_index <= T_index_Omega_ramp_stop ) ){
			Omega += dOmega_by_dT_index;
		}
		star_rotation_angle += Omega * delta_T * T_factor;
		
		apply_force_free_conditions_within_force_free_region();
		
		VSH_decomposition( B_r, B_t, B_p, B_r_L, B_1_L, B_2_L );
		VSH_decomposition( E_r, E_t, E_p, E_r_L, E_1_L, E_2_L );
		calculate_stdev();
		calculate_force_free_conditions();
		
		calculate_radial_derivatives_of_VSH_coeffs();
		calculate_B_1_L_and_B_1_L_dr();
		
		calculate_time_derivatives( T_index );
		
		calculate_Cartesian_vector_components();
		evaluate_radial_derivatives();
		
		calculate_magnetic_energy();
		
		count_nans();
		output_to_screen        ( T_index );
		output_to_csv_profiles  ( T_index );
		output_to_csv_history   ( T_index );
		output_to_csv_VSH_coeffs( T_index );
		
		integrate_vectors_wrt_time();
		apply_inner_boundary_conditions( T_index );
		apply_outer_boundary_conditions();
		
	}
	
	
	
	
	//----- Output execution time and finish code -----
	output_execution_time();
	
	return 0;
}