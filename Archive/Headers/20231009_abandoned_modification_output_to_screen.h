/*
ABANDONED ATTEMPT
Estimate the time remaining depending whether we are before the ramp onset or after it.
This is because we have different ell_max for either case.
It's proving too overcomplicated and in the future I don't want too much time during the initial phase anyway, so abandon and focus on things that make the code actually work.
There is a compilation error due to unmatched parentheses.



*/




void output_to_screen( int T_index, int n ){
	
	if( T_index == T_index_Omega_ramp_start ){
		double time_now_seconds = std::chrono::high_resolution_clock::now().time_since_epoch().count() * 1e-9;
		time_taken_before_ramp  = time_now_seconds - time_start_seconds;
	}
	
	if( ( T_index % cout_freq_T == 0 ) or ( T_index == n_timesteps ) ){
	
		double percent_done = round( (double)1000*T_index/n_timesteps )/10; // This format makes it exactly one decimal place.
		
		std::cout <<std::left<<std::setw(w)<< T_index <<std::left<<std::setw(w)<< T <<"|\t"
				  <<std::left<<std::setw(w)<< B_r[n][cout_i][cout_j] <<std::left<<std::setw(w)<< B_t[n][cout_i][cout_j] <<std::left<<std::setw(w)<< B_p[n][cout_i][cout_j] <<"|\t"
				  <<std::left<<std::setw(w)<< E_r[n][cout_i][cout_j] <<std::left<<std::setw(w)<< E_t[n][cout_i][cout_j] <<std::left<<std::setw(w)<< E_p[n][cout_i][cout_j] <<"|\t"
				  <<std::left<<std::setw(w)<< E_dot_B[cout_i][cout_j] <<std::left<<std::setw(w)<< Bsq_minus_Esq[cout_i][cout_j]
				  <<std::left<<std::setw(w)<< percentage_E_points_points_changed_for_second_force_free_condition <<"|\t"
				  <<std::left<<std::setw(w)<< nans_tot <<std::left<<std::setw(w)<< percent_done;
		
		if( T_index == 0 ){
			std::cout << std::endl;
		} else {
			double time_now_seconds = std::chrono::high_resolution_clock::now().time_since_epoch().count() * 1e-9;
			double runtime = time_now_seconds - time_start_seconds;
			
			double time_left = 0;
			if( T_index < T_index_Omega_ramp_start ){
				double time_left_before_ramp = ( (double) ( T_index_Omega_ramp_start + 1.0 ) / T_index - 1.0 ) * runtime;
				double time_left_after_ramp  = (double) ( T_index_Omega_ramp_stop - T_index_Omega_ramp_start ) / ( ( T_index_Omega_ramp_start + 1.0 ) * runtime * ( ell_max_rotation_on + 1.0 ) / ( (double) ell_max_rotation_off + 1.0 );
				std::cout << "before\t" << time_left_before_ramp <<"\tafter"<< time_left_after_ramp << std::endl;
				time_left = ( time_left_before_ramp * ( ell_max_rotation_off + 1.0 ) + time_left_after_ramp * ( ell_max_rotation_on + 1.0 ) ) / ( (double) ell_max_rotation_on + 1.0 );
			}
			else{
				time_left = ( (double) n_timesteps / T_index - 1.0 ) * ( runtime - time_taken_before_ramp );
			}
			int time_left_h = time_left / 3600;
			int time_left_m = time_left / 60 - time_left_h * 60;
			int time_left_s = time_left - time_left_h * 3600 - time_left_m * 60;
			
			if( time_left_m < 10 ){
				std::cout << time_left_h << ":0" << time_left_m;
			} else {
				std::cout << time_left_h << ":"  << time_left_m;
			}
			if( time_left_s < 10 ){
				std::cout << ":0" << time_left_s << std::endl;
			} else {
				std::cout << ":"  << time_left_s << std::endl;
			}
		}
	}
}