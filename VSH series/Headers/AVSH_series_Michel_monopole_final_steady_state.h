/*
AVSH_series_Michel_monopole_final_steady_state.h
Parfrey PhD thesis s4.2.1
*/



double f_0   = 1.432;
double Omega = 0.33;
	
	
	
double A_r_function( double r, double theta ){
	return 0;
}


double A_theta_function( double r, double theta ){
	return 0;
}

double A_phi_function( double r, double theta ){
	return f_0 * Omega * sin(theta) * pow( r, -1 );
}








double A_r_ell_guess( double r, int ell ){
	return 0;
}

double A_1_ell_guess( double r, int ell ){
	return 0;
}

double A_2_ell_guess( double r, int ell ){
	if( ell == 1 ){
		return - 2.0 * sqrt( pi / 3.0 ) * f_0 * Omega * pow( r, -1 );
	}
	return 0;
}