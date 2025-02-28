/*
AVSH_series_Michel_monopole.h
Parfrey PhD thesis s4.2.1
*/



double f_0 = 1.432;
	
	
	
double A_r_function( double r, double theta ){
	return f_0 * pow( r, -2 );
}


double A_theta_function( double r, double theta ){
	return 0;
}

double A_phi_function( double r, double theta ){
	return 0;
}








double A_r_ell_guess( double r, int ell ){
	if( ell == 0 ){
		return 2.0 * sqrt(pi) * f_0 * pow( r, -2 );
	}
	return 0;
}

double A_1_ell_guess( double r, int ell ){
	return 0;
}

double A_2_ell_guess( double r, int ell ){
	return 0;
}