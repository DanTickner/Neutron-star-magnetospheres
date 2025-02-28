/*
AVSH_series_dipole.h
*/
	
	
	
double A_r_function( double r, double theta ){
	return 2.0 * cos( theta ) * pow( r, -3 );
}


double A_theta_function( double r, double theta ){
	return sin( theta ) * pow( r, -3 );
}

double A_phi_function( double r, double theta ){
	return 0;
}








double A_r_ell_guess( double r, int ell ){
	if( ell == 1 ){
		return 4.0 * sqrt( pi / 3.0 ) * pow( r, -3 );
	}
	return 0;
}

double A_1_ell_guess( double r, int ell ){
	if( ell == 1 ){
		return -2.0 * sqrt( pi / 3.0 ) * pow( r, -3 );
	}
	return 0;
}

double A_2_ell_guess( double r, int ell ){
	return 0;
}