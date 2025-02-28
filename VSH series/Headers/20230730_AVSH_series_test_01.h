/*
20230730_AVSH_series_test_01.h

// The vector B used in my specific example for fAxB.

*/


double A_r_function( double r, double theta ){
	return pow( cos(theta), 2 );
}


double A_theta_function( double r, double theta ){
	return pow( sin(theta), 2 );
}

double A_phi_function( double r, double theta ){
	return cos(theta);
}



double A_r_ell_guess( double r, int ell ){
	if( ell == 0 ){
		return ( 2.0 / 3.0 ) * sqrt(pi);
	}
	if( ell == 2 ){
		return ( 4.0 / 3.0 ) * sqrt( pi / 5.0 );
	}
	return 0;
}

double A_1_ell_guess( double r, int ell ){
	if( ell == 2 ){
		return r * - sqrt( pi/5.0 ) * 2.0 / 3.0;
	}
	return 0;
}

double A_2_ell_guess( double r, int ell ){
	if( ell == 0 ){
		return r * 2.0/3.0;
	}
	if( ell == 2 ){
		return r * -2.0/3.0;
	}
	return 0;
}