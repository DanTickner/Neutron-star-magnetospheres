/*
20230728_AVSH_series_test_01.h

*/


double A_r_function( double r, double theta ){
	return r * cos(theta);
}


double A_theta_function( double r, double theta ){
	return r * sin(theta) * cos(theta);
}

double A_phi_function( double r, double theta ){
	//return r * pow(sin(theta),2);
	return r * sin(theta) * cos(theta);
}



double A_r_ell_guess( double r, int ell ){
	if( ell == 1 ){
		return 2.0 * sqrt( pi/3.0 ) * r;
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