/*
20230728_AVSH_series_test_fA_01_hardcoded.h

The function f is hardcoded into the vector component expressions, so we can use this with the normal AVSH series code in order to check if things are going wrong.

*/

double f_function( double r, double theta ){
	//return r * pow(theta,2);
	return 1.0;
}


double A_r_function( double r, double theta ){
	return f_function(r,theta) * r * cos(theta);
}


double A_theta_function( double r, double theta ){
	return f_function(r,theta) * r * sin(theta) * cos(theta);
}

double A_phi_function( double r, double theta ){
	return f_function(r,theta) * r * sin(theta) * cos(theta);
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