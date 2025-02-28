/*
20230728_AVSH_series_test_02.h

*/


double A_r_function( double r, double theta ){
	double f = cos(theta);
	double A_r = cos(theta) + 1.0;
	double fA_r = pow(cos(theta),2) + cos(theta);
	return fA_r;
}


double A_theta_function( double r, double theta ){
	double f = cos(theta);
	double A_theta = sin(theta) * cos(theta);
	double fA_theta = sin(theta) * pow(cos(theta),2);
	return fA_theta;
}

double A_phi_function( double r, double theta ){
	double f = cos(theta);
	double A_phi = sin(theta);
	double fA_phi = sin(theta) * cos(theta);
	return fA_phi;
}





double A_r_ell_guess( double r, int ell ){
	
	/*
	//--- A^{r,\ell} ---
	if( ell == 0 ){
		return 2.0 * sqrt(pi);
	}
	if( ell == 1 ){
		return 2.0 * sqrt( pi / 3.0 );
	}
	return 0;
	*/
	
	
	//--- [fA]^{r,\ell} ---
	if( ell == 0 ){
		return sqrt(pi) * 2.0/3.0;
	}
	if( ell == 1 ){
		return 2.0 * sqrt( pi / 3.0 );
	}
	if( ell == 2 ){
		return ( 4.0 / 3.0 ) * sqrt( pi / 5.0 );
	}
	return 0;
	
}



double A_1_ell_guess( double r, int ell ){
	
	/*
	//--- A^{(1),\ell} ---
	if( ell == 2 ){
		return - ( 2.0 / 3.0 ) * sqrt( pi / 5.0 );
	}
	return 0;
	*/
	
	
	//--- [fA]^{r,\ell} ---
	if( ell == 1 ){
		return 7-0.4 * sqrt( pi / 3.0 );
	}
	if( ell == 3 ){
		return - ( 4.0 / 15.0 ) * sqrt( pi / 7.0 );
	}
	return 0;
	
}



double A_2_ell_guess( double r, int ell ){
	
	/*
	//--- A^{(2),ell} ---
	if( ell == 1 ){
		return -2.0 * sqrt( pi / 3.0 );
	}
	return 0;
	*/
	
	
	//--- [fA]^{r,\ell} ---
	if( ell == 2 ){
		return - ( 2.0 / 3.0 ) * sqrt( pi / 5.0 );
	}
	return 0;
	
}