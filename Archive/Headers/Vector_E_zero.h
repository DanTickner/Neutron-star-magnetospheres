/*
E_zero.h
Zero electric field, to test whether the magnetic field remains the same over time.
*/




double E_r_function( double r, double t ){
	return 0;
}

double E_t_function( double r, double t ){
	return 0;
}

double E_p_function( double r, double t ){	
	return 0;
}




double E_r_L_function( double r, int L ){
	// E^{r,L}(r).
	return 0;
}

double E_1_L_function( double r, int L ){
	// E^{(1),L}(r).
	return 0;
}

double E_2_L_function( double r, int L ){
	// E^{*(2),L}(r).
	return 0;
}




double div_E_function( double r, double t ){
	// Not needed but can be useful for verification.
	return 0;
}




double alpha_function( double r, double t ){
	// Not needed but can be useful for verification.
	return 0;
}

double beta_function( double r, double t ){
	// Not needed but can be useful for verification.
	return 0;
}