/*
E_induced_by_constant_rotation_dipole.h
Electric field induced by a rotating magnetic dipole with constant angular velocity in the phi-direction.
Book 16 p58.
*/




double E_r_function( double r, double t ){
	return Omega * pow( r, -2 ) * pow( sin(t), 2 );
}

double E_t_function( double r, double t ){
	return Omega * -2.0 * pow( r, -2 ) * sin(t) * cos(t);
}

double E_p_function( double r, double t ){	
	return 0;
}




double E_r_L_function( double r, int L ){
	// E^{r,L}(r).
	if( L == 0 ){
		return (4.0/3.0) * sqrt( pi ) * Omega * pow( r, -2 );
	}
	if( L == 2 ){
		return -(4.0/3.0) * sqrt( 0.2*pi ) * Omega * pow( r, -2 );
	}
	return 0;
}

double E_1_L_function( double r, int L ){
	// E^{(1),L}(r).
	if( L == 2 ){
		return (4.0/3.0) * sqrt( 0.2*pi ) * Omega * pow( r, -2 );
	}
	return 0;
}

double E_2_L_function( double r, int L ){
	// E^{*(2),L}(r).
	return 0;
}




double div_E_function( double r, double t ){
	// Not needed but can be useful for verification.
	return Omega * 2.0 * pow(r,-3) * ( 1.0 - 3.0*pow(cos(t),2) );
}




double alpha_function( double r, double t ){
	// Not needed but can be useful for verification.
	return 0;
}

double beta_function( double r, double t ){
	// Not needed but can be useful for verification.
	return 2.0 * Omega * pow(r,3) * ( 3.0*pow(cos(t),2)-1.0 ) / ( 3.0*pow(cos(t),2) + 1.0 );
}