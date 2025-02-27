/* 
Vector_B_static_dipole_A_equals_minus2Psi.h
Exact expressions for the magnetic field of an electric dipole, with an electric field that ramps up linearly from zero to some maximum value.
B = mu_0 m / ( 4pi ) ( 2cos(theta)/r^3 e_r + sin(theta)/r^3 e_theta  - 2 Psi e_r),
where Psi = 1/r sin^2(theta).
We are working in code units, so the leading factors are set to 1, and r is normalised by the stellar radius.
*/




double B_r_function( double r, double t ){
	return pow( r, -3 ) * 2.0 * cos( t );
}

double B_t_function( double r, double t ){
	return pow( r, -3 ) * sin( t );
}

double B_p_function( double r, double t ){
	double Psi = pow( sin(t), 2 ) / r;
	return -2.0 * Psi;
}




double B_r_ell_function( double r, int ell ){
	// B^{r,L}(r).
	if( ell == 1 ){
		return 4.0 * sqrt( pi / 3.0 ) * pow( r, -3 );
	}
	return 0;
}

double B_1_ell_function( double r, int ell ){
	// B^{(1),L}(r).
	return -0.5 * B_r_ell_function( r, ell );
}

double B_2_ell_function( double r, int ell ){
	// B^{*(2),L}(r).
	return B_r_ell_function( r, ell );
}