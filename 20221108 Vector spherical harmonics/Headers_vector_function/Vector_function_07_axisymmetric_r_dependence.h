// Vector_function_07_axisymmetric_r_dependence.h
// Exact expressions for the three components of a spherical vector independent of r and its spatial derivatives, and the exact VSH series coefficients.
// A = cos(theta) e_r + sin(theta) e_theta + P_2(cos(theta)) e_phi, where P_2 is the second Legendre polynomial.
// In this version, the functions take r as an argument even though they don't depend on r.

double A_r_function( double r, double t, double p ){
	return cos(t);
}




double A_t_function( double r, double t, double p ){
	return sin(t);
}




double A_p_function( double r, double t, double p ){
	return 0.5*(3*pow(cos(t),2)-1);
}




double div_A_function( double r, double t, double p ){
	// Verified.
	return 4.0/r * cos(t);
}




std::vector<double> curl_A_function( double r, double t, double p ){
	// Verified.
	std::vector<double> ret (3);
	ret[0] = -1.0/r * 3.0 * cos(t) * sin(t) + ( 3.0*pow(cos(t),2)-1 ) * cos(t) / ( 2.0 * r * sin(t) );
	ret[1] = -1.0/(2.0*r) * ( 3.0*pow(cos(t),2)-1 );
	ret[2] = 2.0/r * sin(t);
	return ret;
}