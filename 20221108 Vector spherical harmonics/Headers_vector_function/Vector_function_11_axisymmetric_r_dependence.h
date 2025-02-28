/* Vector_function_11_axisymmetric_r_dependence.h
Exact expressions for the three components of a spherical vector independent of r and its spatial derivatives, and the exact VSH series coefficients.
A = ( 1 + e^(-r) cos(theta) ) e_r
  + ( 0.9 + r sin(theta) ) e_theta
  + ( 0.7 + r^2 P_2(cos(theta)) ) e_phi,
where P_2 is the second Legendre polynomial.
*/

double A_r_function( double r, double t, double p ){
	return 1.0 + exp(-r) * cos(t);
}




double A_t_function( double r, double t, double p ){
	return 0.9 + r * sin(t);
}




double A_p_function( double r, double t, double p ){
	return 0.7 + r*r * 0.5*(3*pow(cos(t),2)-1);
}




double div_A_function( double r, double t, double p ){
	// Not yet verified.
	return ( -1.0 + 2.0/r ) * exp(-r) * cos(t) + 2.0*cos(t);
}




std::vector<double> curl_A_function( double r, double t, double p ){
	// Verified.
	std::vector<double> ret (3);
	ret[0] = -3.0 * r * cos(t)*sin(t) + 0.5 * r * cos(t) * pow(sin(t),-1) * ( 3.0*pow(cos(t),2)-1 );
	ret[1] = -1.5 * r * ( 3.0*pow(cos(t),2)-1 );
	ret[2] = 2.0 * sin(t) + pow(r,-1) * exp(-r) * sin(t);
	return ret;
}