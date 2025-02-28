/* Vector_function_09_axisymmetric_r_dependence.h
Exact expressions for the three components of a spherical vector independent of r and its spatial derivatives, and the exact VSH series coefficients.
A = e^(-r) cos(theta) e_r
  + r sin(theta) e_theta
  + r^2 P_2(cos(theta)) e_phi,
where P_2 is the second Legendre polynomial.
*/

double A_r_function( double r, double t, double p ){
	return r * cos(t);
}




double A_t_function( double r, double t, double p ){
	return exp(-r) * sin(t);
}




double A_p_function( double r, double t, double p ){
	return r*r * cos(2.0*t);
}