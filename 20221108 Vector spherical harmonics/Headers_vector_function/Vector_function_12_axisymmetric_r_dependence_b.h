/* Vector_function_12_axisymmetric_r_dependence_b.h
Exact expressions for the three components of a spherical vector independent of r and its spatial derivatives, and the exact VSH series coefficients.
A = ( 0.6 + e^(-r) cos(theta) ) e_r
  + ( 1.2 + r sin(theta) ) e_theta
  + ( 1.3 + r^2 P_2(cos(theta)) ) e_phi,
where P_2 is the second Legendre polynomial.
Name as B so that we can import it along with another code and perform vector calculations.
*/

double B_r_function( double r, double t, double p ){
	return 0.6 + r * cos(t);
}




double B_t_function( double r, double t, double p ){
	return 1.2 + exp(-r) * sin(t);
}




double B_p_function( double r, double t, double p ){
	return 1.3 * r*r * cos(2.0*t);
}