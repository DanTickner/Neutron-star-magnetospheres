/* Vector_function_06.h
Exact expressions for the three components of a spherical vector independent of r and its spatial derivatives, and the exact VSH series coefficients.
B = [ sin(theta) + 1 + i cos(theta) ] e_r
  + [ cos(theta) cos(phi) + i ( theta + cos(phi) ) ] e_theta
  + [ sin(theta) * P_2(cos(phi)) + i cos(theta) cos(phi) ] e_phi,
where P_2 is the second Legendre polynomial.
Use different name to _05 so that both header files can be called at once.
*/

std::complex<double> B_r_function( double t, double p ){
	double re = sin(t) + 1;
	double im = cos(t);
	return std::complex<double>{ re, im };
}




std::complex<double> B_t_function( double t, double p ){
	double re = cos(t) * cos(p);
	double im = t + cos(p);
	return std::complex<double> { re, im };
}




std::complex<double> B_p_function( double t, double p ){
	double re = sin(t) * 0.5*(3*pow(cos(p),2)-1);
	double im = cos(t)*cos(p);
	return std::complex<double> { re, im };
}