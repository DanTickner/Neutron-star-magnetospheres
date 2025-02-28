// Vector_function_10_dipole_log.h
// Exact expressions for the magnetic field of an electric dipole.
// B = mu_0 m / ( 4pi ) * ( 2cos(theta)/r^3 e_r + sin(theta)/r^3 e_theta ),
// mu_0 is the permeability of free space and m is the magnitude of the magnetic dipole moment. Set both to one.
// For ease of implementation, just set mu_0 and m within each function definition by hardcoding it.
// Replace the 1/r^3 with ln(2-x), which is a well-behaved monotonically decreasing function.
// Do this because you can't create a Chebyshev series for f(x)=1/x^n, due to the appearance of singularities and non-existent limits in the integrand.

double A_r_function( double r, double t, double p ){
	double mu_0 = 1;
	double m    = 1;
	return mu_0*m * pow(4*pi,-1) * log(2.0-r) * 2.0 * cos(t);
}




double A_t_function( double r, double t, double p ){
	double mu_0 = 1;
	double m    = 1;
	return mu_0*m * pow(4*pi,-1) * log(2.0-r) * sin(t);
}




double A_p_function( double r, double t, double p ){
	return 0;
}




double div_A_function( double r, double t, double p ){
	// NOT Verified!
	return 0;
}




std::vector<double> curl_A_function( double r, double t, double p ){
	// See probationary review. Typo in writeup so incomplete.
	std::vector<double> ret (3);
	double mu_0 = 1;
	double m    = 1;
	
	ret[0] = 0;
	ret[1] = 0;
	ret[2] = 123;
	return ret;
}