// Vector_function_03.h
// Exact expressions for the three components of a spherical vector independent of r and its spatial derivatives, and the exact VSH series coefficients.
// A = Psi_1^1, where Psi_ell^m = r grad Y_ell^m is the second vector spherical harmonic.

std::complex<double> A_r_function( double t, double p ){
	return 0;
}




std::complex<double> A_t_function( double t, double p ){
	return -sqrt(3.0/(8*pi)) * exp(std::complex<double>{0,p}) * cos(t);
}




std::complex<double> A_p_function( double t, double p ){
	return -sqrt(3.0/(8*pi)) * exp(std::complex<double>{0,p}) * std::complex<double>{0,1};
}