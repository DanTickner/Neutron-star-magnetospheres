// Vector_function_04.h
// Exact expressions for the three components of a spherical vector independent of r and its spatial derivatives, and the exact VSH series coefficients.
// A = Phi_1^1, where Phi_ell^m = e_r cross ( r grad Y_ell^m ) is the third vector spherical harmonic.

std::complex<double> A_r_function( double t, double p ){
	return 0;
}




std::complex<double> A_t_function( double t, double p ){
	return -sqrt(3.0/(8*pi)) * exp(std::complex<double>{0,p}) * std::complex<double>{0,-1};
}




std::complex<double> A_p_function( double t, double p ){
	return -sqrt(3.0/(8*pi)) * exp(std::complex<double>{0,p}) * cos(t);
}