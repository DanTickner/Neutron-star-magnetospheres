/*
B_r_1.h

VSH coefficients for magnetic field of a dipole

B^{r,1}(r) = 4 sqrt(pi/3) 1/r^3

Can only get numerical results for the Chebyshev coefficients of this function, but useful to compare them to what my evolutionary code outputs.

*/

double f( double x ){
	return 4.0 * sqrt(pi/3.0) * pow( x, -3 );
}




double dfdx( double x ){
	return -12.0 * sqrt(pi/3.0) * pow( x, -4 );
}




double d2fdx2( double x ){
	return 48.0 * sqrt(pi/3.0) * pow( x, -5 );
}



double a_n_guess( int n ){
	return 123;
}