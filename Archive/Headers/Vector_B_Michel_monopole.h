/* 
Vector_B_Michel_monopole.h

Parfrey PhD thesis section 4.2.1
Michel F. C., 1973, ApJL, 180, https://ui.adsabs.harvard.edu/abs/1973ApJ...180L.133M/abstract

B = f_0/r^2 e_r with f_0 an arbitrary scalar.
*/


double f_0 = 1;

double B_r_function( double r, double t ){
	return f_0 * pow( r, -2 );
}

double B_t_function( double r, double t ){
	return 0;
}

double B_p_function( double r, double t ){
	return 0;
}




double B_r_L_function( double r, int L ){
	// B^{r,L}(r).
	if( L == 0 ){
		return 2.0 * sqrt(pi) * f_0 * pow( r, -2 );
	}
	return 0;
}

double B_1_L_function( double r, int L ){
	// B^{(1),L}(r).
	return 0;
}

double B_2_L_function( double r, int L ){
	// B^{*(2),L}(r).
	return 0;
}