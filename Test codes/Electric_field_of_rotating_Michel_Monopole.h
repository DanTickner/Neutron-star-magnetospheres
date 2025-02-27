/* 
Electric_field_of_rotating_Michel_Monopole.h

Parfrey S4.2.1

B = f_0/r^2 e_r

Results in E = 
Exact expressions for the magnetic field of an electric dipole,
B = mu_0 m / ( 4pi ) ( 2cos(theta)/r^3 e_r + sin(theta)/r^3 e_theta ).
We are working in code units, so the leading factors are set to 1, and r is normalised by the stellar radius.
The leading factor is unity because, no matter the magnetic dipole moment or stellar radius chosen, we have normalised by those values.

V03: Allow for the header file ../Initial_conditions_01.h to have already been included.
     No longer need a constant value for some small radius because we will always set r_min >= 1.
     Name change; previous version was Expression_B_dipole_with_min_r_axisymmetric_02.h
	 div(B) and curl(B) functions no longer needed. See Archive for previous version containing them.
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
	if( ell == 0 ){
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