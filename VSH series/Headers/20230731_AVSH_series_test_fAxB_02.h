/*
20230731_AVSH_series_test_fAxB_02.h

Scalar function f multiplying vector product A cross B. All functions are axisymmetric.
*/

double f_function( double r, double theta ){
	return cos(theta);
}



double A_r_function( double r, double theta ){
	return cos(theta) + 1.0;
}

double A_theta_function( double r, double theta ){
	return sin(theta) * cos(theta);
	//return sin(theta) * cos(theta) * sin(theta);	// Use this version to check that [fAxB]^{r,\ell} works, since it's zero with the main choice of vector components in this header file.
}

double A_phi_function( double r, double theta ){
	return sin(theta);
}



double B_r_function( double r, double theta ){
	return pow( cos(theta), 2 );
}

double B_theta_function( double r, double theta ){
	return pow( cos(theta), 2 ) * sin(theta);
}

double B_phi_function( double r, double theta ){
	return sin(theta) * cos(theta);
}




double fAxB_r_function( double r, double theta ){
	return 0;
}

double fAxB_theta_function( double r, double theta ){
	return - pow( cos(theta), 2 ) * sin(theta);
}

double fAxB_phi_function( double r, double theta ){
	return pow( cos(theta), 3 ) * sin(theta);
}