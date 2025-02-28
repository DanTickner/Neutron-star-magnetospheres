/*
20230730_AVSH_series_test_fAxB_01.h

Scalar function f multiplying vector product A cross B. All functions are axisymmetric.
*/

double f_function( double r, double theta ){
	return 10 * sin(theta);
}



double A_r_function( double r, double theta ){
	return r * cos(theta);
}

double A_theta_function( double r, double theta ){
	return r * sin(2.0*theta);
}

double A_phi_function( double r, double theta ){
	return r * cos(theta);
}



double B_r_function( double r, double theta ){
	return 5 * r * pow(sin(theta),2);
}

double B_theta_function( double r, double theta ){
	return 4 * r * cos(theta);
}

double B_phi_function( double r, double theta ){
	return 3 * r * sin(theta);
}