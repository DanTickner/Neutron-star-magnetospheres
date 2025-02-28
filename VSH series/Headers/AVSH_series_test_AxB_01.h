/*
AVSH_series_test_AxB_01.h

Scalar function f multiplying vector product A cross B. All functions are axisymmetric.

Based on 20230731_AVSH_series_test_fAxB_02.h.
*/



double A_r_function( double r, double theta ){
	return cos(theta) + 1.0;
}

double A_theta_function( double r, double theta ){
	return sin(theta) * cos(theta);
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




double AxB_r_function( double r, double theta ){
	return 0;
}

double AxB_theta_function( double r, double theta ){
	return - cos(theta) * sin(theta);
}

double AxB_phi_function( double r, double theta ){
	return pow( cos(theta), 2 ) * sin(theta);
}