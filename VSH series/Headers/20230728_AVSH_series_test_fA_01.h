/*
20230728_AVSH_series_test_fA_01.h

*/

double f_function( double r, double theta ){
	//return r * pow(theta,2) * log(1.0+theta);
	//return r;
	return cos(theta);
	//return sin(theta);
}



double A_r_function( double r, double theta ){
	//return r * cos(theta) + r;
	return cos(theta) + 1.0;
}


double A_theta_function( double r, double theta ){
	return sin(theta) * cos(theta);
}

double A_phi_function( double r, double theta ){
	//return r * pow(sin(theta),2);
	return sin(theta);
}