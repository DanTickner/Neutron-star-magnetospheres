/*
20230727_VSH_series_test_fAxB_01.h

Scalar function f multiplying vector product A cross B.
*/

std::complex<double> f_function( double r, double theta, double phi ){
	return 10 * sin(theta);
}



std::complex<double> A_r_function( double r, double theta, double phi ){
	return r * cos(theta) * sin(phi);
}

std::complex<double> A_theta_function( double r, double theta, double phi ){
	return r * sin(2.0*theta)* pow(cos(phi),2);
}

std::complex<double> A_phi_function( double r, double theta, double phi ){
	return r * cos(theta) * exp(std::complex<double>{0,-2*phi});
}



std::complex<double> B_r_function( double r, double theta, double phi ){
	return 5 * r * pow(sin(theta),2) * cos(phi);
}

std::complex<double> B_theta_function( double r, double theta, double phi ){
	return 4 * r * cos(theta)* pow(cos(phi),3);
}

std::complex<double> B_phi_function( double r, double theta, double phi ){
	return 3 * r * sin(theta) * exp(std::complex<double>{0,phi});
}