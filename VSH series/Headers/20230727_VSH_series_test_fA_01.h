/*
20230727_VSH_series_test_fA_01.h

Scalar function f multiplying vector A.
*/

std::complex<double> f_function( double r, double theta, double phi ){
	return sin(theta);
}



std::complex<double> A_r_function( double r, double theta, double phi ){
	return r * cos(theta) * sin(phi);
}

std::complex<double> A_theta_function( double r, double theta, double phi ){
	return r * sin(2.0*theta)* pow(cos(phi),2);
}

std::complex<double> A_phi_function( double r, double theta, double phi ){
	return r * cos(theta) * exp(-phi);
}