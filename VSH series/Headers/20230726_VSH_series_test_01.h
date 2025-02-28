/*
20230726_VSH_series_test_01.h

*/


std::complex<double> A_r_function( double r, double theta, double phi ){
	return r * cos(theta) * sin(phi);
}


std::complex<double> A_theta_function( double r, double theta, double phi ){
	return r * exp(-sin(theta)) * pow(cos(phi),4);
}

std::complex<double> A_phi_function( double r, double theta, double phi ){
	return r * pow(cos(theta),5) * exp(-phi);
}