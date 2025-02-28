/*
20230726_VSH_series_test_02.h

Reverse-engineer a vector with a nice VSH series expansion, to test where the code is producing inaccurate results.

*/


std::vector< std::complex<double> > A_true ( double theta, double phi ){
	
	std::vector< std::complex<double> > ret { 0, 0, 0 };

	ret = vector_sum( ret, scalar_times_vector( 1, VSH_Y  (theta,phi,1,0) ) );
	ret = vector_sum( ret, scalar_times_vector( 1, VSH_Psi(theta,phi,2,1) ) );
	ret = vector_sum( ret, scalar_times_vector( 1, VSH_Phi(theta,phi,3,2) ) );
	
	return ret;
}


std::complex<double> A_r_function( double r, double theta, double phi ){
	return A_true(theta,phi)[0];
}


std::complex<double> A_theta_function( double r, double theta, double phi ){
	return A_true(theta,phi)[1];
}

std::complex<double> A_phi_function( double r, double theta, double phi ){
	return A_true(theta,phi)[2];
}