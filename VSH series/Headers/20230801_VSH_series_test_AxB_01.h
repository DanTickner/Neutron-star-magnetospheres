/*
20230801_VSH_series_test_AxB_01.h

Scalar function f multiplying vector product A cross B.
*/

/*
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
*/


//--- Copying AVSH functions over, to check the code works ---
std::complex<double> A_r_function( double r, double theta, double phi ){
	return cos(theta) + 1.0;
}

std::complex<double> A_theta_function( double r, double theta, double phi ){
	//return sin(theta) * cos(theta);
	return sin(theta) * cos(theta) * sin(theta);	// Use this version to check that [fAxB]^{r,\ell} works, since it's zero with the main choice of vector components in this header file.
}

std::complex<double> A_phi_function( double r, double theta, double phi ){
	return sin(theta);
}



std::complex<double> B_r_function( double r, double theta, double phi ){
	return pow( cos(theta), 2 );
}

std::complex<double> B_theta_function( double r, double theta, double phi ){
	return pow( cos(theta), 2 ) * sin(theta);
}

std::complex<double> B_phi_function( double r, double theta, double phi ){
	return sin(theta) * cos(theta);
}






std::complex<double> A_r_ell_guess( double r, int ell, int m ){
	if( ( ell == 0 ) and ( m == 0 ) ){
		return 2.0 * sqrt( pi );
	}
	if( ( ell == 1 ) and ( m == 0 ) ){
		return 2.0 * sqrt( pi/3.0 );
	}
	return 0;
}

std::complex<double> A_1_ell_guess( double r, int ell, int m ){
	if( ( ell == 2 ) and ( m == 0 ) ){
		return - sqrt( pi/5.0 ) * 2.0 / 3.0;
	}
	return 0;
}

std::complex<double> A_2_ell_guess( double r, int ell, int m ){
	if( ( ell == 1 ) and ( m == 0 ) ){
		return -2.0 * sqrt( pi / 3.0 );
	}
	return 0;
}