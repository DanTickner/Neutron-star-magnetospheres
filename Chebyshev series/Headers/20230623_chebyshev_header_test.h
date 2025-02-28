/*
20230623_chebyshev_header_test.h

Function for which we wish to determine a Chebyshev series.

f(x) = 3x^4 + 4x^3 - 2x^2 - x + 1

*/

double f( double x ){
	return 3.0*pow(x,4) + 4.0*pow(x,3) - 2.0*pow(x,2) - x + 1.0;
}




double dfdx( double x ){
	return 12.0*pow(x,3) + 12.0*pow(x,2) - 4.0*x - 1.0;
}




double d2fdx2( double x ){
	return 36.0*pow(x,2) + 24.0*x - 4.0;
}




double a_n_guess( int n ){
	if( n == 0 ){
		return 0.123;
	}
	if( n == 1 ){
		return 0.234;
	}
	if( n == 2 ){
		return 0.432;
	}
	if( n == 3 ){
		return 0.543;
	}
	if( n == 4 ){
		return 0.210;
	}
	return 0;
}