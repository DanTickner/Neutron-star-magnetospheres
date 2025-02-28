/*
20230628_one_over_x.h

Function for which we wish to determine a Chebyshev series.

f(x) = 1/x

*/

double f( double x ){
	return pow( x, -1 );
}




double dfdx( double x ){
	return - pow( x, -2 );
}




double d2fdx2( double x ){
	return 2.0 * pow( x, -3 );
}



double delta( int a, int b ){
	return ( a == b );
}

double h_n( int n ){
	return 1.0 - 0.5 * delta( n, 0 );
}

double a_n_guess( int n ){
	//return 2.0 * pow( -1, 0.5*(n-1) ) * h_n(n) * delta( n%2, 1 );
	return 2 * ( delta( n%4, 1 ) - delta( n%4, 3 ) );
}