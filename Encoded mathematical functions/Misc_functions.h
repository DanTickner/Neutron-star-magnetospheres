/* Misc_functions.h

Miscellaneous mathematical functions which I have defined.

*/


double choose( int n, int k ){
	return tgamma( n +1 ) / ( tgamma( k +1 ) * tgamma( n-k +1 ) );
}




double double_factorial( int n ){
	// Test mathematical expressions/Test_double_factorial.cpp
	// case n=-1 added for this code. Probably a nice general rule to use for -n but this is good enough for our purposes.
	// This function not needed in final version of integral_guess_definite().
	if( n == -1 ){
		return 1;
	}
	if( n % 2 == 0 ){
		return pow( 2, 0.5*n ) * tgamma( 0.5*n + 1 );
	}
	return tgamma( n + 1 ) / ( pow( 2, 0.5*(n-1) ) * tgamma( 0.5*(n+1) ) );
}




int delta( int a, int b ){
	// Kronecker delta symbol.
	if( a == b ){
		return 1;
	}
	return 0;
}