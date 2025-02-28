/*
20230705_one_over_x_times_f_other_functions.h

Functions f and g for which we wish to determine a Chebyshev series of fg in terms of the series of f and of g.

f(x) = 1/x
g(x) = e^(1.3*x) + 3/2

Haven't bothered to find the exact Chebyshev coefficients of f, g or fg.
This code is more a test that it can be applied to unknown, not-exactly-decomposable-in-finitely-many-terms functions, and still get accurate results.
This will only work for intervals [a,b] with a>0.

*/

double f( double x ){
	return pow(x,-1);
}




double g( double x ){
	return exp( 1.3 * x ) + 1.5;
}




double fg( double x ){
	return f(x) * g(x);
}




double fg_explicit( double x ){
	return pow(x,-1) * ( exp( 1.3 * x ) + 1.5 );
}






double f_n_guess( int n ){
	
	switch( n ){
		
		case 0: return 0;
		
		default: return 0;
		
	}
	
}




double g_n_guess( int n ){
	
	switch( n ){
		
		case 0: return 0;
		
		default: return 0;
		
	}
	
}




double fg_n_guess( int n ){
	
	switch( n ){
		
		case 0: return 0;
		
		default: return 0;
		
	}
	
}