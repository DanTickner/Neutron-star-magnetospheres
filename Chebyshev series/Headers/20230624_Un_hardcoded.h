/*
20230624_Un_hardcoded.h

Function for which we wish to determine a Chebyshev series.

U_n(x) is the Chebyshev polynomial of the second kind.
Choose a value of the order by changing the variable assignment below. Avoid "n" because that's an index used in the main code.

*/

int order = 2;

double f( double x ){
	
	switch( order ){
		case 0:
			return 1.0;
			break;
		case 1:
			return 2.0*x;
			break;
		case 2:
			return 4.0*x*x - 1.0;
			break;
		case 3:
			return 8.0*pow(x,3) - 4.0*x;
			break;
		
		
		default:
			return 0;
	}
	
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
	