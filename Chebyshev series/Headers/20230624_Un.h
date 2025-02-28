/*
20230624_Un.h

Function for which we wish to determine a Chebyshev series.

U_n(x) is the Chebyshev polynomial of the second kind.
Choose a value of the order by changing the variable assignment below. Avoid "n" because that's an index used in the main code.
The singularity at x=\pm1 due to 1/sqrt(1-x^2) is handled by two if-statements and a tolerance, which can be adjusted if required.

*/

int m = 5;							// Chosen order of U_m.
double endpoint_tolerance = 1e-5;

double f( double x ){
	
	if( abs( x - 1 ) < endpoint_tolerance ){
		//U_m(1) = m+1
		return (double) m+1.0;
	}
	
	if( abs( x - (-1) ) < endpoint_tolerance ){
		// U_m(-1) = (-1)^m (m+1).
		return pow(-1,m) * (double) m+1.0;
	}
	
	return sin( (double)(m+1.0) * acos(x) ) / sqrt( 1.0 -x*x );
	
}




double choose( int n, int k ){
	return tgamma( n +1 ) / ( tgamma( k +1 ) * tgamma( n-k +1 ) );
}




double a_n_guess( int n ){
	
	/*
	//--- V0 (doesn't work) ---
	double coeff = 0;
	
	for( int k=(m-n)/2; k<=m/2; k++ ){
		coeff += 2.0 * pow(-1,k) * tgamma( m-k +1 ) / ( tgamma( k +1 ) * tgamma( 0.5*(m-2*k+n) +1 ) * tgamma( 0.5*(m-2*k-n) +1 ) );
	}
	
	if( n == 0 ){
		coeff *= 0.5;
	}
	
	return coeff;
	*/
	
	
	/*
	//--- V1 ---
	double a_n = 0;
	
	for( int k=0; k<=m/2; k++ ){
		
		double term = pow(-1,k) * pow( 2, m-2*k ) * tgamma( m-k +1 ) / ( tgamma( k +1 ) * tgamma( m-2*k +1 ) );
		
		double integral = 0;
		
		if( ( m-2*k >= n ) and ( (m-2*k+n)%2 == 0 ) ){
			integral = pi * pow( 2, -(m-2*k) ) * choose( m-2*k, 0.5*(m-2*k-n) );
		}
		
		a_n += term * integral;
		
	}
	
	a_n *= 2.0 / pi;
	
	if( n == 0 ){
		a_n *= 0.5;
	}
	
	return a_n;
	*/
	
	
	/*
	//--- V2 (Let pi and 2^(m-2k) cancel, baby steps) ---
	double a_n = 0;
	
	
	for( int k=0; k<=m/2; k++ ){
		
		double term = pow(-1,k) * tgamma( m-k +1 ) / ( tgamma( k +1 ) * tgamma( m-2*k +1 ) );
		
		double integral = 0;
		
		if( ( m-2*k >= n ) and ( (m-2*k+n)%2 == 0 ) ){
			integral = choose( m-2*k, 0.5*(m-2*k-n) );
			//std::cout << n <<"\t"<< m <<"\t"<< k << std::endl;
		}
		
		a_n += term * integral;
		
	}
	
	a_n *= 2.0;
	
	if( n == 0 ){
		a_n *= 0.5;
	}
	
	return a_n;
	*/
	
	
	/*
	//--- V3 (absorb k-condition into bounds of for-loop and remove if-statement) ---
	if( (m+n)%2 == 1 ){
		return 0;
	}
	
	double a_n = 0;
	
	
	for( int k=0; k<=(m-n)/2; k++ ){
		
		double term = pow(-1,k) * tgamma( m-k +1 ) / ( tgamma( k +1 ) * tgamma( m-2*k +1 ) );
		
		double integral = choose( m-2*k, 0.5*(m-2*k-n) );
		
		a_n += term * integral;
		
	}
	
	a_n *= 2.0;
	
	if( n == 0 ){
		a_n *= 0.5;
	}
	
	return a_n;
	*/
	
	
	//--- V4 (Fix V0) ---
	if( (m+n)%2 == 1 ){
		return 0;
	}
	
	double coeff = 0;
	
	for( int k=0; k<=(m-n)/2; k++ ){
		coeff += 2.0 * pow(-1,k) * tgamma( m-k +1 ) / ( tgamma( k +1 ) * tgamma( 0.5*(m-2*k+n) +1 ) * tgamma( 0.5*(m-2*k-n) +1 ) );
	}
	
	if( n == 0 ){
		coeff *= 0.5;
	}
	
	return coeff;
	
}
	