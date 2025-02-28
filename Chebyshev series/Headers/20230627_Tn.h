/*
20230627_Tn.h

Function for which we wish to determine a Chebyshev series.

T_n(x) is the Chebyshev polynomial of the first  kind.
Choose a value of the order by changing the variable assignment below. Avoid "n" because that's an index used in the main code.
The singularity at x=\pm1 due to 1/sqrt(1-x^2) is handled by two if-statements and a tolerance, which can be adjusted if required.

*/

int m = 7;							// Chosen order of U_m.
double endpoint_tolerance = 1e-5;

double Un( double x, int n ){
	// Used for dfdx.
	
	if( abs( x - 1 ) < endpoint_tolerance ){
		//U_n(1) = n+1
		return (double) n+1.0;
	}
	
	if( abs( x - (-1) ) < endpoint_tolerance ){
		// U_n(-1) = (-1)^n (n+1).
		return pow(-1,n) * (double) n+1.0;
	}
	
	return sin( (double)(n+1.0) * acos(x) ) / sqrt( 1.0 -x*x );
	
}




double f( double x ){
	
	return cos( m * acos( x ) );
	
}




double dfdx( double x ){
	
	return (double) m * Un( x, m-1 );
	
}




double d2fdx2( double x ){
	
	/*
	//--- V1 ---
	std::cout << "\nd2fdx2" << std::endl;
	
	double ret = 0;
	
	for( int k=0; k<=m-2; k++ ){
		std::cout << "k=" << k << std::endl;
		
		double coeff = 0;
		
		for( int ell=k+1; ell<=m-1; ell++ ){
			if( ( (k+ell-1)%2 == 0 ) and (ell+m-1)%2 == 0 ){
				coeff += ell;
				std::cout << m <<"\t"<< k <<"\t"<< ell << std::endl;
			}
		}
		
		if( k == 0 ){
			coeff *= 0.5;
		}
		
		ret += coeff * cos( k * acos( x ) );
		
	}
	
	ret *= 4.0 * m;
	
	std::cout << std::endl;
	
	return ret;
	*/
	
	
	/*
	//--- V2 ---
	
	double ret = 0;
	
	for( int k=0; k<=m-2; k++ ){
		if( (k+m)%2 == 0 ){
		
			double coeff = 0;
			
			for( int ell=k+1; ell<=m-1; ell++ ){
				if( (ell+m)%2 == 1 ){
					coeff += ell;
				}
			}
			
			if( k == 0 ){
				coeff *= 0.5;
			}
			
			ret += coeff * cos( k * acos( x ) );
			
		}
		
	}
	
	ret *= 4.0 * m;
	
	return ret;
	*/
	
	
	//--- V3 ---
	double ret = 0;
	
	for( int k=m%2; k<=m-2; k+=2 ){
		
			double coeff = 0;
			
			for( int ell=k+1; ell<=m-1; ell+=2 ){
				coeff += ell;
			}
			
			if( k == 0 ){
				coeff *= 0.5;
			}
			
			ret += coeff * cos( k * acos( x ) );
			
		
	}
	
	ret *= 4.0 * m;
	
	return ret;

}




double a_n_guess( int n ){
	
	return (double) n == m;
	
}