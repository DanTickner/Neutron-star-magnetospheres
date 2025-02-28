/*
cd OneDrive\PhD\Codes\Check mathematical expressions
g++ 20230624_Chebyshev_polynomial_second_kind_in_terms_of_first_kind.cpp -o 20230624_Chebyshev_polynomial_second_kind_in_terms_of_first_kind
20230624_Chebyshev_polynomial_second_kind_in_terms_of_first_kind


U_n(x) = 2 * sum_{k=n(mod2)//k+n even}^n h_k T_k

h_k = 0.5 if k=0 else 1

Book 19, p8
*/

#include <iostream>
#include <iomanip>
#include <math.h>

double Un( double x, int n );
double Tn( double x, int n );
double Un_hardcoded( double x, int n );
double Un_in_terms_of_Tn( double x, int n );


int main(){
	
	//----- Define variables -----
	int n_max = 12;
	double x = 0.345;
	int w = 16;			// Fixed width.
	
	
	std::cout <<std::left<<std::setw(w)<< "n" <<std::left<<std::setw(w)<< "LHS" <<std::left<<std::setw(w)<< "RHS" <<std::left<<std::setw(w)<< "Difference" << std::endl;
	
	
	for( int n=0; n<=n_max; n++ ){
		
		double LHS = Un( x, n );
		double RHS = Un_in_terms_of_Tn( x, n );
		
		std::cout <<std::left<<std::setw(w)<< n <<std::left<<std::setw(w)<< LHS <<std::left<<std::setw(w)<< RHS <<std::left<<std::setw(w)<< abs( LHS - RHS ) << std::endl;
		
	}
	
	
	return 0;
	
}




double Un( double x, int n ){
	
	double endpoint_tolerance = 1e-5;	// Freely adjustable.
	
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

double Tn( double x, int n ){
	
	return cos( n * acos( x ) );
	
}




double Un_hardcoded( double x, int n ){
	
	switch( n ){
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




double Un_in_terms_of_Tn( double x, int n ){
	
	
	/*
	//--- V1 ---
	double Un = 0;
	
	for( int k=n%2; k<=n; k+=2 ){
		double term = 2.0 * Tn( x, k );
		if( k == 0 ){
			term *= 0.5;
		}
		Un += term;
	}
	*/
	
	
	//--- V2 ---
	//This is the best version. V3 is just to check something.
	/*
	double Un = 0;
	
	for( int k=n%2; k<=n; k+=2 ){
		Un += Tn( x, k );
	}
	
	Un *= 2.0;
	
	if( n%2 == 0 ){
		Un -= 1.0;
	}
	
	return Un;
	*/
	
	
	/*
	//--- V3 ---
	// Write as infinite sum multiplied by a condition g(k). This will help when changing indices of double sums. Book 19, p15.
	
	double Un = 0;
	
	int infinity = 3*n;	// "infinity" for upper limit
	
	for( int k=0; k<=infinity; k++ ){
		
		//--- Evaluate g(k) ---
		// If g(k)=0, skip to next k. This just saves things like list indices overflowing, and is equivalent to multiplying that term by 0.
		bool condition_1 = k >= ( n % 2 );
		bool condition_2 = k <= n;
		bool condition_3 = ( k + n ) % 2 == 0;
		
		bool g_k = condition_1 and condition_2 and condition_3;
		
		//
		if ( ! g_k ){
			continue;
		}
		//
		
		double h_k = 1.0;
		if( k == 0 ){
			h_k = 0.5;
		}
		double T_k = cos( k * acos( x ) );
		
		Un += 2.0 * h_k * T_k * g_k;
		
		
	}
	
	
	return Un;
	*/
	
	
	//--- V4 (sum to infinity with condition that k is even, not needing mod(2) function in the summation but as a condition instead) ---
	double Un = 0;
	
	for( int k=0; k<=n; k++ ){
		if( (k+n)%2 == 0 ){
			double term = 2.0 * Tn( x, k );
			if( k == 0 ){
				term *= 0.5;
			}
			Un += term;
		}
	}
	
	return Un;
}