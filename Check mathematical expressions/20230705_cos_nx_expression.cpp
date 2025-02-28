/*
cd OneDrive\PhD\Codes\Check mathematical expressions
g++ 20230705_cos_nx_expression.cpp -o 20230705_cos_nx_expression
20230705_cos_nx_expression


cos(nx) = sum_{k=0}^{floor(n/2)} {n choose 2k} (-1)^k sin^{2k}(x) cos^{n-2k}(x)

Book 19, p57B.
*/

#include <iostream>
#include <iomanip>
#include <math.h>

double choose( int n, int k );
double cos_nx_expression( double x, int n );


int main(){
	
	//----- Define variables -----
	int n_max = 12;
	double x = 0.345;
	int w = 16;			// Fixed width.
	
	
	std::cout <<std::left<<std::setw(w)<< "n" <<std::left<<std::setw(w)<< "LHS" <<std::left<<std::setw(w)<< "RHS" <<std::left<<std::setw(w)<< "Difference" << std::endl;
	
	
	for( int n=0; n<=n_max; n++ ){
		
		double LHS = cos( n * x );
		double RHS = cos_nx_expression( x, n );
		
		std::cout <<std::left<<std::setw(w)<< n <<std::left<<std::setw(w)<< LHS <<std::left<<std::setw(w)<< RHS <<std::left<<std::setw(w)<< abs( LHS - RHS ) << std::endl;
		
	}
	
	
	return 0;
	
}




double choose( int n, int k ){
	return tgamma( n +1 ) / ( tgamma( k +1 ) * tgamma( n-k +1 ) );
}

double cos_nx_expression( double x, int n ){
	
	//--- V0 (check code compiles) ---
	//return cos( n * x );
	
	
	//--- V1 ---
	double ret = 0;
	
	for( int k=0; k<=n/2; k++ ){
		//std::cout << k << std::endl;	// Check that k goes up to floor(n/2).
		ret += choose( n, 2*k ) * pow( -1, k ) * pow( sin(x), 2*k ) * pow( cos(x), n-2*k );
		
	}
	
	return ret;
	
}