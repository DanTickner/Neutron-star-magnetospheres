// cd OneDrive\PhD\Codes\Test mathematical expressions
// g++ 20220829_sum_for_integral_1.cpp -o 20220829_sum_for_integral_1 -std=c++11

// -x^{2n+2} - (2n+2)! sum_{k=0}^n (-1)^k/(2k)! x^{2k} ?=? (2n+2)! sum_{k=0}^{n+1} (-1)^k/(2k)! x^{2k}

#include <iostream>
#include <math.h>


//----- Function declarations -----
double LHS( int n, double x );
double RHS( int n, double x );


int main(){
	
	//----- Define variables -----
	int n = 3;
	double x = 0.4;
	
	std::cout << LHS( n, x ) << "\t" << RHS( n, x ) << std::endl;
	
	return 0;
}


double LHS( int n, double x ){
	
	double ret = 0;
	
	for( int k=0; k<=x; k++ ){
		ret += pow(x,2*k) * pow(-1,k) / tgamma(2*k+1);
	}
	
	return - pow(x,2*n+2) - tgamma(2*n+2+1) * ret;
}


double RHS( int n, double x ){
	
	double ret = 0;
	
	for( int k=0; k<=x+1; k++ ){
		ret += pow(x,2*k) * pow(-1,k) / tgamma(2*k+1);
	}
	
	return tgamma(2*n+2+1) * ret;
}