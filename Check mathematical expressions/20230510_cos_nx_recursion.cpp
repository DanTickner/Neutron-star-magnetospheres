/*
cd OneDrive\PhD\Codes\Check mathematical expressions
g++ 20230510_cos_nx_recursion.cpp -o 20230510_cos_nx_recursion
20230510_cos_nx_recursion


cos(nx) = 2 cos(x) cos[(n-1)x] - cos[(n-2)x]

Book 16, p60.
*/

#include <iostream>
#include <iomanip>
#include <math.h>


int main(){
	
	//----- Define variables -----
	int n_max = 12;
	double x = 0.345;
	int w = 16;			// Fixed width.
	
	
	std::cout <<std::left<<std::setw(w)<< "n" <<std::left<<std::setw(w)<< "LHS" <<std::left<<std::setw(w)<< "RHS" <<std::left<<std::setw(w)<< "Difference" << std::endl;
	
	
	for( int n=0; n<=n_max; n++ ){
		
		double LHS = cos( (double) n * x );
		double RHS = 2.0 * cos(x) * cos( (double)(n-1.0) * x ) - cos( (double)(n-2.0) * x );
		
		std::cout <<std::left<<std::setw(w)<< n <<std::left<<std::setw(w)<< LHS <<std::left<<std::setw(w)<< RHS <<std::left<<std::setw(w)<< abs( LHS - RHS ) << std::endl;
		
	}
	
	
	return 0;
	
}