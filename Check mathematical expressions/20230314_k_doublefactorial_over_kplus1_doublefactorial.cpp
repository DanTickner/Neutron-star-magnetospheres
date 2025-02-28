/*
cd OneDrive\PhD\Codes\Test mathematical expressions
g++ 20230314_k_doublefactorial_over_kplus1_doublefactorial.cpp -o 20230314_k_doublefactorial_over_kplus1_doublefactorial
20230314_k_doublefactorial_over_kplus1_doublefactorial

Let U be the sequence with terms
u_k = k!!/(k+1)!!
Then, we have the recursion relation
u_k = k/(k+1) u_{k-2}
and so only need to generate u_0 and u_1.
This avoids terms blowing up when calculating u_k for k>19. See the below file.
20230116 Time-evolution of force-free Maxwell equations axisymmetric/CPVSH_series_E_induced_by_constant_rotation_01.cpp

Also test the following relation
k!!/(k+1)!! = (k-1)! / ( (k+1) 2^k (k/2)! if k odd
Book 14, p44.
*/

#include <iostream>
#include <iomanip>
#include <math.h>
#include <vector>


//----- Global variables and function declarations -----
int double_factorial( int n );


int main(){
	
	//----- Define variables -----
	int k_max = 30;	// Maximum value of k
	std::vector<double> sequence_explicit ( k_max+1 );
	std::vector<double> sequence_recursion( k_max+1 );
	std::vector<double> sequence_explicit2( k_max+1 );
	
	
	//----- Calculate the sequence by the recursion method -----
	sequence_recursion[0] = 1.0;
	sequence_recursion[1] = 0.5;
	
	for( int k=2; k<=k_max; k++ ){
		sequence_recursion[k] = ( (double) k ) / ( (double) k+1.0 ) * sequence_recursion[k-2];
	}
	
	
	//----- Calculate the sequence by the explicit method -----
	for( int k=0; k<=k_max; k++ ){
		sequence_explicit[k] = double_factorial(k) / ( (double) double_factorial( k+1 ) );
		
		if( k % 2 == 0 ){
			sequence_explicit2[k] = pow(2,k) * pow( tgamma( 0.5*k +1 ), 2 ) / ( (double) tgamma( k+1 +1 ) );
		} else {
			sequence_explicit2[k] = tgamma( k+1 +1 ) / ( (double) pow(2,k+1) * pow( tgamma( 0.5*(k+1) +1 ), 2 ) );
		}
	}
	
	
	//----- Output results -----
	std::cout <<std::left<<std::setw( 5)<< "k" <<std::left<<std::setw(15)<< "Recursion" <<std::left<<std::setw(15)<< "Explicit"
              <<std::left<<std::setw(15)<< "Ratio"<< std::endl;
	
	for( int k=0; k<=k_max; k++ ){
		std::cout <<std::left<<std::setw( 5)<< k <<std::left<<std::setw(15)<< sequence_recursion[k] <<std::left<<std::setw(15)<< sequence_explicit[k]
		          <<std::left<<std::setw(15)<< sequence_explicit2[k] << std::endl;
                  //<<std::left<<std::setw(15)<< sequence_recursion[k] / sequence_explicit[k] << std::endl;
	}
	
	std::cout << "Done" << std::endl;
	
	return 0;
}




int double_factorial( int n ){
	// Test mathematical expressions/Test_double_factorial.cpp
	if( n % 2 == 0 ){
		return pow( 2, 0.5*n ) * tgamma( 0.5*n + 1 );
	}
	return tgamma( n + 1 ) / ( pow( 2, 0.5*(n-1) ) * tgamma( 0.5*(n+1) ) );
}