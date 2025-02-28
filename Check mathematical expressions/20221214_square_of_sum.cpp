// cd OneDrive\PhD\Codes\Test mathematical expressions
// g++ 20221214_square_of_sum.cpp -o 20221214_square_of_sum

/* ( sum_{i=1}^N a_i )^2
= sum_{i=1}^N sum_{j=1}^N a_i a_j
= sum_{i=1}^N a_i^2 + 2 sum_{i=1}^N sum_{j=i+1}^N a_i a_j
= sum_{i=1}^N ( a_i^2 + 2 sum_{j=i+1}^n a_i a_j
*/

#include <iostream>
#include <stdlib.h>					// srand, rand
#include <vector>
#include "../Vector_Operations.h"	// print_vector();


//----- Global variables and function declarations -----


int main(){
	
	//----- Define variables -----
	int N = 0;					// Number of terms
	std::vector<double> A (N);
	
	
	//----- Initialise random seed -----
	// Program output will be different every time. You don't even need to compile the code again for different results - just run it again.
	// Comment-out this line to get the same output every time.
	// Numbers won't change too much. Probably a better random-seed can be chosen.
	srand( time(0) );
	
	
	//----- Populuate the lists with random numbers -----
	for( int n=0; n<N; n++ ){
		A[n] = rand() % 10;
	}
	
	
	//----- Calculate the product of the sums of the two lists -----
	double sum_A = 0;
	for( int n=0; n<N; n++ ){
		sum_A += A[n];
	}
	double square_exact = pow( sum_A, 2 );
	
	double square_guess_1 = 0;
	for( int n1=0; n1<N; n1++ ){
		for( int n2=0; n2<N; n2++ ){
			square_guess_1 += A[n1] * A[n2] ;
		}
	}
	
	double square_guess_2 = 0;
	for( int n=0; n<N; n++ ){
		square_guess_2 += pow( A[n], 2 );
	}
	for( int n1=0; n1<N; n1++ ){
		for( int n2=n1+1; n2<N; n2++ ){
			square_guess_2 += 2.0 * A[n1] * A[n2];
		}
	}
	
	double square_guess_3 = 0;
	for( int n1=0; n1<N; n1++ ){
		square_guess_3 += pow( A[n1], 2 );
		for( int n2=n1+1; n2<N; n2++ ){
			square_guess_3 += 2.0 * A[n1] * A[n2];
		}
	}
	
	
	//----- Output results -----
	print_vector( A, "A" );
	std::cout << "sum A:\t" << sum_A << std::endl;
	std::cout << "(sum A)^2\t:" << square_exact <<"\t"<< square_guess_1 <<"\t"<< square_guess_2 <<"\t"<< square_guess_3 << std::endl;
	
	return 0;
}