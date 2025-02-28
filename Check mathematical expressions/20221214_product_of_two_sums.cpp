// cd OneDrive\PhD\Codes\Test mathematical expressions
// g++ 20221214_product_of_two_sums.cpp -o 20221214_product_of_two_sums

/* ( sum_{i=1}^N a_i ) * ( sum_{i=1}^N b_i )
= sum_{i=1}^N sum_{j=1}^N a_i b_j
= sum_{i=1}^N a_i b_i + sum_{i=1}^N sum_{j=i+1}^N ( a_i b_j + a_j b_i )
= sum_{i=1}^N ( a_i b_i = sum_{j=i+1}^N [ a_i b_j + a_j b_i ] )
*/

#include <iostream>
#include <stdlib.h>					// srand, rand
#include <vector>
#include "../Vector_Operations.h"	// print_vector();


//----- Global variables and function declarations -----


int main(){
	
	//----- Define variables -----
	int N = 3;					// Number of terms
	std::vector<double> A (N);
	std::vector<double> B (N);
	
	
	//----- Initialise random seed -----
	// Program output will be different every time. You don't even need to compile the code again for different results - just run it again.
	// Comment-out this line to get the same output every time.
	// Numbers won't change too much. Probably a better random-seed can be chosen.
	srand( time(0) );
	
	
	//----- Populuate the lists with random numbers -----
	for( int n=0; n<N; n++ ){
		A[n] = rand() % 10;
		B[n] = rand() % 10;
	}
	
	
	//----- Calculate the product of the sums of the two lists -----
	double sum_A = 0;
	double sum_B = 0;
	for( int n=0; n<N; n++ ){
		sum_A += A[n];
		sum_B += B[n];
	}
	double product_exact = sum_A * sum_B;
	
	double product_guess_1 = 0;
	for( int n1=0; n1<N; n1++ ){
		for( int n2=0; n2<N; n2++ ){
			product_guess_1 += A[n1] * B[n2] ;
		}
	}
	
	double product_guess_2 = 0;
	for( int n=0; n<N; n++ ){
		product_guess_2 += A[n] * B[n];
	}
	for( int n1=0; n1<N; n1++ ){
		for( int n2=n1+1; n2<N; n2++ ){
			//std::cout << "Indices:\t" << n1+1 <<"\t"<< n2+1 << std::endl;
			product_guess_2 += A[n1] * B[n2] + A[n2] * B[n1];
		}
	}
	
	double product_guess_3 = 0;
	for( int n1=0; n1<N; n1++ ){
		product_guess_3 += A[n1] * B[n1];
		for( int n2=n1+1; n2<N; n2++ ){
			product_guess_3 += A[n1] * B[n2] + A[n2] * B[n1];
		}
	}
	
	//--- 4: sum_{i=1}^n ( a_ib_i + a_i sum_{j=i+1}^n b_j + b_i sum_{j=i+1}^n a_j ---
	double product_guess_4 = 0;
	for( int i=0; i<N; i++ ){
		
		product_guess_4 += A[i] * B[i];
		
		double product_guess_4a = 0;
		for( int j=i+1; j<N; j++ ){
			product_guess_4a += B[j];
		}
		product_guess_4 += A[i] * product_guess_4a;
		
		double product_guess_4b = 0;
		for( int j=i+1; j<N; j++ ){
			product_guess_4b += A[j];
		}
		product_guess_4 += B[i] * product_guess_4b;
	}
	
	
	double product_guess_5 = 0;
	for( int n1=0; n1<N; n1++ ){
		
		double product_guess_5a = 0;
		for( int n2=n1; n2<N; n2++ ){
			product_guess_5a += B[n2];
		}
		product_guess_5 += A[n1] * product_guess_5a;
		
		double product_guess_5b = 0;
		for( int n2=n1+1; n2<N; n2++ ){
			product_guess_5b += A[n2];
		}
		product_guess_5 += B[n1] * product_guess_5b;
	}
	
	
	double product_guess_6 = 0;
	for( int n1=0; n1<N; n1++ ){
		
		double product_guess_6a = 0;
		for( int n2=n1+1; n2<N; n2++ ){
			product_guess_6a += B[n2];
		}
		product_guess_6 += A[n1] * product_guess_6a;
		
		double product_guess_6b = 0;
		for( int n2=n1; n2<N; n2++ ){
			product_guess_6b += A[n2];
		}
		product_guess_6 += B[n1] * product_guess_6b;
	}
	
	
	//----- Output results -----
	print_vector( A, "A" );
	print_vector( B, "B" );
	std::cout << "sum A:\t" << sum_A << std::endl;
	std::cout << "sum B:\t" << sum_B << std::endl;
	std::cout << "(sum A) times (sum B)\t:" << product_exact <<"\t"<< product_guess_1 <<"\t"<< product_guess_2 <<"\t"<< product_guess_3 <<"\t"<< product_guess_4 << "\t"<< product_guess_5 << std::endl;
	
	//std::cout << RAND_MAX <<std::endl;
	
	
	return 0;
}