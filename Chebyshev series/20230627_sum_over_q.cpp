/*

cd OneDrive\PhD\Codes\Chebyshev series
g++ 20230627_sum_over_q.cpp -o 20230627_sum_over_q
20230627_sum_over_q

Can we write
sum_{q=n+1//n+q odd}^{k-1} q
as some number depending only on n,k?

*/


#include <iostream>
#include <math.h>


int main(){
	
	//----- Define variables -----
	int n_max = 17;						// Maximum order for Chebyshev series.
	
	
	//----- Perform the summation for the relevant n, k -----
	std::cout << "n\tk\tsum\tguess\tdiff." << std::endl;
	
	for( int n=0; n<=n_max; n++ ){
		for( int k=n+2; k<=n_max; k++ ){
			if( (n+k)%2 == 0 ){
				
				double sum = 0;
				for( int q=n+1; q<=k-1; q++ ){
					if( (n+q)%2 == 1 ){
						sum += q;
					}
				}
				
				double sum_guess = 0.25 * ( k*k - n*n );
				
				std::cout << n <<"\t"<< k <<"\t"<< sum <<"\t"<< sum_guess <<"\t"<< sum-sum_guess << std::endl;
				
			}
		}
	}
	
	
	std::cout << "\nDone" << std::endl;
	
	return 0;
}