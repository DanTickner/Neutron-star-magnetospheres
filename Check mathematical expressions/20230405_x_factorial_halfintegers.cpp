/*
cd OneDrive\PhD\Codes\Check mathematical expressions
g++ 20230405_x_factorial_halfintegers.cpp -o 20230405_x_factorial_halfintegers
20230405_x_factorial_halfintegers

Recursive relation for
( x )!

where x is an integer or half-integer, going from -x_max to x_max.

Book 15, p19.
*/

#include <iostream>
#include <iomanip>
#include <math.h>
#include <vector>

#include "../Vector_Operations.h"	// cout_vector()


int main(){
	
	//----- Define variables -----
	int w1 = 5;		// Fixed width for text output, indices.
	int w2 = 15;	// Fixed width for text output, factorials.
	int x_max = 10;
	
	const double pi = acos(-1);
	
	std::vector<double> factorials_exact ( 4*x_max+1 );
	std::vector<double> factorials_guess ( 4*x_max+1 );
	
	
	//----- Create list of x-values -----
	// Just to be sure that indices are correct. Not used.
	std::vector<double> x_list (4*x_max+1 );
	
	for( int i=0; i<=4*x_max; i++ ){
		double x = -x_max + 0.5*i;
		//std::cout << x << std::endl;
		x_list[i] = x;
	}
	//cout_vector( x_list );
	//std::cout << "Number of values of x:\t"<< x_list.size() << std::endl;
	
	
	
	
	//----- Calculate exact expression -----
	for( int i=0; i<=4*x_max; i++ ){
		double x = -x_max + 0.5*i;
		factorials_exact[i] = tgamma( x +1 );
	}
	//cout_vector( factorials_exact );
	
	
	//----- Alternative version using just one list -----
	factorials_guess[2*x_max] = 1.0;
	factorials_guess[2*x_max-1] = sqrt( pi );
	
	for( int i=2*x_max+1; i<=4*x_max; i++ ){
		double x = -x_max + 0.5*i;
		factorials_guess[i] = x * factorials_guess[i-2];
	}
	for( int i=2*x_max-3; i>=-x_max; i-=2 ){
		double x = -x_max + 0.5*i;
		factorials_guess[i] = 1.0/(x+2.0)    * factorials_guess[i+2];
		std::cout << i <<"\t"<< x <<"\t"<< 1.0/(x+2.0) <<"\t"<< factorials_guess[i+2] << std::endl;
	}
	
	
	
	//----- Output results -----
	std::cout << std::endl;
	std::cout <<std::left<<std::setw(w1)<< "i" <<std::left<<std::setw(w1)<< "x" << "  |  "
			  <<std::left<<std::setw(w2)<< "xby2f_guess" <<std::left<<std::setw(w2)<< "xby2f_exact" << "  |  "
			  <<std::left<<std::setw(w2)<< "xby2f_guess - xby2f_exact" << std::endl;
				  
				  
				  
	
	for( int i=0; i<=4*x_max; i++ ){
		double x = -x_max + 0.5*i;
		std::cout <<std::left<<std::setw(w1)<< i <<std::left<<std::setw(w1)<< x << "  |  "
		          <<std::left<<std::setw(w2)<< factorials_guess[i] <<std::left<<std::setw(w2)<< factorials_exact[i] << "  |  "
		          <<std::left<<std::setw(w2)<< factorials_guess[i] - factorials_exact[i] << std::endl;
	
	}
	/*
	
	int i = 0;	// List index
	for( int x=-x_max; x<=x_max; x++ ){
		std::cout <<std::left<<std::setw(w1)<< i <<std::left<<std::setw(w1)<< x << "  |  "
		          <<std::left<<std::setw(w2)<< factorials_guess[x_max+x] <<std::left<<std::setw(w2)<< factorials_exact[x_max+x] << "  |  "
		          <<std::left<<std::setw(w2)<< factorials_guess[x_max+x] - factorials_exact[x_max+x] << std::endl;
		i += 1;
	}
	/*
	
	//----- Test -----
	std::cout << std::endl;
	//int x_test = 0;
	
	for( int x_test=-x_max; x_test<=x_max; x_test++ ){
		int i_test = x_max + x_test;
		std::cout << i_test <<"\t"<< x_test <<"\t"<< xby2f_guess[i_test] <<"\t"<< tgamma( 0.5*x_test +1 ) << std::endl;
	}
	
	
	
	
	
	*/
	
	//----- Done -----
	
	std::cout << "\nDone" << std::endl;
	
	return 0;
}