/*
cd OneDrive\PhD\Codes\Check mathematical expressions
g++ 20230405_x_by_2_factorial.cpp -o 20230405_x_by_2_factorial
20230405_x_by_2_factorial

Recursive relation for
( x/2 )!

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
	
	std::vector<double> xby2f_exact ( 2*x_max );
	std::vector<double> xby2f_guess ( 2*x_max );
	
	
	//----- Create list of x-values -----
	// Just to be sure that indices are correct. Not used.
	std::vector<double> x_list (2*x_max+1 );
	
	for( int x=-x_max; x<=x_max; x++ ){
		x_list[x_max+x] = x;
	}
	cout_vector( x_list );
	std::cout << "Number of values of x:\t"<< x_list.size() << std::endl;
	
	
	
	
	//----- Calculate exact expression -----
	for( int x=-x_max; x<=x_max; x++
	){
		xby2f_exact[x_max+x] = tgamma( 0.5*x +1 );
	}
	
	
	//----- Alternative version using just one list -----
	xby2f_guess[x_max] = 1.0;
	xby2f_guess[x_max-1] = sqrt( pi );
	
	for( int x=1; x<=x_max; x++ ){
		xby2f_guess[x_max+x] = (double) 0.5*x * xby2f_guess[x_max+x-2];
	}
	for( int x=-3; x>=-x_max; x-=2 ){
		xby2f_guess[x_max+x] = 2.0/(x+2.0)    * xby2f_guess[x_max+x+2];
	}
	
	
	
	
	
	//----- Output results -----
	std::cout << std::endl;
	std::cout <<std::left<<std::setw(w1)<< "i" <<std::left<<std::setw(w1)<< "x" << "  |  "
			  <<std::left<<std::setw(w2)<< "xby2f_guess" <<std::left<<std::setw(w2)<< "xby2f_exact" << "  |  "
			  <<std::left<<std::setw(w2)<< "xby2f_guess - xby2f_exact" << std::endl;
				  
	int i = 0;	// List index
	for( int x=-x_max; x<=x_max; x++ ){
		std::cout <<std::left<<std::setw(w1)<< i <<std::left<<std::setw(w1)<< x << "  |  "
		          <<std::left<<std::setw(w2)<< xby2f_guess[x_max+x] <<std::left<<std::setw(w2)<< xby2f_exact[x_max+x] << "  |  "
		          <<std::left<<std::setw(w2)<< xby2f_guess[x_max+x] - xby2f_exact[x_max+x] << std::endl;
		i += 1;
	}
	
	
	//----- Test -----
	std::cout << std::endl;
	//int x_test = 0;
	
	for( int x_test=-x_max; x_test<=x_max; x_test++ ){
		int i_test = x_max + x_test;
		std::cout << i_test <<"\t"<< x_test <<"\t"<< xby2f_guess[i_test] <<"\t"<< tgamma( 0.5*x_test +1 ) << std::endl;
	}
	
	
	
	
	
	
	
	//----- Done -----
	
	std::cout << "\nDone" << std::endl;
	
	return 0;
}