/*
cd OneDrive\PhD\Codes\Check mathematical expressions
g++ 20230404_a_minus_b_factorial_and_b_minus_a_factorial.cpp -o 20230404_a_minus_b_factorial_and_b_minus_a_factorial
20230404_a_minus_b_factorial_and_b_minus_a_factorial

Recursive relations for
(a-b)!
and
(b-a)!
for a,b both increasing from 0 to some maximum.

Book 15, p19.
*/

#include <iostream>
#include <iomanip>
#include <math.h>
#include <vector>


int main(){
	
	//----- Define variables -----
	int w1 = 5;		// Fixed width for text output, indices.
	int w2 = 15;	// Fixed width for text output, factorials.
	int a_min = 0;
	int a_max = 4;
	int b_min = 0;
	int b_max = 4;
	
	std::vector< std::vector<double> > ambf_exact ( a_max-a_min+1, std::vector<double> ( b_max-b_min+1 ) );
	std::vector< std::vector<double> > bmaf_exact ( a_max-a_min+1, std::vector<double> ( b_max-b_min+1 ) );
	std::vector< std::vector<double> > ambf_guess ( a_max-a_min+1, std::vector<double> ( b_max-b_min+1 ) );
	std::vector< std::vector<double> > bmaf_guess ( a_max-a_min+1, std::vector<double> ( b_max-b_min+1 ) );
	
	//------ Calculate exact values -----
	for( int a=a_min; a<=a_max; a++ ){
		for( int b=b_min; b<=b_max; b++ ){
			ambf_exact[a][b] = tgamma( a-b +1 );
			bmaf_exact[a][b] = tgamma( b-a +1 );
		}
	}
	
	
	//----- Calculate guessed values by recursive relation -----
	
	for( int a=0; a<=a_max; a++ ){	
		ambf_guess[a][a] = 1.0;
		for( int b=a-1; b>=0; b-- ){
			ambf_guess[a][b] = ambf_guess[a][b+1] * ((double)a-b);
		}
	}
	
	
	
	//----- Output results -----
	for( int a=a_min; a<=a_max; a++ ){
		for( int b=b_min; b<=b_max; b++ ){
			std::cout <<std::left<<std::setw(w1)<< a <<std::left<<std::setw(w1)<< b <<"  |  "
			          <<std::left<<std::setw(w2)<< ambf_guess[a][b] <<std::left<<std::setw(w2)<< ambf_exact[a][b]
					  <<std::left<<std::setw(w2)<< ambf_exact[a][b] - ambf_guess[a][b] <<"  |  "
					  <<std::left<<std::setw(w2)<< bmaf_guess[a][b] <<std::left<<std::setw(w2)<< bmaf_exact[a][b]
					  <<std::left<<std::setw(w2)<< bmaf_exact[a][b] - bmaf_guess[a][b] << std::endl;
		}
	}
	
	std::cout << "Done" << std::endl;
	
	return 0;
}