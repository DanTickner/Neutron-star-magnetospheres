/*
cd OneDrive\PhD\Codes\Check mathematical expressions
g++ 20230404_a_pm_b_m1_over2_factorial.cpp -o 20230404_a_pm_b_m1_over2_factorial
20230404_a_pm_b_m1_over2_factorial

Recursive relations for
( (a+b-1)/2 )!
and
( (a-b-1)/2 )!
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
	int a_max = 7;
	int b_min = 0;
	int b_max = 7;
	
	const double pi = acos(-1);
	
	std::vector< std::vector<double> > apb_exact ( a_max-a_min+1, std::vector<double> ( b_max-b_min+1 ) );
	std::vector< std::vector<double> > amb_exact ( a_max-a_min+1, std::vector<double> ( b_max-b_min+1 ) );
	std::vector< std::vector<double> > apb_guess ( a_max-a_min+1, std::vector<double> ( b_max-b_min+1 ) );
	std::vector< std::vector<double> > amb_guess ( a_max-a_min+1, std::vector<double> ( b_max-b_min+1 ) );
	std::vector<double> xby2f   (a_max-a_min+1);	// ( x/2 )!
	
	//------ Calculate exact values -----
	for( int a=a_min; a<=a_max; a++ ){
		for( int b=b_min; b<=b_max; b++ ){
			apb_exact[a][b] = tgamma( 0.5*(a+b-1.0) +1 );
			amb_exact[a][b] = tgamma( 0.5*(a-b-1.0) +1 );
		}
	}
	
	
	//----- Calculate guessed values by recursive relation -----
	
	//--- Generate ( x/2 )! ---
	xby2f[0] = 1;
	xby2f[1] = 0.5 * sqrt( pi );
	for( int x=2; x<=a_max; x++ ){
		xby2f[x] = 0.5*x * xby2f[x-2];
	}
	
	
	//--- ( (a+b-1)/2 )! ---
	/*	for( int a=0; a<=a_max; a++ ){	
		
		//apb_guess[a][0] = tgamma( 0.5*(a-1.0) +1 );
		apb_guess[a][0] = xby2f[a-1];
		
		//apb_guess[a][1] = tgamma( 0.5*a       +1 );
		apb_guess[a][1] = xby2f[a];
		
		//amb_guess[a][0] = xby2f[a-1];
		//amb_guess[a][1] = xby2f[a-2];
		
		for( int b=2; b<=b_max; b++ ){
			apb_guess[a][b] = apb_guess[a][b-2] * 0.5*(a+b-1.0);
			//amb_guess[a][b] = amb_guess[a][b-2] * 0.5*(a-b-1.0);
			//std::cout << a <<"\t"<< b <<"\t"<< apb_guess[a][b-2] <<"\t"<< 0.5*(a+b-1.0) << std::endl;
		}
	}*/
	
	// a=0.
	apb_guess[0][0] = sqrt(pi);
	apb_guess[0][1] = 1.0;
	for( int b=2; b<=b_max; b++ ){
		apb_guess[0][b] = apb_guess[0][b-2] * 0.5*(b-1.0);
	}
	
	// a>0.
	for( int a=1; a<=a_max; a++ ){	
		apb_guess[a][0] = xby2f[a-1];
		apb_guess[a][1] = xby2f[a];
		
		for( int b=2; b<=b_max; b++ ){
			apb_guess[a][b] = apb_guess[a][b-2] * 0.5*(a+b-1.0);
		}
	}
	
	//--- ( (a-b-1)/2 )! ---
	/*for( int a=0; a<=a_max; a++ ){
		amb_guess[a][a-1] = 1.0;			// WHAT HAPPENS IF A=0 HERE? DO AGAIN
		amb_guess[a][a  ] = sqrt(pi);
		
		//--- Calculate from b=a-2 down to b=0 ---
		for( int b=a-2; b>=0; b-- ){
			amb_guess[a][b] = (a-b-1.0)/2.0 * amb_guess[a][b+2];
		}
		
		//--- Calculate from b=a+1 up to b=b_max ---
		for( int b=a+1; b<=b_max; b++ ){
			amb_guess[a][b] = 2.0/(a-b+1.0) * amb_guess[a][b-2];
			if( isinf( amb_guess[a][b] ) ){
				amb_guess[a][b] = 0;
			}
		}
	}*/
	
	// a=0 must be considered separately, or we hit an indexing issue with amb_guess[a][a-1].
	amb_guess[0][0] = sqrt( pi );
	for( int b=2; b<=b_max; b+=2 ){
		amb_guess[0][b] = 2.0/(-b+1.0) * amb_guess[0][b-2];
	}
	
	// Remainder of the cases a>0.
	for( int a=1; a<=a_max; a++ ){
		amb_guess[a][a-1] = 1.0;
		amb_guess[a][a  ] = sqrt(pi);
		
		//--- Calculate from b=a-2 down to b=0 ---
		for( int b=a-2; b>=0; b-- ){
			amb_guess[a][b] = (a-b-1.0)/2.0 * amb_guess[a][b+2];
		}
		
		//--- Calculate from b=a+1 up to b=b_max ---
		for( int b=a+1; b<=b_max; b++ ){
			amb_guess[a][b] = 2.0/(a-b+1.0) * amb_guess[a][b-2];
			if( isinf( amb_guess[a][b] ) ){
				amb_guess[a][b] = 0;
			}
		}
	}
	
	
	
	//----- Output results -----
	for( int a=a_min; a<=a_max; a++ ){
		for( int b=b_min; b<=b_max; b++ ){
			std::cout <<std::left<<std::setw(w1)<< a <<std::left<<std::setw(w1)<< b <<"  |  "
			          <<std::left<<std::setw(w2)<< apb_guess[a][b] <<std::left<<std::setw(w2)<< apb_exact[a][b]
					  <<std::left<<std::setw(w2)<< apb_exact[a][b] - apb_guess[a][b] <<"  |  "
					  <<std::left<<std::setw(w2)<< amb_guess[a][b] <<std::left<<std::setw(w2)<< amb_exact[a][b]
					  <<std::left<<std::setw(w2)<< amb_exact[a][b] - amb_guess[a][b] << std::endl;
		}
		std::cout << std::endl;
	}
	
	std::cout << "Done" << std::endl;
	
	return 0;
}