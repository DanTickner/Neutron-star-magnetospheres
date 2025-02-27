/*
cd OneDrive\PhD\Codes\20240311 Time evolution code\Test codes
g++ Calculate_Finite_Difference_Coefficients_General_Second_Derivative.cpp -o Calculate_Finite_Difference_Coefficients_General_Second_Derivative --std=c++11
Calculate_Finite_Difference_Coefficients_General_Second_Derivative

Solve matrix equations to calculate the coefficients in the finite difference expressions for N_B gridpoints to the left of x and N_F gridpoints to the right, all equally spaced.

Recovers the symmetric and forward expressions.
Doesn't seem to work for backward expressions, but anything with N_B > 1 is good.
*/

#include <iostream>
#include <iomanip>
#include <vector>
#include <math.h>
#include <fstream>
#include <chrono>


//----- Declare functions -----
std::vector<double> gaussian_elimination( std::vector< std::vector<double> > A, std::vector<double> b );



//----- Main code -----
int main(){
	
	
	//----- Extra variables useful only for this testing code -----
	int N_B    = 0;	// number of gridpoints backward we look.
	int N_F    = 8;	// number of gridpoints forward  we look.
	int N_eqns = N_B + N_F - 1;	// The number of simultaneous equations.
	
	double a_N_F = 1;	// change this result to guarantee integers.
	
	std::string output_filename_test_base = "20240502_Finite_difference_coefficients_general_second_derivative";	// Don't add ".csv" as the code will append _r, _t, etc to the filename for each csv.
	
	std::string output_filename_test_csv = "../CSV/"  + output_filename_test_base + "_" + std::to_string( N_B ) + "_" + std::to_string( N_F ) +  ".csv";
	
	std::ofstream output_file_test_csv;
	
	
	//----- Build the matrices -----
	// Here, we use i to index the row and j to index the column, starting from 1, which is the standard mathematical convention.
	// C++ starts from 0, so we subtract 1 to keep the mathematical form of the equation visible for easier checking.
	
	std::vector< std::vector<double> > M ( N_eqns, std::vector<double> ( N_eqns ) );
	std::vector<std::string> A ( N_eqns );	// This is for output only, since the gaussian_elimination function provides us with a vector.
	std::vector<double> N ( N_eqns );
	
	
	//--- M ---
	//- i = 1 -
	for( int j=1; j<=N_eqns; j++ ){
		M[1-1][j-1] = ( j <= N_B ) ? ( - N_B + j - 1 ) : ( - N_B + j );
	}
	
	//- i > 1 -
	for( int i=2; i<=N_eqns; i++ ){
		for( int j=1; j<=N_eqns; j++ ){
			double base = ( j <= N_B ) ? ( - N_B + j - 1 ) : ( - N_B + j );
			M[i-1][j-1] = pow( base, i + 1 );
		}
	}
	
	//--- A ---
	for( int i=1; i<=N_eqns; i++ ){
		int a_subscript = ( i <= N_B ) ? ( - N_B + i - 1 ) : ( - N_B + i );
		A[i-1] = "\\tilde a_{" + std::to_string( a_subscript ) + "}";
	}
	
	//--- N ---
	N[1-1] = - N_F;
	for( int i=2; i<=N_eqns; i++ ){
		N[i-1] = - pow( N_F, i + 1 );
	}
	
	
	
	//----- Output to CSV -----
	// Now we can revert to i,j starting from 0 since they're only indices.
	
	output_file_test_csv.open( output_filename_test_csv );
	for( int i=0; i<N_eqns; i++ ){
		
		output_file_test_csv << "(,";
		
		for( int j=0; j<N_eqns; j++ ){
			output_file_test_csv << M[i][j] <<",";
		}
		
		output_file_test_csv << ") dot (," << A[i] <<",) = (," << N[i] << ",)\n";
		
	}
	
	
	//----- Calculate the coefficients -----
	std::vector<double> tilde_a_n = gaussian_elimination( M, N );
	
	for( int i=1; i<=N_eqns; i++ ){
		std::cout << A[i-1] << ":\t"<< tilde_a_n[i-1] << std::endl;
		//std::cout << A[i-1] << ":\t"<<std::setprecision(15)<< tilde_a_n[i-1] << std::endl;	// Use when the conversion of the stated decimal to a fraction isn't obvious.
	}
	std::cout << std::endl;
	
	
	//----- Calculate the expression for the derivative -----
	std::cout << "Setting a_{" << N_F << "} = " << a_N_F << ". Change this value and run the code if this doesn't match expectations.\n" << std::endl;
	std::vector<double> a_n ( N_eqns + 1 );
	
	for( int i=0; i<N_eqns; i++ ){
		a_n[i] = tilde_a_n[i] * a_N_F;
	}
	
	a_n.back() = a_N_F;
	
	double f_of_x_coeff = 0;
	double denominator  = 0;
	
	for( int n=-N_B; n<=N_F; n++ ){
		
		if( n == 0 ){
			continue;
		}
		
		int index = ( n < 0 ) ? ( N_B + n ) : ( N_B + n - 1 );	// converts a_n to the correct index of the vector a.
		
		f_of_x_coeff += a_n[ index ];
		denominator  += a_n[ index ] * n * n;
		std::cout << a_n[index] <<"\t"<< n <<"\t"<< a_n[ index ] * n * n <<"\t"<< denominator << std::endl;
	}
	
	a_n.insert( a_n.begin()+N_B, f_of_x_coeff );	// Need a_n.begin() because it turns N_B from an integer into an iterator, which is the variable type required by std::vector.insert().
	
	std::cout << "\nFinal a_n (except absorbing sign and evenness of denominator). Index then value:" << std::endl;
	for( int i=-N_B; i<=N_F; i++ ){
		std::cout << i << "\t";
	}
	std::cout << std::endl;
	for( int i=0; i<a_n.size(); i++ ){
		std::cout << a_n[i] << "\t";
	}
	std::cout << std::endl;
	
	std::cout << "\nf_of_x_coeff:\t" << f_of_x_coeff << std::endl;
	std::cout <<   "denominator :\t" << denominator  << std::endl;
	
	//----- Absorb the sign of the denominator -----
	
	if( denominator < 0 ){
		
		denominator  = - denominator;
		f_of_x_coeff = - f_of_x_coeff;
		
		for( int i=0; i<a_n.size(); i++ ){
			a_n[i] = - a_n[i];
		}
		
	}
	
	//----- Absorb the half in the denominator; otherwise double the numerator -----
	// I think this works but I don't want any more uncertainty in the expressions I find, so comment-out for now and do this manually.
	/*
	bool is_denominator_integer = ( abs( denominator - int( denominator ) ) < 1e-10 );
	if( ( is_denominator_integer ) and ( int( denominator ) % 2 == 0 ) ){
		std::cout << "Denominator is even" << std::endl;
		denominator  = denominator / 2;
	}
	else {
		for( int i=0; i<a_n.size(); i++ ){
			a_n[i] = a_n[i] * 2;
		}
	}*/
	
	//----- Output the mathematical expression -----
	// I prefer the convention where you have the expression going from most-forward to most-backward.
	
	std::cout << "\nf''(x) = 1 / ( " << denominator << "* 0.5 * h^2 ) * [" << std::endl;	// Case when the half is not absorbed into the denominator or the coeffs. Comment-out if letting the code check for evenness.
	//std::cout << "\nf''(x) = 1 / ( " << denominator << "* 0.5 * h^2 ) * [" << std::endl;	// General case.
	
	//--- n > 0 ---
	for( int n=N_F; n>0; n-- ){
		
		int index = N_B + n;
		
		( a_n[ index ] >= 0 ) ? ( std::cout << " + " ) : ( std::cout << " - " );
		
		std::cout << abs( a_n[ index ] ) << " f(x+";
		
		if( n > 1 ){
			std::cout << n;
		}
		
		std::cout << "h)";
	}
	
	//--- n = 0 ---
	// a_0 is subtracted so the sign should be opposite.	
	( a_n[ N_B ] >= 0 ) ? ( std::cout << " - " ) : ( std::cout << " + " );
	
	std::cout << abs( a_n[ N_B ] ) << " f(x)";
	
	//--- n < 0 --
	for( int n=-1; n>=-N_B; n-- ){
		
		int index = N_B + n;
		
		( a_n[ index ] >= 0 ) ? ( std::cout << " + " ) : ( std::cout << " - " );
		
		std::cout << abs( a_n[ index ] ) << " f(x-";
		
		if( n < -1 ){
			std::cout << abs( n );
		}
		
		std::cout << "h)";
	}
	
	std::cout << "\n] + Order( h^" << N_B + N_F - 1 << " )." << std::endl;
	
	
	//----- Code finished -----
	
	std::cout << "\nCSV file saved:\t" << output_filename_test_csv << std::endl;
	
	std::cout << "\nCode finished" << std::endl;
	
	return 0;
	
}


//----- Define functions -----
std::vector<double> gaussian_elimination( std::vector< std::vector<double> > A, std::vector<double> b ){
	// solve the system of equations written in matrix form as Ax=b,
	// where x are the coefficients to be found.
	
	//----- Elimination phase -----
	for( int j=0; j<b.size(); j++ ){
		for( int i=j+1; i<b.size(); i++ ){
			if( A[i][j] != 0 ){
				double lambda = A[i][j] / A[j][j];
				for( int k=j+1; k<b.size(); k++ ){
					A[i][k] -= lambda * A[j][k];
				}
				b[i] -= lambda * b[j];
			}
		}
	}
	
	//----- Back-substitution phase -----
	for( int i=b.size()-1; i>-1; i-- ){
		double dotproduct = 0;
		for( int j=i+1; j<b.size(); j++ ){
			dotproduct += A[i][j] * b[j];
		}
		b[i] = ( b[i] - dotproduct ) / A[i][i];
	}
	
	return b;
	
}