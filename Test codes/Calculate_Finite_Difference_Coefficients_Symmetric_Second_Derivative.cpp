/*
cd OneDrive\PhD\Codes\20240311 Time evolution code\Test codes
g++ Calculate_Finite_Difference_Coefficients_Symmetric_Second_Derivative.cpp -o Calculate_Finite_Difference_Coefficients_Symmetric_Second_Derivative --std=c++11
Calculate_Finite_Difference_Coefficients_Symmetric_Second_Derivative

Solve matrix equations to calculate the coefficients in the finite difference expressions for N gridpoints to the left of x and to the right of x, all equally spaced.
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
	int N = 3;	// number of gridpoints backward and forward we look.
	
	double a_N = 2;	// change this result to guarantee integers.
	
	std::string output_filename_test_base = "20240501_Finite_difference_coefficients_symmetric_second_derivative";	// Don't add ".csv" as the code will append _r, _t, etc to the filename for each csv.
	
	std::string output_filename_test_csv = "../CSV/"  + output_filename_test_base + "_" + std::to_string( N ) +  ".csv";
	
	std::ofstream output_file_test_csv;
	
	
	//----- Build the matrices -----
	// Here, we use i to index the row and j to index the column, starting from 1, which is the standard mathematical convention.
	// C++ starts from 0, so we subtract 1 to keep the mathematical form of the equation visible for easier checking.
	
	std::vector< std::vector<double> > M ( N - 1, std::vector<double> ( N - 1 ) );
	std::vector<std::string> A ( N - 1 );	// This is for output only, since the gaussian_elimination function provides us with a vector.
	std::vector<double> N_vector ( N - 1 );	// Add _vector to avoid confusion with integer N.
	
	
	//--- M ---
	for( int i=1; i<=N-1; i++ ){
		for( int j=1; j<=N-1; j++ ){
			M[i-1][j-1] = pow( j, 2.0*i+2.0 );
		}
	}
	
	//--- A ---
	for( int i=1; i<=N-1; i++ ){
		int a_subscript = i;
		A[i-1] = "a_{" + std::to_string( a_subscript ) + "}";
	}
	
	//--- N ---
	for( int i=1; i<=N-1; i++ ){
		N_vector[i-1] = - pow( N, 2.0*i+2.0 );
	}
	
	
	
	//----- Output to CSV -----
	// Now we can revert to i,j starting from 0 since they're only indices.
	
	output_file_test_csv.open( output_filename_test_csv );
	for( int i=0; i<N-1; i++ ){
		
		output_file_test_csv << "(,";
		
		for( int j=0; j<N-1; j++ ){
			output_file_test_csv << M[i][j] <<",";
		}
		
		output_file_test_csv << ") dot (," << A[i] <<",) = (," << N_vector[i] << ",)\n";
		
	}
	
	
	//----- Calculate the coefficients -----
	std::vector<double> tilde_a_n = gaussian_elimination( M, N_vector );
	
	for( int i=1; i<=N-1; i++ ){
		std::cout << A[i-1] << ":\t"<< tilde_a_n[i-1] << std::endl;
		//std::cout << A[i-1] << ":\t"<<std::setprecision(15)<< tilde_a_n[i-1] << std::endl;	// Use when the conversion of the stated decimal to a fraction isn't obvious.
	}
	std::cout << std::endl;
	
	
	//----- Calculate the expression for the derivative -----
	std::cout << "Setting a_{" << N << "} = " << a_N << ". Change this value and run the code if this doesn't match expectations.\n" << std::endl;
	std::vector<double> a_n ( N );
	
	for( int i=0; i<N-1; i++ ){
		a_n[i] = tilde_a_n[i] * a_N;
	}
	
	a_n.back() = a_N;
	
	double f_of_x_coeff = 0;
	double denominator  = 0;
	
	for( int n=1; n<=N; n++ ){
		
		int index = n - 1;	// Making space for a_0 to be added once its calculation is complete after this loop.
		
		f_of_x_coeff += a_n[ index ] * 2.0;
		denominator  += a_n[ index ] * n * n;
	}
	
	a_n.insert( a_n.begin(), f_of_x_coeff );	// Need a_n.begin() because it turns N_B from an integer into an iterator, which is the variable type required by std::vector.insert().
	
	std::cout << "\nFinal a_n (except absorbing sign). a_{-n}=a_n not shown. Index then value:" << std::endl;
	for( int i=0; i<=N; i++ ){
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
	
	
	
	//----- Output the mathematical expression -----
	// I prefer the convention where you have the expression going from most-forward to most-backward.
	
	std::cout << "\nf''(x) = 1 / ( " << denominator << "* h^2 ) * [" << std::endl;	// General case.
	
	//--- n > 0 ---
	for( int n=N; n>0; n-- ){
		
		int index = n;
		
		( a_n[ index ] >= 0 ) ? ( std::cout << " + " ) : ( std::cout << " - " );
		
		std::cout << abs( a_n[ index ] ) << " f(x+";
		
		if( n > 1 ){
			std::cout << n;
		}
		
		std::cout << "h)";
	}
	
	//--- n = 0 ---
	// a_0 is subtracted so the sign should be opposite.	
	( a_n[ 0 ] >= 0 ) ? ( std::cout << " - " ) : ( std::cout << " + " );
	
	std::cout << abs( a_n[ 0 ] ) << " f(x)";
	
	//--- n < 0 --
	for( int n=-1; n>=-N; n-- ){
		
		int index = -n;
		
		( a_n[ index ] >= 0 ) ? ( std::cout << " + " ) : ( std::cout << " - " );
		
		std::cout << abs( a_n[ index ] ) << " f(x-";
		
		if( n < -1 ){
			std::cout << abs( n );
		}
		
		std::cout << "h)";
	}
	
	std::cout << "\n] + Order( h^" << 2 * N << " )." << std::endl;
	
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