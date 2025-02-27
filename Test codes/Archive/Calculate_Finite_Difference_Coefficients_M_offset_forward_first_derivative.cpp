/*
cd OneDrive\PhD\Codes\20240311 Time evolution code\Test codes
g++ Calculate_Finite_Difference_Coefficients_M_offset_forward_first_derivative.cpp -o Calculate_Finite_Difference_Coefficients_M_offset_forward_first_derivative --std=c++11
Calculate_Finite_Difference_Coefficients_M_offset_forward_first_derivative

Solve matrix equations to calculate the coefficients in the finite difference expressions for N_L gridpoints to the left and N_R gridpoints to the right of x.
Derived with N_L < N_R in mind, but might work for general N_L, N_R.
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
	int N_L    = 1;	// number of gridpoints to the left  we look.
	int N_R    = 2;	// number of gridpoints to the right we look. Must be greater than M for the logic in this code to work.
	int N_eqns = N_L + N_R - 1;	// The number of simultaneous equations.
	
	std::string output_filename_test_base = "20240426_Finite_difference_coefficients_M_offset_forward_first_derivative";	// Don't add ".csv" as the code will append _r, _t, etc to the filename for each csv.
	
	std::string output_filename_test_csv = "../CSV/"  + output_filename_test_base + "_" + std::to_string( N_L ) + "_" + std::to_string( N_R ) +  ".csv";
	
	std::ofstream output_file_test_csv;
	
	/*
	//----- Create CSV here so that we can output as we go along -----
	
	
	output_file_test_csv << "order,method,n_points";
	
	for( int i=0; i<N_max; i++ ){
		output_file_test_csv << ",tilde_a_" << i;
	}
	
	output_file_test_csv << "\n";
	*/
	
	
	//----- Build the matrices -----
	// Here, we use i to index the row and j to index the column, starting from 1, which is the standard mathematical convention.
	// C++ starts from 0, so we subtract 1 to keep the mathematical form of the equation visible for easier checking.
	
	std::vector< std::vector<double> > M ( N_eqns, std::vector<double> ( N_eqns ) );
	std::vector<std::string> A ( N_eqns );	// This is for output only, since the gaussian_elimination function provides us with a vector.
	std::vector<double> N ( N_eqns );
	
	for( int i=1; i<=N_eqns; i++ ){
		
		//--- M ---
		for( int j=1; j<=N_eqns; j++ ){
			
			if( j <= N_L ){
				M[i-1][j-1] = pow( j - N_L - 1, i );
			}
			else{
				M[i-1][j-1] = pow( j - N_L, i );
			}
		}
			
		//--- A ---
		int a_subscript = 0;
		
		if( i <= N_L ){
			a_subscript = i - N_L - 1;
		}
		else{
			a_subscript = i - N_L;
		}
		
		A[i-1] = "\\tilde a_{" + std::to_string( a_subscript ) + "}";
		
		//--- N ---
		N[i-1] = - pow( N_R, i );
			
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
	}
	
	
	//----- Calculate the expression for the derivative -----
	double a_N_R = 1.0;	// change this result to guarantee integers.
	
	std::vector<double> a_n ( N_eqns + 1 );
	
	
	for( int n=-N_L; n<0; n++ ){
		a_n[N_L+n] = tilde_a_n[n] * a_N_R;
	}
	a_n[N_L] = 0;	// This is already true from the declaration of a_n so not raelly needed.
	for( int n=1; n<N_R; n++ ){
		a_n[N_L+n] = tilde_a_n[n] * a_N_R;
	}
	a_n.back() = a_N_R;
	
	for( int n=-N_L; n<=N_R; n++ ){
		if( N_L+n == 0 ){ continue; }
		std::cout << "a_{" << n << "}:\t"<< a_n[N_L+n] << std::endl;
	}
	
	double f_of_x_coeff = 0;
	
	for( int n=-N_L; n<=N_R; n++ ){
		
		if( n == 0 ){
			continue;
		}
		
		std::cout << N_L + n <<"\t"<< a_n[ N_L + n ] << std::endl;
		
		f_of_x_coeff += a_n[ N_L + n ];
	}
	
	std::cout << "f_of_x_coeff:\t" << f_of_x_coeff << std::endl;
		
		
		
	
	
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