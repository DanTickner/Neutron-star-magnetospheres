/*
cd OneDrive\PhD\Codes\20240311 Time evolution code\Test codes
g++ Calculate_Finite_Difference_Coefficients.cpp -o Calculate_Finite_Difference_Coefficients --std=c++11
Calculate_Finite_Difference_Coefficients

Solve matrix equations to calculate the coefficients in the finite difference expressions for n gridpoints, using a forward method, a backward method and a symmetric method.
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
	std::string output_filename_test_base = "20240621_Finite_difference_coefficients";	// Don't add ".csv" as the code will append _r, _t, etc to the filename for each csv.
	
	std::string output_filename_test_csv = "../CSV/"  + output_filename_test_base + ".csv";
	std::string output_filename_test_csv_symmetric_matrix = "../CSV/"  + output_filename_test_base + "_symmetric_matrix.csv";
	
	std::ofstream output_file_test_csv;
	std::ofstream output_file_test_csv_symmetric_matrix;
	
	int N_max = 7;	// Stop generating coefficients once we have considered this many gridpoints. This value is N+1, where N is the maximum distance away from x in a forward or backward difference equation.
	
	
	//----- Create CSV here so that we can output as we go along -----
	output_file_test_csv.open( output_filename_test_csv );
	output_file_test_csv_symmetric_matrix.open( output_filename_test_csv_symmetric_matrix );
	
	output_file_test_csv << "order,method,n_points";
	
	for( int i=0; i<N_max; i++ ){
		output_file_test_csv << ",tilde_a_" << i;
	}
	
	output_file_test_csv << "\n";
	
	
	
	
	//----- Forward and backward difference methods -----
	
	for( int N=1; N<=N_max; N++ ){
		// N is the maximum deviation in grid spacing from x. Then, we have n+1 gridpoints in a forward/backward difference model and ___ gridpoints in a symmetric difference model.
		
		std::cout << "\nFirst-order, forward and backward, N = " << N << std::endl;
		
		output_file_test_csv << "first,forward_and_backward," << N;
		
		
		std::vector< std::vector<double> > A ( N-1, std::vector<double> ( N-1 ) );
		std::vector<double> b ( N-1 );
		
		for( int i=1; i<=N-1; i++ ){
			
			for( int j=1; j<=N-1; j++ ){
				A[i-1][j-1] = pow( (double) j, i+1.0 );
			}
			
			b[i-1] = - pow( (double) N, i+1.0 );
		}
		
		//--- Print matrix and vector to the screen. Optional. ---
		for( int i=1; i<=N-1; i++ ){
			
			std::cout << "A_" << i << ":";
			
			for( int j=1; j<=N-1; j++ ){
				std::cout << "\t"<< A[i-1][j-1];
			}
			
			std::cout << std::endl;
		}
		
		std::cout << "\nb:";
		
		for( int j=1; j<=N-1; j++ ){
			std::cout << "\t"<< b[j-1];
		}
		
		std::cout <<"\n"<< std::endl;
		
		/*
		// Check for case N=3.
		A[0][0] = 1.0;
		A[0][1] = 4.0;
		
		A[1][0] = 1.0;
		A[1][1] = 8.0;
		
		b[0] = -9.0;
		b[1] = -27.0;
		*/
		
		//--- Output to screen and CSV ---
		std::vector<double> tilde_a_n = gaussian_elimination( A, b );
		
		for( int n=1; n<=N-1; n++ ){
			
			std::cout << "tilde a_" << n << ":\t" << tilde_a_n[n-1] << std::endl;
			
			output_file_test_csv <<","<< tilde_a_n[n-1];
			
		}
		output_file_test_csv << "\n";
		
	}
	
	
	
	
	//----- Symmetric difference method -----
	for( int N=1; N<=N_max; N++ ){
		// N is the maximum deviation in grid spacing from x. Then, we have n+1 gridpoints in a forward/backward difference model and ___ gridpoints in a symmetric difference model.
		std::cout << "code gets here: symmmetric N=" << N << std::endl;
		
		std::cout << "\nFirst-order, symmetric, N = " << N << std::endl;
		
		output_file_test_csv << "first,symmetric," << N;
		
		
		std::vector< std::vector<double> > A ( N-1, std::vector<double> ( N-1 ) );
		std::vector<double> b ( N-1 );
		
		for( int i=1; i<=N-1; i++ ){
			
			for( int j=1; j<=N-1; j++ ){
				A[i-1][j-1] = pow( (double) j, 2.0*i+1.0 );
			}
			
			b[i-1] = - pow( (double) N, 2.0*i+1.0 );
		}
		
		//--- Print matrix and vector to the screen. Optional. ---
		for( int i=1; i<=N-1; i++ ){
			
			std::cout << "A_" << i << ":";
			
			for( int j=1; j<=N-1; j++ ){
				std::cout <<"\t"<< A[i-1][j-1];
			}
			
			std::cout << std::endl;
		}
		
		std::cout << "\nb:";
		
		for( int j=1; j<=N-1; j++ ){
			std::cout << "\t"<< b[j-1];
		}
		
		std::cout <<"\n"<< std::endl;
		
		
		//--- Output to screen and CSV ---
		std::vector<double> tilde_a_n = gaussian_elimination( A, b );
		
		for( int n=1; n<=N-1; n++ ){
			
			std::cout << "tilde a_" << n << ":\t" << tilde_a_n[n-1] << std::endl;
			
			output_file_test_csv <<","<< tilde_a_n[n-1];
			
		}
		output_file_test_csv << "\n";
		
		//--- Output matrix to CSV ---
		for( int i=0; i<N-1; i++ ){
		
			output_file_test_csv_symmetric_matrix << "(,";
			
			for( int j=0; j<N-1; j++ ){
				output_file_test_csv_symmetric_matrix << A[i][j] <<",";
			}
			
			output_file_test_csv_symmetric_matrix << ") dot (," << tilde_a_n[i] <<",) = (," << b[i] << ",)\n";
		
		}
		
	}
	
	
	
	
	//----- 1-offset forward method -----
	for( int N=1; N<=N_max; N++ ){
		// N is the maximum deviation in grid spacing from x. Then, we have n+1 gridpoints in a forward/backward difference model and ___ gridpoints in a symmetric difference model.
		
		std::cout << "\nFirst-order, 1-offset forward, N = " << N << std::endl;
		
		output_file_test_csv << "first,1_offset_forward," << N;
		
		
		std::vector< std::vector<double> > A ( N-1, std::vector<double> ( N-1 ) );
		std::vector<double> b ( N-1 );
		
		for( int i=1; i<=N-1; i++ ){
			
			for( int j=1; j<=N-1; j++ ){
				
				if( j == 1 ){
					A[i-1][j-1] = pow( -1, i+1 );
				}
				
				else {
					A[i-1][j-1] = pow( (double) j-1.0, i+1.0 );
				}
				
			}
			
			b[i-1] = - pow( (double) N, i+1.0 );
		}
		
		//--- Print matrix and vector to the screen. Optional. ---
		for( int i=1; i<=N-1; i++ ){
			
			std::cout << "A_";
			if( i == 1 ){
				std::cout << "-1";
			}
			else {
				std::cout << i-1;
			}
			std::cout << ":";
			
			for( int j=1; j<=N-1; j++ ){
				std::cout << "\t"<< A[i-1][j-1];
			}
			
			std::cout << std::endl;
		}
		
		std::cout << "\nb:";
		
		for( int j=1; j<=N-1; j++ ){
			std::cout << "\t"<< b[j-1];
		}
		
		std::cout <<"\n"<< std::endl;
		
		
		//--- Output to screen and CSV ---
		std::vector<double> tilde_a_n = gaussian_elimination( A, b );
		
		for( int n=1; n<=N-1; n++ ){
			
			std::cout << "tilde a_" << n << ":\t" << tilde_a_n[n-1] << std::endl;
			
			output_file_test_csv <<","<< tilde_a_n[n-1];
			
		}
		output_file_test_csv << "\n";
		
	}
	
	
	
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