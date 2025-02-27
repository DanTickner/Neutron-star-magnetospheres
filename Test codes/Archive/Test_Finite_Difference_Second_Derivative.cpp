/*
cd OneDrive\PhD\Codes\20240311 Time evolution code\Test codes
g++ Test_Finite_Difference_Second_Derivative.cpp -o Test_Finite_Difference_Second_Derivative --std=c++11
Test_Finite_Difference_Second_Derivative

Test the implementation of combinations of finite difference schemes so that all gridpoints have an evaluated first derivative up to required accuracy.

REASON ARCHIVED: The symmetric derivatives carried the wrong error estimate.
*/

#include <iostream>
#include <iomanip>
#include <vector>
#include <math.h>
#include <fstream>
#include <chrono>


//----- Declare functions -----
double f_function         ( double x );
double f_dx2_exact_function( double x );

std::vector<double> f_dx2_FD_function_order_1( std::vector<double> f, double h );
std::vector<double> f_dx2_FD_function_order_2( std::vector<double> f, double h );
std::vector<double> f_dx2_FD_function_order_3( std::vector<double> f, double h );
std::vector<double> f_dx2_FD_function_order_4( std::vector<double> f, double h );
std::vector<double> f_dx2_FD_function_order_5( std::vector<double> f, double h );



//----- Main code -----
int main(){
	
	
	//----- Extra variables useful only for this testing code -----
	
	double x_0     = 2.0;	// Minimum gridpoint.
	double h       = 0.1;	// Constant grid spacing. Probably doesn't matter but keep positive.
	int    N       = 10;	// Number of gridpoints to use for the entire domain. Keep low so that the endpoints are easily diagnosed.
	int    h_order = 5;		// The code will chose the scheme with accuracy O( h^N ). This is the N.
	
	std::string output_filename_test_base = "Test_finite_difference_second_derivative_20240501";	// Don't add ".csv" as the code will append _r, _t, etc to the filename for each csv.
	
	std::string output_filename_test_csv = "../CSV/"  + output_filename_test_base + "_" + std::to_string( h_order ) + ".csv";
	
	std::ofstream output_file_test_csv;
	
	std::vector<double> x                      ( N );
	std::vector<double> f                      ( N );
	std::vector<double> f_dx2_exact            ( N );
	std::vector<double> f_dx2_FD               ( N );
	std::vector<double> f_dx2_ARE              ( N );	// Absolute relative error of finite-difference derivative at each gridpoint.
	std::vector<double> f_dx2_abs_err          ( N );	// Absolute error.
	std::vector<double> f_dx2_abs_err_expected ( N );	// Expected absolute error ( equal to h^N ).
	
	
	//----- Build vectors of gridpoints and exact values -----
	for( int n=0; n<N; n++ ){
		x[n] = x_0 + n * h;
		f[n] = f_function( x[n] );
		f_dx2_exact[n] = f_dx2_exact_function( x[n] );
	}
	
	
	//----- Calculate derivatives to various accuracies -----
	std::cout << "Calculating second derivatives to order ( h^" << h_order << " )." << std::endl;
	
	switch( h_order ){
		
		case 1:
			f_dx2_FD = f_dx2_FD_function_order_1( f, h );
			break;
		case 2:
			f_dx2_FD = f_dx2_FD_function_order_2( f, h );
			break;
		case 3:
			f_dx2_FD = f_dx2_FD_function_order_3( f, h );
			break;
		case 4:
			f_dx2_FD = f_dx2_FD_function_order_4( f, h );
			break;
		case 5:
			f_dx2_FD = f_dx2_FD_function_order_5( f, h );
			break;
		default:
			std::cout << "No finite-difference method has been defined for order " << h_order << ". Stopping the code." << std::endl;
			return 1;
	}
	
	for( int n=0; n<N; n++ ){
		f_dx2_ARE             [n] = ( f_dx2_exact[n] != 0 ) ? ( abs( 1.0 - f_dx2_FD[n] / f_dx2_exact[n] ) ) : 0.0;
		f_dx2_abs_err         [n] = abs( f_dx2_FD[n] - f_dx2_exact[n] );
		f_dx2_abs_err_expected[n] = pow( h, h_order );
	}
	
	
	
	//----- Output to CSV -----
	output_file_test_csv.open( output_filename_test_csv );
	output_file_test_csv << "n,f_dx2_exact,f_dx2_FD,ARE,absolute_error,abs_err_expected\n";
	
	for( int n=0; n<N; n++ ){
		output_file_test_csv << n <<","<< f_dx2_exact[n] <<","<< f_dx2_FD[n] <<","<< f_dx2_ARE[n] <<","<< f_dx2_abs_err[n] <<","<< f_dx2_abs_err_expected[n] << "\n";
	}
	
	
	
	//----- Code finished -----
	
	std::cout << "\nCSV file saved:\t" << output_filename_test_csv << std::endl;
	
	std::cout << "\nCode finished" << std::endl;
	
	return 0;
	
}


//----- Define functions -----
double f_function         ( double x ){
	return pow( x, -3.0 );
}


double f_dx2_exact_function( double x ){
	return 12.0 * pow( x, -5.0 );
}




//----- Odd-order h -----
// Can have all datapoints covered, but only if we use forward and backward expressions for *all* points outside the range where the symmetric can be applied.

std::vector<double> f_dx2_FD_function_order_1( std::vector<double> f, double h ){
	int N = f.size();
	std::vector<double> f_dx2_approx ( N );
	
	f_dx2_approx[0] = ( f[2] - 2.0*f[1] + f[0] ) / ( h*h );
	
	for( int n=1; n<N-1; n++ ){
		f_dx2_approx[n] = ( f[n+1] - 2.0*f[n] + f[n-1] ) / ( h*h );
	}
	
	f_dx2_approx[N-1] = ( f[N-1] - 2.0*f[N-2] + f[N-3] ) / ( h*h );
	
	return f_dx2_approx;
}

std::vector<double> f_dx2_FD_function_order_3( std::vector<double> f, double h ){
	int N = f.size();
	std::vector<double> f_dx2_approx ( N );
	
	for( int n=0; n<=1; n++ ){
		f_dx2_approx[n] = ( 11.0*f[n+4] - 56.0*f[n+3] + 114.0*f[n+2] - 104.0*f[n+1] + 35.0*f[n] ) / ( 12.0*h*h );
	}
	
	for( int n=2; n<=N-3; n++ ){
		f_dx2_approx[n] = ( -f[n+2] + 16.0*f[n+1] - 30.0*f[n] + 16.0*f[n-1] - f[n-2] ) / ( 12.0*h*h );
	}
	
	for( int n=N-2; n<=N-1; n++ ){
		f_dx2_approx[n] = ( 11.0*f[n-4] - 56.0*f[n-3] + 114.0*f[n-2] - 104.0*f[n-1] + 35.0*f[n] ) / ( 12.0*h*h );
	}
	
	f_dx2_approx[N-1] = ( f[N-1] - 2.0*f[N-2] + f[N-3] ) / ( h*h );
	
	return f_dx2_approx;
}

std::vector<double> f_dx2_FD_function_order_5( std::vector<double> f, double h ){
	int N = f.size();
	std::vector<double> f_dx2_approx ( N );
	
	for( int n=0; n<=2; n++ ){
		f_dx2_approx[n] = ( 137.0*f[n+6] - 972.0*f[n+5] + 2970.0*f[n+4] - 5080.0*f[n+3] + 5265.0*f[n+2] - 3132.0*f[n+1] + 812.0*f[n] ) / ( 180.0*h*h );
	}
	
	for( int n=3; n<=N-4; n++ ){
		f_dx2_approx[n] = ( 2.0*f[n+3] - 27.0*f[n+2] + 270.0*f[n+1] - 490.0*f[n] + 270.0*f[n-1] - 27.0*f[n-2] + 2.0*f[n-3] ) / ( 180.0*h*h );
	}
	
	for( int n=N-3; n<=N-1; n++ ){
		f_dx2_approx[n] = ( 137.0*f[n-6] - 972.0*f[n-5] + 2970.0*f[n-4] - 5080.0*f[n-3] + 5265.0*f[n-2] - 3132.0*f[n-1] + 812.0*f[n] ) / ( 180.0*h*h );
	}
	
	return f_dx2_approx;
}




//----- Even-order h -----
/*
Can have all datapoints covered, but we can't get a symmetric expression.
We could use e.g. forward for all but the last few and backward for the last, but this probably isn't a good idea.
Instead, just do the two endpoints so that we can check that expression works, and don't use an even-order h second derivative method in the main code.
*/


std::vector<double> f_dx2_FD_function_order_2( std::vector<double> f, double h ){
	int N = f.size();
	std::vector<double> f_dx2_approx ( N );
	
	f_dx2_approx[0  ] = ( -f[3  ] + 4.0*f[2  ] - 5.0*f[1  ] + 2.0*f[0  ] ) / ( h*h );
	f_dx2_approx[N-1] = ( -f[N-4] + 4.0*f[N-3] - 5.0*f[N-2] + 2.0*f[N-1] ) / ( h*h );
	
	return f_dx2_approx;
}


std::vector<double> f_dx2_FD_function_order_4( std::vector<double> f, double h ){
	int N = f.size();
	std::vector<double> f_dx2_approx ( N );
	
	f_dx2_approx[0  ] = ( -10.0*f[5  ] + 61.0*f[4  ] - 156.0*f[3  ] + 214.0*f[2  ] - 154.0*f[1  ] + 45.0*f[0  ] ) / ( 12.0*h*h );
	f_dx2_approx[N-1] = ( -10.0*f[N-6] + 61.0*f[N-5] - 156.0*f[N-4] + 214.0*f[N-3] - 154.0*f[N-2] + 45.0*f[N-1] ) / ( 12.0*h*h );
	
	return f_dx2_approx;
}