/*
Trapezium_rule.h

Trapezium rule for numerical integration of a function. The function should be fed as a list of precalculated values, and the coordinates should also be fed as lists of values.

*/

#include <vector>

//----- One dimension -----
double integral_trapezium_1D_Cartesian_with_cout( std::vector<double> f, std::vector<double> x ){
	// 1D integral of f(x) given pre-calculated lists of coordinates and function values, by the trapezium rule.
	// In this version, we print all intermediate values to the screen and make it explicit how the summation index works.
	// This is useful for verifying if any indices are missed.
	
	int    N   = x.size();	// Number of datapoints
	double ret = 0;
	
	std::cout << "\ni\tterm\tresult" << std::endl;
	
	for( int i=0; i<=N-2; i++ ){
		// Use i<=N-2 not i<N-1 to be explicitly the same as the formula written in the report.
		double term = 0.5 * ( f[i+1] + f[i] ) * ( x[i+1] - x[i] );
		ret += term;
		std::cout << i <<"\t"<< term <<"\t"<< ret << std::endl;
	}
	
	return ret;

}


double integral_trapezium_1D_Cartesian( std::vector<double> f, std::vector<double> x ){
	// 1D integral of f(x) given pre-calculated lists of coordinates and function values, by the trapezium rule.
	
	int    N   = x.size();	// Number of datapoints
	double ret = 0;
	
	for( int i=0; i<=N-2; i++ ){
		ret += ( f[i+1] + f[i] ) * ( x[i+1] - x[i] );
	}
	
	return 0.5 * ret;

}


double integral_trapezium_1D_Cartesian_efficient( std::vector<double> f, std::vector<double> x ){
	// 1D integral of f(x) given pre-calculated lists of coordinates and function values, by the trapezium rule.
	
	int N = x.size();	// Number of datapoints
	
	double ret = f[0] * ( x[1] - x[0] ) + f[N-1] * ( x[N-1] - x[N-2] );
	
	for( int i=1; i<=N-2; i++ ){
		ret += f[i] * ( x[i+1] - x[i-1] );
	}
	
	return 0.5 * ret;
	
}


double integral_trapezium_1D_Cartesian_constant_spacing( std::vector<double> f, double h ){
	// 1D integral of f(x) given pre-calculated list of function values, with coordinate spacing h,	by the trapezium rule.
	
	int    N   = f.size();	// Number of datapoints
	double ret = 0;
	
	for( int i=0; i<=N-2; i++ ){
		ret += f[i+1] + f[i];
	}
	
	return 0.5 * h * ret;

}


double integral_trapezium_1D_Cartesian_constant_spacing_efficient( std::vector<double> f, double h ){
	// 1D integral of f(x) given pre-calculated list of function values, with coordinate spacing h,	by the trapezium rule.
	
	int N = f.size();	// Number of datapoints
	
	double ret = 0.5 * ( f[0] + f[N-1] );
	
	for( int i=1; i<=N-2; i++ ){
		ret += f[i];
	}
	
	return h * ret;
	
}




double integral_trapezium_1D( std::vector<double> f, std::vector<double> du ){
	// 1D integral of f(u) given pre-calculated list of function values, in a general coordinate system with line elements du.
	
	int    N   = f.size();	// Number of datapoints
	double ret = 0;
	
	for( int i=0; i<=N-2; i++ ){
		ret += ( f[i+1] + f[i] ) * du[i];
	}
	
	return 0.5 * ret;

}




//----- Two dimensions -----
double integral_trapezium_2D_with_cout( std::vector< std::vector<double> > f, std::vector< std::vector<double> > dA ){
	
	int N1 = f   .size();
	int N2 = f[0].size();
	
	std::cout << "\nN1 = " << N1 <<"\tN2 = " << N2 << std::endl;
	std::cout << "\ni\tj\tterm\tresult" << std::endl;
	
	double ret = 0;
	
	for( int i=0; i<=N1-2; i++ ){
		for( int j=0; j<=N2-2; j++ ){
			
			double term = 0.25 * ( f[i][j] + f[i][j+1] + f[i+1][j] + f[i+1][j+1] ) * dA[i][j];
			ret += term;
			std::cout << i <<"\t"<< j <<"\t"<< term <<"\t" << ret << std::endl;
			
		}
	}
	
	return ret;
	
}


double integral_trapezium_2D( std::vector< std::vector<double> > f, std::vector< std::vector<double> > dA ){
	
	int N1 = f   .size();
	int N2 = f[0].size();
	
	double ret = 0;
	
	for( int i=0; i<=N1-2; i++ ){
		for( int j=0; j<=N2-2; j++ ){
			ret += ( f[i][j] + f[i][j+1] + f[i+1][j] + f[i+1][j+1] ) * dA[i][j];
		}
	}
	
	return 0.25 * ret;
	
}




//----- Three dimensions -----
double integral_trapezium_3D( std::vector< std::vector< std::vector<double> > > f, std::vector< std::vector< std::vector<double> > > dV ){
	
	int N1 = f      .size();
	int N2 = f[0]   .size();
	int N3 = f[0][0].size();
	
	double ret = 0;
	
	for( int i=0; i<=N1-2; i++ ){
		for( int j=0; j<=N2-2; j++ ){
			for( int k=0; k<=N3-2; k++ ){
				ret += ( f[i][j][k] + f[i][j][k+1] + f[i][j+1][k] + f[i][j+1][k+1] + f[i+1][j][k] + f[i+1][j][k+1] + f[i+1][j+1][k] + f[i+1][j+1][k+1] ) * dV[i][j][k];
			}
		}
	}
	
	return 0.125 * ret;
	
}




//----- One-dimensional Newton-Cotes methods -----
double integral_1D_Newton_Cotes_2( std::vector<double> f, double h ){
	// Newton-Cotes method for 2-gridpoint subinterval. Equivalent to the trapezium rule.
	
	int N = f.size();
	
	double ret = 0;
	
	for( int i=0; i<=N-2; i++ ){
		ret += 0.5 * ( f[i] + f[i+1] );
	}
	
	return ret * h;
	
}


double integral_1D_Newton_Cotes_3( std::vector<double> f, double h ){
	// Newton-Cotes method for 3-gridpoint subinterval. Equivalent to Simpson's rule. Use trapezium rule at endpoints.
	// DOESN'T YET WORK 20240521.
	std::cout << "WARNING: This function isn't yet working as of 20240521." << std::endl;
	
	int N = f.size();
	
	double ret = 0;
	
	ret += 0.5 * ( f[0  ] + f[1  ] );
	ret += 0.5 * ( f[N-2] + f[N-1] );
	
	for( int i=1; i<=N-3; i+=2 ){
		ret += (1.0/3.0) * ( f[i-1] + f[i+1] + 4.0*f[i] );
	}
	
	return ret * h;
	
}


double integral_1D_Newton_Cotes_4( std::vector<double> f, double h ){
	// Newton-Cotes method for 4-gridpoint subinterval. Equivalent to Simpson's 3/8 rule. Only works if number of data points is one more than a multiple of 3.
	http://www.holoborodko.com/pavel/numerical-methods/numerical-integration/overlapped-newton-cotes-quadratures/
	
	int N = f.size();
	
	if( N % 3 != 1 ){
		std::cout << "This method is most accurate when you have 3m+1 data points. Consider a different method or adding one or two extra points." << std::endl;
	}
	
	double ret = 0;
	
	for( int k=1; k<=N/3; k++ ){
		ret += f[3*k-3] + f[3*k] + 3.0 * ( f[3*k-2] + f[3*k-1] );
	}
	
	return ret * h * 0.375;
	
}

double integral_1D_Newton_Cotes_5( std::vector<double> f, double h ){
	// Newton-Cotes method for 5-gridpoint subinterval. This is equivalent to Boole's rule. Only works if number of data points is one more than a multiple of 4.
	http://www.holoborodko.com/pavel/numerical-methods/numerical-integration/overlapped-newton-cotes-quadratures/
	
	int N = f.size();
	
	if( N % 4 != 1 ){
		std::cout << "This method is most accurate when you have 5m+1 data points. Consider a different method or adding one or two extra points." << std::endl;
	}
	
	double ret = 0;
	
	for( int k=1; k<=N/4; k++ ){
		ret += 7.0 * ( f[4*k-4] + f[4*k] ) + 32.0 * ( f[4*k-3] + f[4*k-1] ) + 12.0 * f[4*k-2];
	}
	
	return ret * h * 2.0 / 45.0;
	
}




double integral_1D_Newton_Cotes_6( std::vector<double> f, double h ){
	// Newton-Cotes method for 6-gridpoint subinterval. Only works if number of data points is one more than a multiple of 5.
	http://www.holoborodko.com/pavel/numerical-methods/numerical-integration/overlapped-newton-cotes-quadratures/
	
	int N = f.size();
	
	if( N % 5 != 1 ){
		std::cout << "This method is most accurate when you have 9m+1 data points. Consider a different method or adding one or two extra points." << std::endl;
	}
	
	double ret = 0;
	
	for( int k=1; k<=N/5; k++ ){
		ret += 19.0 * ( f[5*k-5] + f[5*k] ) + 75.0 * ( f[5*k-4] + f[5*k-1] ) + 50.0 * ( f[5*k-3] + f[5*k-2] );
	}
	
	return ret * h * 5.0 / 288.0;
	
}