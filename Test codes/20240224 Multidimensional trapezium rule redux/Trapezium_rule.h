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




//----- Two dimensions -----
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