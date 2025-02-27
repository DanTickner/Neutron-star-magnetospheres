/*
Finite_Difference_Expressions.h

First and second derivatives of a list of function values h, with constant gridpoint spacing h.
Expressions given up to various accuracies, which use various numbers of gridpoints.
Keep in header file to avoid copying between codes, which is vulnerable to divergence between expressions that are not updated.
Once a final expression decided, copy it into the main evolution code and use it for all derivatives.

*/



//----- First derivatives -----
std::vector<double> f_dx_FD_function_order_1( std::vector<double> &f, double h ){
	
	int N = f.size();
	std::vector<double> f_dx ( N );
	
	for( int n=0; n<N-1; n++ ){
		f_dx[n] = ( f[n+1] - f[n] ) / h;
	}
	
	f_dx[N-1] = ( f[N-1] - f[N-2] ) / h;
	
	return f_dx;
}


std::vector<double> f_dx_FD_function_order_2( std::vector<double> &f, double h ){
	
	int N = f.size();
	std::vector<double> f_dx ( N );
	
	f_dx[0  ] = ( -f[2  ] + 4.0*f[1  ] - 3.0*f[0  ] ) / ( 2.0*h );
	f_dx[N-1] = (  f[N-3] - 4.0*f[N-2] + 3.0*f[N-1] ) / ( 2.0*h );
	
	for( int n=1; n<=N-2; n++ ){
		f_dx[n] = ( f[n+1] - f[n-1] ) / ( 2.0*h );
	}
	
	return f_dx;
}


std::vector<double> f_dx_FD_function_order_4( std::vector<double> &f, double h ){
	
	int N = f.size();
	std::vector<double> f_dx ( N );
	
	f_dx[0  ] = ( -3.0*f[4  ] + 16.0*f[3  ] - 36.0*f[2  ] + 48.0*f[1  ] - 25.0*f[0  ] ) / ( 12.0*h );
	f_dx[N-1] = (  3.0*f[N-5] - 16.0*f[N-4] + 36.0*f[N-3] - 48.0*f[N-2] + 25.0*f[N-1] ) / ( 12.0*h );
	
	f_dx[1  ] = (      f[4  ] -  6.0*f[3  ] + 18.0*f[2  ] - 10.0*f[1  ] -  3.0*f[0  ] ) / ( 12.0*h );
	f_dx[N-2] = ( -    f[N-5] +  6.0*f[N-4] - 18.0*f[N-3] + 10.0*f[N-2] +  3.0*f[N-1] ) / ( 12.0*h );
	
	for( int n=2; n<=N-3; n++ ){
		f_dx[n] = ( -f[n+2]+f[n-2] + 8.0*(f[n+1]-f[n-1])  ) / ( 12.0*h );
	}
	
	return f_dx;
}

std::vector<double> f_dx_FD_function_order_6( std::vector<double> &f, double h ){
	
	int N = f.size();
	std::vector<double> f_dx ( N );
	
	f_dx[0  ] = ( -10.0*f[6  ] + 72.0*f[5  ] - 225.0*f[4  ] + 400.0*f[3  ] - 450.0*f[2  ] + 360.0*f[1  ] - 147.0*f[0  ] ) / ( 60.0*h );
	f_dx[N-1] = (  10.0*f[N-7] - 72.0*f[N-6] + 225.0*f[N-5] - 400.0*f[N-4] + 450.0*f[N-3] - 360.0*f[N-2] + 147.0*f[N-1] ) / ( 60.0*h );
	
	f_dx[1  ] = (   2.0*f[6  ] - 15.0*f[5  ] +  50.0*f[4  ] - 100.0*f[3  ] + 150.0*f[2  ] -  77.0*f[1  ] -  10.0*f[0  ] ) / ( 60.0*h );
	f_dx[N-2] = ( - 2.0*f[N-7] + 15.0*f[N-6] -  50.0*f[N-5] + 100.0*f[N-4] - 150.0*f[N-3] +  77.0*f[N-2] +  10.0*f[N-1] ) / ( 60.0*h );
	
	f_dx[2  ] = (      -f[6  ] +  8.0*f[5  ] -  30.0*f[4  ] +  80.0*f[3  ] -  35.0*f[2  ] -  24.0*f[1  ] +   2.0*f[0  ] ) / ( 60.0*h );
	f_dx[N-3] = (       f[N-7] -  8.0*f[N-6] +  30.0*f[N-5] -  80.0*f[N-4] +  35.0*f[N-3] +  24.0*f[N-2] -   2.0*f[N-1] ) / ( 60.0*h );
	
	for( int n=3; n<=N-4; n++){
		f_dx[n] = ( f[n+3]-f[n-3] - 9.0*(f[n+2]-f[n-2]) + 45.0*(f[n+1]-f[n-1]) ) / ( 60.0*h );
	}
	
	return f_dx;
}


std::vector<double> f_dx_FD_function_order_6_with_8_near_endpoints( std::vector<double> &f, double h ){
	
	int N = f.size();
	std::vector<double> f_dx ( N );
	
	f_dx[0  ] = ( -105.0*f[8  ] + 960.0*f[7  ] - 3920.0*f[6  ] + 9408.0*f[5  ] - 14700.0*f[4  ] + 15680.0*f[3  ] - 11760.0*f[2  ] + 6720.0*f[1  ] - 2283.0*f[0  ] ) / ( 840.0*h );
	f_dx[N-1] = (  105.0*f[N-9] - 960.0*f[N-8] + 3920.0*f[N-7] - 9408.0*f[N-6] + 14700.0*f[N-5] - 15680.0*f[N-4] + 11760.0*f[N-3] - 6720.0*f[N-2] + 2283.0*f[N-1] ) / ( 840.0*h );
	
	f_dx[1  ] = (   15.0*f[8  ] - 140.0*f[7  ] +  588.0*f[6  ] - 1470.0*f[5  ] +  2450.0*f[4  ] -  2940.0*f[3  ] +  2940.0*f[2  ] - 1338.0*f[1  ] -  105.0*f[0  ] ) / ( 840.0*h );
	f_dx[N-2] = (  -15.0*f[N-9] + 140.0*f[N-8] -  588.0*f[N-7] + 1470.0*f[N-6] -  2450.0*f[N-5] +  2940.0*f[N-4] -  2940.0*f[N-3] + 1338.0*f[N-2] +  105.0*f[N-1] ) / ( 840.0*h );
	
	f_dx[2  ] = (   -5.0*f[8  ] +  48.0*f[7  ] -  210.0*f[6  ] +  560.0*f[5  ] -  1050.0*f[4  ] +  1680.0*f[3  ] -   798.0*f[2  ] -  240.0*f[1  ] +   15.0*f[0  ] ) / ( 840.0*h );
	f_dx[N-3] = (    5.0*f[N-9] -  48.0*f[N-8] +  210.0*f[N-7] -  560.0*f[N-6] +  1050.0*f[N-5] -  1680.0*f[N-4] +   798.0*f[N-3] +  240.0*f[N-2] -   15.0*f[N-1] ) / ( 840.0*h );
	
	for( int n=3; n<=N-4; n++){
		f_dx[n] = ( f[n+3]-f[n-3] - 9.0*(f[n+2]-f[n-2]) + 45.0*(f[n+1]-f[n-1]) ) / ( 60.0*h );
	}
	
	return f_dx;
	
}

std::vector<double> f_dx_FD_function_order_10( std::vector<double> &f, double h ){
	
	int N = f.size();
	std::vector<double> f_dx ( N );
	
	f_dx[0  ] = ( -252.0*f[10  ] + 2800.0*f[9   ] - 14175.0*f[8  ] + 43200.0*f[7  ] - 88200.0*f[6  ] + 127008.0*f[5  ] - 132300.0*f[4  ] + 100800.0*f[3  ] - 56700.0*f[2  ] + 25200.0*f[1  ] - 7381.0*f[0  ] ) / ( 2520.0*h );
	f_dx[N-1] = (  252.0*f[N-11] - 2800.0*f[N-10] + 14175.0*f[N-9] - 43200.0*f[N-8] + 88200.0*f[N-7] - 127008.0*f[N-6] + 132300.0*f[N-5] - 100800.0*f[N-4] + 56700.0*f[N-3] - 25200.0*f[N-2] + 7381.0*f[N-1] ) / ( 2520.0*h );
	
	f_dx[1  ] = (   28.0*f[10  ] -  315.0*f[9   ] +  1620.0*f[8  ] -  5040.0*f[7  ] + 10584.0*f[6  ] -  15876.0*f[5  ] +  17640.0*f[4  ] -  15120.0*f[3  ] + 11340.0*f[2  ] -  4609.0*f[1  ] -  252.0*f[0  ] ) / ( 2520.0*h );
	f_dx[N-2] = (  -28.0*f[N-11] +  315.0*f[N-10] -  1620.0*f[N-9] +  5040.0*f[N-8] - 10584.0*f[N-7] +  15876.0*f[N-6] -  17640.0*f[N-5] +  15120.0*f[N-4] - 11340.0*f[N-3] +  4609.0*f[N-2] +  252.0*f[N-1] ) / ( 2520.0*h );
	
	f_dx[2  ] = (   -7.0*f[10  ] +   80.0*f[9   ]  -  420.0*f[8  ] +  1344.0*f[7  ] -  2940.0*f[6  ] +   4704.0*f[5  ] -   5880.0*f[4  ] +   6720.0*f[3  ] -  3069.0*f[2  ] -   560.0*f[1  ] +   28.0*f[0  ] ) / ( 2520.0*h );
	f_dx[N-3] = (    7.0*f[N-11] -   80.0*f[N-10]  +  420.0*f[N-9] -  1344.0*f[N-8] +  2940.0*f[N-7] -   4704.0*f[N-6] +   5880.0*f[N-5] -   6720.0*f[N-4] +  3069.0*f[N-3] +   560.0*f[N-2] -   28.0*f[N-1] ) / ( 2520.0*h );
	
	f_dx[3  ] = (    3.0*f[10  ] -   35.0*f[9   ] +   189.0*f[8  ] -   630.0*f[7  ] +  1470.0*f[6  ] -   2646.0*f[5  ] +   4410.0*f[4  ] -   1914.0*f[3  ] -   945.0*f[2  ] +   105.0*f[1  ] -    7.0*f[0  ] ) / ( 2520.0*h );
	f_dx[N-4] = (   -3.0*f[N-11] +   35.0*f[N-10] -   189.0*f[N-9] +   630.0*f[N-8] -  1470.0*f[N-7] +   2646.0*f[N-6] -   4410.0*f[N-5] +   1914.0*f[N-4] +   945.0*f[N-3] -   105.0*f[N-2] +    7.0*f[N-1] ) / ( 2520.0*h );
	
	f_dx[4  ] = (   -2.0*f[10  ] +   24.0*f[9   ] -   135.0*f[8  ] +   480.0*f[7  ] -  1260.0*f[6  ] +   3024.0*f[5  ] -    924.0*f[4  ] -   1440.0*f[3  ] +   270.0*f[2  ] -    40.0*f[1  ] +    3.0*f[0  ] ) / ( 2520.0*h );
	f_dx[N-5] = (    2.0*f[N-11] -   24.0*f[N-10] +   135.0*f[N-9] -   480.0*f[N-8] +  1260.0*f[N-7] -   3024.0*f[N-6] +    924.0*f[N-5] +   1440.0*f[N-4] -   270.0*f[N-3] +    40.0*f[N-2] -    3.0*f[N-1] ) / ( 2520.0*h );
	
	for( int n=5; n<=N-6; n++ ){
		f_dx[n] = ( 2.0*(f[n+5]-f[n-5]) - 25.0*(f[n+4]-f[n-4]) + 150.0*(f[n+3]-f[n-3]) - 600.0*(f[n+2]-f[n-2]) + 2100.0*(f[n+1]-f[n-1]) ) / ( 2520.0*h );
	}
	
	return f_dx;
	
}





//----- Second derivatives -----
//--- Odd-order h ---
/*
Can have all datapoints covered, but we can't get a symmetric expression.
We could use e.g. forward for all but the last few and backward for the last, but this probably isn't a good idea.
Instead, just do the two end regions so that we can check that expression works, and don't use an even-order h second derivative method in the main code.
*/

std::vector<double> f_dx2_FD_function_order_1( std::vector<double> &f, double h ){
	
	int N = f.size();
	std::vector<double> f_dx2 ( N );
	
	f_dx2[0  ] = ( f[2  ] - 2.0*f[1  ] + f[0  ] ) / ( h*h );
	f_dx2[N-1] = ( f[N-1] - 2.0*f[N-2] + f[N-3] ) / ( h*h );
	
	return f_dx2;
}

std::vector<double> f_dx2_FD_function_order_3( std::vector<double> &f, double h ){
	
	int N = f.size();
	std::vector<double> f_dx2 ( N );
	
	for( int n=0; n<=1; n++ ){
		f_dx2[n] = ( 11.0*f[n+4] - 56.0*f[n+3] + 114.0*f[n+2] - 104.0*f[n+1] + 35.0*f[n] ) / ( 12.0*h*h );
	}
	
	for( int n=N-2; n<=N-1; n++ ){
		f_dx2[n] = ( 11.0*f[n-4] - 56.0*f[n-3] + 114.0*f[n-2] - 104.0*f[n-1] + 35.0*f[n] ) / ( 12.0*h*h );
	}
	
	f_dx2[N-1] = ( f[N-1] - 2.0*f[N-2] + f[N-3] ) / ( h*h );
	
	return f_dx2;
}

std::vector<double> f_dx2_FD_function_order_5( std::vector<double> &f, double h ){
	
	int N = f.size();
	std::vector<double> f_dx2 ( N );
	
	for( int n=0; n<=2; n++ ){
		f_dx2[n] = ( 137.0*f[n+6] - 972.0*f[n+5] + 2970.0*f[n+4] - 5080.0*f[n+3] + 5265.0*f[n+2] - 3132.0*f[n+1] + 812.0*f[n] ) / ( 180.0*h*h );
	}
	
	for( int n=3; n<=N-4; n++ ){
		f_dx2[n] = ( 2.0*f[n+3] - 27.0*f[n+2] + 270.0*f[n+1] - 490.0*f[n] + 270.0*f[n-1] - 27.0*f[n-2] + 2.0*f[n-3] ) / ( 180.0*h*h );	// CHANGED DENOMINATOR FROM 180 ON 20240501
	}
	
	for( int n=N-3; n<=N-1; n++ ){
		f_dx2[n] = ( 137.0*f[n-6] - 972.0*f[n-5] + 2970.0*f[n-4] - 5080.0*f[n-3] + 5265.0*f[n-2] - 3132.0*f[n-1] + 812.0*f[n] ) / ( 180.0*h*h );
	}
	
	return f_dx2;
}




//--- Even-order h ---
// Can have all datapoints covered, but only if we use forward and backward expressions for *all* points outside the range where the symmetric can be applied.


std::vector<double> f_dx2_FD_function_order_2( std::vector<double> &f, double h ){
	
	int N = f.size();
	std::vector<double> f_dx2 ( N );
	
	f_dx2[0  ] = ( -f[3  ] + 4.0*f[2  ] - 5.0*f[1  ] + 2.0*f[0  ] ) / ( h*h );
	f_dx2[N-1] = ( -f[N-4] + 4.0*f[N-3] - 5.0*f[N-2] + 2.0*f[N-1] ) / ( h*h );
	
	for( int n=1; n<=N-2; n++ ){
		f_dx2[n] = ( f[n+1] - 2.0*f[n] + f[n-1] ) / ( h*h );
	}
	
	return f_dx2;
}


std::vector<double> f_dx2_FD_function_order_4( std::vector<double> &f, double h ){
	
	int N = f.size();
	std::vector<double> f_dx2 ( N );
	
	for( int n=0; n<=1; n++ ){
		f_dx2[n] = ( -10.0*f[n+5] + 61.0*f[n+4] - 156.0*f[n+3] + 214.0*f[n+2] - 154.0*f[n+1] + 45.0*f[n] ) / ( 12.0*h*h );
	}
	
	for( int n=N-2; n<=N-1; n++ ){
		f_dx2[n] = ( -10.0*f[n-5] + 61.0*f[n-4] - 156.0*f[n-3] + 214.0*f[n-2] - 154.0*f[n-1] + 45.0*f[n] ) / ( 12.0*h*h );
	}
	
	for( int n=2; n<=N-3; n++ ){
		f_dx2[n] = ( -(f[n+2]+f[n-2]) + 16.0*(f[n+1]+f[n-1]) - 30.0*f[n] ) / ( 12.0*h*h );
	}
	
	return f_dx2;
}


std::vector<double> f_dx2_FD_function_order_6( std::vector<double> &f, double h ){
	
	int N = f.size();
	std::vector<double> f_dx2 ( N );
	
	for( int n=0; n<=2; n++ ){
		f_dx2[n] = ( -126.0*f[n+7] + 1019.0*f[n+6] - 3618.0*f[n+5] + 7380.0*f[n+4] - 9490.0*f[n+3] + 7911.0*f[n+2] - 4014.0*f[n+1] + 938.0*f[n] ) / ( 180.0*h*h );
	}
	
	for( int n=N-3; n<=N-1; n++ ){
		f_dx2[n] = ( -126.0*f[n-7] + 1019.0*f[n-6] - 3618.0*f[n-5] + 7380.0*f[n-4] - 9490.0*f[n-3] + 7911.0*f[n-2] - 4014.0*f[n-1] + 938.0*f[n] ) / ( 180.0*h*h );
	}
	
	for( int n=3; n<=N-4; n++ ){
		f_dx2[n] = ( 2.0*(f[n+3]+f[n-3]) - 27.0*(f[n+2]+f[n-2]) + 270.0*(f[n+1]+f[n-1]) - 490.0*f[n] ) / ( 180.0*h*h );
	}
	
	return f_dx2;
}


std::vector<double> f_dx2_FD_function_order_6_with_8_near_endpoints( std::vector<double> &f, double h ){
	
	int N = f.size();
	std::vector<double> f_dx2 ( N );
	
	for( int n=0; n<=3; n++ ){
		f_dx2[n] = ( 3267.0*f[n+8] - 29664.0*f[n+7] + 120008.0*f[n+6] - 284256.0*f[n+5] + 435330.0*f[n+4] - 448672.0*f[n+3] + 312984.0*f[n+2] - 138528.0*f[n+1] + 29531.0*f[n] ) / ( 5040.0*h*h );
	}
	
	for( int n=N-4; n<=N-1; n++ ){
		f_dx2[n] = ( 3267.0*f[n-8] - 29664.0*f[n-7] + 120008.0*f[n-6] - 284256.0*f[n-5] + 435330.0*f[n-4] - 448672.0*f[n-3] + 312984.0*f[n-2] - 138528.0*f[n-1] + 29531.0*f[n] ) / ( 5040.0*h*h );
	}
	
	for( int n=4; n<=N-5; n++ ){
		f_dx2[n] = ( 2.0*(f[n+3]+f[n-3]) - 27.0*(f[n+2]+f[n-2]) + 270.0*(f[n+1]+f[n-1]) - 490.0*f[n] ) / ( 180.0*h*h );
	}
	
	return f_dx2;
}





//----- First partial derivatives of 2D functions -----
std::vector< std::vector<double> > f_dx_FD_function_order_1( std::vector< std::vector<double> > &f, double h_x ){
	// Partial derivative wrt x of a function f(x,y).
	
	int N_x = f   .size();
	int N_y = f[0].size();
	std::vector< std::vector<double> > f_dx ( N_x, std::vector<double> ( N_y ) );
	
	for( int j=0; j<N_y; j++ ){
		
		for( int i=0; i<N_x-1; i++ ){
			f_dx[i][j] = ( f[i+1][j] - f[i][j] ) / h_x;
		}
		
		f_dx[N_x-1][j] = ( f[N_x-1][j] - f[N_x-2][j] ) / h_x;
	}
	
	return f_dx;
	
}

std::vector< std::vector<double> > f_dy_FD_function_order_1( std::vector< std::vector<double> > &f, double h_y ){
	// Partial derivative wrt y of a function f(x,y).
	
	int N_x = f   .size();
	int N_y = f[0].size();
	std::vector< std::vector<double> > f_dy ( N_x, std::vector<double> ( N_y ) );
	
	for( int i=0; i<N_x; i++ ){
		
		for( int j=0; j<N_y-1; j++ ){
			f_dy[i][j] = ( f[i][j+1] - f[i][j] ) / h_y;
		}
		
		f_dy[i][N_y-1] = ( f[i][N_y-1] - f[i][N_y-2] ) / h_y;
	}
	
	return f_dy;
	
}

std::vector< std::vector<double> > f_dx_FD_function_order_6_with_8_near_endpoints( std::vector< std::vector<double> > &f, double h_x ){
	
	int N_x = f   .size();
	int N_y = f[0].size();
	std::vector< std::vector<double> > f_dx ( N_x, std::vector<double> ( N_y ) );
	
	for( int j=0; j<N_y; j++ ){
	
		f_dx[0    ][j] = ( -105.0*f[8    ][j] + 960.0*f[7    ][j] - 3920.0*f[6    ][j] + 9408.0*f[5    ][j] - 14700.0*f[4    ][j] + 15680.0*f[3    ][j] - 11760.0*f[2    ][j] + 6720.0*f[1    ][j] - 2283.0*f[0    ][j] ) / ( 840.0*h_x );
		f_dx[N_x-1][j] = (  105.0*f[N_x-9][j] - 960.0*f[N_x-8][j] + 3920.0*f[N_x-7][j] - 9408.0*f[N_x-6][j] + 14700.0*f[N_x-5][j] - 15680.0*f[N_x-4][j] + 11760.0*f[N_x-3][j] - 6720.0*f[N_x-2][j] + 2283.0*f[N_x-1][j] ) / ( 840.0*h_x );
		
		f_dx[1    ][j] = (   15.0*f[8    ][j] - 140.0*f[7    ][j] +  588.0*f[6    ][j] - 1470.0*f[5    ][j] +  2450.0*f[4    ][j] -  2940.0*f[3    ][j] +  2940.0*f[2    ][j] - 1338.0*f[1    ][j] -  105.0*f[0    ][j] ) / ( 840.0*h_x );
		f_dx[N_x-2][j] = (  -15.0*f[N_x-9][j] + 140.0*f[N_x-8][j] -  588.0*f[N_x-7][j] + 1470.0*f[N_x-6][j] -  2450.0*f[N_x-5][j] +  2940.0*f[N_x-4][j] -  2940.0*f[N_x-3][j] + 1338.0*f[N_x-2][j] +  105.0*f[N_x-1][j] ) / ( 840.0*h_x );
		
		f_dx[2    ][j] = (   -5.0*f[8    ][j] +  48.0*f[7    ][j] -  210.0*f[6    ][j] +  560.0*f[5    ][j] -  1050.0*f[4    ][j] +  1680.0*f[3    ][j] -   798.0*f[2    ][j] -  240.0*f[1    ][j] +   15.0*f[0    ][j] ) / ( 840.0*h_x );
		f_dx[N_x-3][j] = (    5.0*f[N_x-9][j] -  48.0*f[N_x-8][j] +  210.0*f[N_x-7][j] -  560.0*f[N_x-6][j] +  1050.0*f[N_x-5][j] -  1680.0*f[N_x-4][j] +   798.0*f[N_x-3][j] +  240.0*f[N_x-2][j] -   15.0*f[N_x-1][j] ) / ( 840.0*h_x );
		
		for( int i=3; i<=N_x-4; i++){
			f_dx[i][j] = ( f[i+3][j]-f[i-3][j] - 9.0*(f[i+2][j]-f[i-2][j]) + 45.0*(f[i+1][j]-f[i-1][j]) ) / ( 60.0*h_x );
		}
	
	}
	
	return f_dx;
	
}

std::vector< std::vector<double> > f_dy_FD_function_order_6_with_8_near_endpoints( std::vector< std::vector<double> > &f, double h_y ){
	
	int N_x = f   .size();
	int N_y = f[0].size();
	std::vector< std::vector<double> > f_dy ( N_x, std::vector<double> ( N_y ) );
	
	for( int i=0; i<N_x; i++ ){
	
		f_dy[i][0    ] = ( -105.0*f[i][8    ] + 960.0*f[i][7    ] - 3920.0*f[i][6    ] + 9408.0*f[i][5    ] - 14700.0*f[i][4    ] + 15680.0*f[i][3    ] - 11760.0*f[i][2    ] + 6720.0*f[i][1    ] - 2283.0*f[i][0    ] ) / ( 840.0*h_y );
		f_dy[i][N_y-1] = (  105.0*f[i][N_y-9] - 960.0*f[i][N_y-8] + 3920.0*f[i][N_y-7] - 9408.0*f[i][N_y-6] + 14700.0*f[i][N_y-5] - 15680.0*f[i][N_y-4] + 11760.0*f[i][N_y-3] - 6720.0*f[i][N_y-2] + 2283.0*f[i][N_y-1] ) / ( 840.0*h_y );
		
		f_dy[i][1    ] = (   15.0*f[i][8    ] - 140.0*f[i][7    ] +  588.0*f[i][6    ] - 1470.0*f[i][5    ] +  2450.0*f[i][4    ] -  2940.0*f[i][3    ] +  2940.0*f[i][2    ] - 1338.0*f[i][1    ] -  105.0*f[i][0    ] ) / ( 840.0*h_y );
		f_dy[i][N_y-2] = (  -15.0*f[i][N_y-9] + 140.0*f[i][N_y-8] -  588.0*f[i][N_y-7] + 1470.0*f[i][N_y-6] -  2450.0*f[i][N_y-5] +  2940.0*f[i][N_y-4] -  2940.0*f[i][N_y-3] + 1338.0*f[i][N_y-2] +  105.0*f[i][N_y-1] ) / ( 840.0*h_y );
		
		f_dy[i][2    ] = (   -5.0*f[i][8    ] +  48.0*f[i][7    ] -  210.0*f[i][6    ] +  560.0*f[i][5    ] -  1050.0*f[i][4    ] +  1680.0*f[i][3    ] -   798.0*f[i][2    ] -  240.0*f[i][1    ] +   15.0*f[i][0    ] ) / ( 840.0*h_y );
		f_dy[i][N_y-3] = (    5.0*f[i][N_y-9] -  48.0*f[i][N_y-8] +  210.0*f[i][N_y-7] -  560.0*f[i][N_y-6] +  1050.0*f[i][N_y-5] -  1680.0*f[i][N_y-4] +   798.0*f[i][N_y-3] +  240.0*f[i][N_y-2] -   15.0*f[i][N_y-1] ) / ( 840.0*h_y );
		
		for( int j=3; j<=N_y-4; j++){
			f_dy[i][j] = ( f[i][j+3]-f[i][j-3] - 9.0*(f[i][j+2]-f[i][j-2]) + 45.0*(f[i][j+1]-f[i][j-1]) ) / ( 60.0*h_y );
		}
	
	}
	
	return f_dy;
	
}

std::vector< std::vector<double> > f_dx_FD_function_order_10( std::vector< std::vector<double> > &f, double h_x ){
	
	int N_x = f   .size();
	int N_y = f[0].size();
	std::vector< std::vector<double> > f_dx ( N_x, std::vector<double> ( N_y ) );
	
	for( int j=0; j<N_y; j++ ){
	
		f_dx[0    ][j] = ( -252.0*f[10    ][j] + 2800.0*f[9     ][j] - 14175.0*f[8    ][j] + 43200.0*f[7    ][j] - 88200.0*f[6    ][j] + 127008.0*f[5    ][j] - 132300.0*f[4    ][j] + 100800.0*f[3    ][j] - 56700.0*f[2    ][j] + 25200.0*f[1    ][j] - 7381.0*f[0    ][j] ) / ( 2520.0*h_x );
		f_dx[N_x-1][j] = (  252.0*f[N_x-11][j] - 2800.0*f[N_x-10][j] + 14175.0*f[N_x-9][j] - 43200.0*f[N_x-8][j] + 88200.0*f[N_x-7][j] - 127008.0*f[N_x-6][j] + 132300.0*f[N_x-5][j] - 100800.0*f[N_x-4][j] + 56700.0*f[N_x-3][j] - 25200.0*f[N_x-2][j] + 7381.0*f[N_x-1][j] ) / ( 2520.0*h_x );
		
		f_dx[1    ][j] = (   28.0*f[10    ][j] -  315.0*f[9     ][j] +  1620.0*f[8    ][j] -  5040.0*f[7    ][j] + 10584.0*f[6    ][j] -  15876.0*f[5    ][j] +  17640.0*f[4    ][j] -  15120.0*f[3    ][j] + 11340.0*f[2    ][j] -  4609.0*f[1    ][j] -  252.0*f[0    ][j] ) / ( 2520.0*h_x );
		f_dx[N_x-2][j] = (  -28.0*f[N_x-11][j] +  315.0*f[N_x-10][j] -  1620.0*f[N_x-9][j] +  5040.0*f[N_x-8][j] - 10584.0*f[N_x-7][j] +  15876.0*f[N_x-6][j] -  17640.0*f[N_x-5][j] +  15120.0*f[N_x-4][j] - 11340.0*f[N_x-3][j] +  4609.0*f[N_x-2][j] +  252.0*f[N_x-1][j] ) / ( 2520.0*h_x );
		
		f_dx[2    ][j] = (   -7.0*f[10    ][j] +   80.0*f[9     ][j] -   420.0*f[8    ][j] +  1344.0*f[7    ][j] -  2940.0*f[6    ][j] +   4704.0*f[5    ][j] -   5880.0*f[4    ][j] +   6720.0*f[3    ][j] -  3069.0*f[2    ][j] -   560.0*f[1    ][j] +   28.0*f[0    ][j] ) / ( 2520.0*h_x );
		f_dx[N_x-3][j] = (    7.0*f[N_x-11][j] -   80.0*f[N_x-10][j] +   420.0*f[N_x-9][j] -  1344.0*f[N_x-8][j] +  2940.0*f[N_x-7][j] -   4704.0*f[N_x-6][j] +   5880.0*f[N_x-5][j] -   6720.0*f[N_x-4][j] +  3069.0*f[N_x-3][j] +   560.0*f[N_x-2][j] -   28.0*f[N_x-1][j] ) / ( 2520.0*h_x );
		
		f_dx[3    ][j] = (    3.0*f[10    ][j] -   35.0*f[9     ][j] +   189.0*f[8    ][j] -   630.0*f[7    ][j] +  1470.0*f[6    ][j] -   2646.0*f[5    ][j] +   4410.0*f[4    ][j] -   1914.0*f[3    ][j] -   945.0*f[2    ][j] +   105.0*f[1    ][j] -    7.0*f[0    ][j] ) / ( 2520.0*h_x );
		f_dx[N_x-4][j] = (   -3.0*f[N_x-11][j] +   35.0*f[N_x-10][j] -   189.0*f[N_x-9][j] +   630.0*f[N_x-8][j] -  1470.0*f[N_x-7][j] +   2646.0*f[N_x-6][j] -   4410.0*f[N_x-5][j] +   1914.0*f[N_x-4][j] +   945.0*f[N_x-3][j] -   105.0*f[N_x-2][j] +    7.0*f[N_x-1][j] ) / ( 2520.0*h_x );
		
		f_dx[4    ][j] = (   -2.0*f[10    ][j] +   24.0*f[9     ][j] -   135.0*f[8    ][j] +   480.0*f[7    ][j] -  1260.0*f[6    ][j] +   3024.0*f[5    ][j] -    924.0*f[4    ][j] -   1440.0*f[3    ][j] +   270.0*f[2    ][j] -    40.0*f[1    ][j] +    3.0*f[0    ][j] ) / ( 2520.0*h_x );
		f_dx[N_x-5][j] = (    2.0*f[N_x-11][j] -   24.0*f[N_x-10][j] +   135.0*f[N_x-9][j] -   480.0*f[N_x-8][j] +  1260.0*f[N_x-7][j] -   3024.0*f[N_x-6][j] +    924.0*f[N_x-5][j] +   1440.0*f[N_x-4][j] -   270.0*f[N_x-3][j] +    40.0*f[N_x-2][j] -    3.0*f[N_x-1][j] ) / ( 2520.0*h_x );
		
		for( int n=5; n<=N_x-6; n++ ){
			f_dx[n][j] = ( 2.0*(f[n+5][j]-f[n-5][j]) - 25.0*(f[n+4][j]-f[n-4][j]) + 150.0*(f[n+3][j]-f[n-3][j]) - 600.0*(f[n+2][j]-f[n-2][j]) + 2100.0*(f[n+1][j]-f[n-1][j]) ) / ( 2520.0*h_x );
		}
	
	}
	
	return f_dx;
	
}

std::vector< std::vector<double> > f_dy_FD_function_order_10( std::vector< std::vector<double> > &f, double h_y ){
	
	int N_x = f   .size();
	int N_y = f[0].size();
	std::vector< std::vector<double> > f_dy ( N_x, std::vector<double> ( N_y ) );
	
	for( int i=0; i<N_x; i++ ){
	
		f_dy[i][0    ] = ( -252.0*f[i][10    ] + 2800.0*f[i][9     ] - 14175.0*f[i][8    ] + 43200.0*f[i][7    ] - 88200.0*f[i][6    ] + 127008.0*f[i][5    ] - 132300.0*f[i][4    ] + 100800.0*f[i][3    ] - 56700.0*f[i][2    ] + 25200.0*f[i][1    ] - 7381.0*f[i][0    ] ) / ( 2520.0*h_y );
		f_dy[i][N_y-1] = (  252.0*f[i][N_y-11] - 2800.0*f[i][N_y-10] + 14175.0*f[i][N_y-9] - 43200.0*f[i][N_y-8] + 88200.0*f[i][N_y-7] - 127008.0*f[i][N_y-6] + 132300.0*f[i][N_y-5] - 100800.0*f[i][N_y-4] + 56700.0*f[i][N_y-3] - 25200.0*f[i][N_y-2] + 7381.0*f[i][N_y-1] ) / ( 2520.0*h_y );
		
		f_dy[i][1    ] = (   28.0*f[i][10    ] -  315.0*f[i][9     ] +  1620.0*f[i][8    ] -  5040.0*f[i][7    ] + 10584.0*f[i][6    ] -  15876.0*f[i][5    ] +  17640.0*f[i][4    ] -  15120.0*f[i][3    ] + 11340.0*f[i][2    ] -  4609.0*f[i][1    ] -  252.0*f[i][0    ] ) / ( 2520.0*h_y );
		f_dy[i][N_y-2] = (  -28.0*f[i][N_y-11] +  315.0*f[i][N_y-10] -  1620.0*f[i][N_y-9] +  5040.0*f[i][N_y-8] - 10584.0*f[i][N_y-7] +  15876.0*f[i][N_y-6] -  17640.0*f[i][N_y-5] +  15120.0*f[i][N_y-4] - 11340.0*f[i][N_y-3] +  4609.0*f[i][N_y-2] +  252.0*f[i][N_y-1] ) / ( 2520.0*h_y );
		
		f_dy[i][2    ] = (   -7.0*f[i][10    ] +   80.0*f[i][9     ] -   420.0*f[i][8    ] +  1344.0*f[i][7    ] -  2940.0*f[i][6    ] +   4704.0*f[i][5    ] -   5880.0*f[i][4    ] +   6720.0*f[i][3    ] -  3069.0*f[i][2    ] -   560.0*f[i][1    ] +   28.0*f[i][0    ] ) / ( 2520.0*h_y );
		f_dy[i][N_y-3] = (    7.0*f[i][N_y-11] -   80.0*f[i][N_y-10] +   420.0*f[i][N_y-9] -  1344.0*f[i][N_y-8] +  2940.0*f[i][N_y-7] -   4704.0*f[i][N_y-6] +   5880.0*f[i][N_y-5] -   6720.0*f[i][N_y-4] +  3069.0*f[i][N_y-3] +   560.0*f[i][N_y-2] -   28.0*f[i][N_y-1] ) / ( 2520.0*h_y );
		
		f_dy[i][3    ] = (    3.0*f[i][10    ] -   35.0*f[i][9     ] +   189.0*f[i][8    ] -   630.0*f[i][7    ] +  1470.0*f[i][6    ] -   2646.0*f[i][5    ] +   4410.0*f[i][4    ] -   1914.0*f[i][3    ] -   945.0*f[i][2    ] +   105.0*f[i][1    ] -    7.0*f[i][0    ] ) / ( 2520.0*h_y );
		f_dy[i][N_y-4] = (   -3.0*f[i][N_y-11] +   35.0*f[i][N_y-10] -   189.0*f[i][N_y-9] +   630.0*f[i][N_y-8] -  1470.0*f[i][N_y-7] +   2646.0*f[i][N_y-6] -   4410.0*f[i][N_y-5] +   1914.0*f[i][N_y-4] +   945.0*f[i][N_y-3] -   105.0*f[i][N_y-2] +    7.0*f[i][N_y-1] ) / ( 2520.0*h_y );
		
		f_dy[i][4    ] = (   -2.0*f[i][10    ] +   24.0*f[i][9     ] -   135.0*f[i][8    ] +   480.0*f[i][7    ] -  1260.0*f[i][6    ] +   3024.0*f[i][5    ] -    924.0*f[i][4    ] -   1440.0*f[i][3    ] +   270.0*f[i][2    ] -    40.0*f[i][1    ] +    3.0*f[i][0    ] ) / ( 2520.0*h_y );
		f_dy[i][N_y-5] = (    2.0*f[i][N_y-11] -   24.0*f[i][N_y-10] +   135.0*f[i][N_y-9] -   480.0*f[i][N_y-8] +  1260.0*f[i][N_y-7] -   3024.0*f[i][N_y-6] +    924.0*f[i][N_y-5] +   1440.0*f[i][N_y-4] -   270.0*f[i][N_y-3] +    40.0*f[i][N_y-2] -    3.0*f[i][N_y-1] ) / ( 2520.0*h_y );
		
		for( int n=5; n<=N_y-6; n++ ){
			f_dy[i][n] = ( 2.0*(f[i][n+5]-f[i][n-5]) - 25.0*(f[i][n+4]-f[i][n-4]) + 150.0*(f[i][n+3]-f[i][n-3]) - 600.0*(f[i][n+2]-f[i][n-2]) + 2100.0*(f[i][n+1]-f[i][n-1]) ) / ( 2520.0*h_y );
		}
	
	}
	
	return f_dy;
	
}

std::vector< std::vector<double> > f_dx_FD_function_order_10_decimal( std::vector< std::vector<double> > &f, double h ){
	// Replace the division of large integers by each other, by decimals, in case this causes numerical errors.
	
	int N_x = f   .size();
	int N_y = f[0].size();
	std::vector< std::vector<double> > f_dx ( N_x, std::vector<double> ( N_y ) );
	double one_over_h = 1.0 / h;
	
	for( int j=0; j<N_y; j++ ){
	
		f_dx[0    ][j] = -0.1*f[10    ][j] + (10.0/9.0)*f[9     ][j] - 5.625*f[8    ][j] + (120.0/7.0)*f[7    ][j] - 35.0*f[6    ][j] + 50.4*f[5    ][j] - 52.5*f[4    ][j] + 40.0*f[3    ][j] - 22.5*f[2    ][j] + 10.0*f[1    ][j] - (7381.0/2520.0)*f[0    ][j];
		f_dx[N_x-1][j] =  0.1*f[N_x-11][j] - (10.0/9.0)*f[N_x-10][j] + 5.625*f[N_x-9][j] - (120.0/7.0)*f[N_x-8][j] + 35.0*f[N_x-7][j] - 50.4*f[N_x-6][j] + 52.5*f[N_x-5][j] - 40.0*f[N_x-4][j] + 22.5*f[N_x-3][j] - 10.0*f[N_x-2][j] + (7381.0/2520.0)*f[N_x-1][j];
		
		f_dx[1    ][j] =  (1.0/90.0)*f[10    ][j] -  0.125*f[9     ][j] +  (9.0/14.0)*f[8    ][j] - 2.0*f[7    ][j] + 4.2*f[6    ][j] -  6.3*f[5    ][j] +  7.0*f[4    ][j] -  6.0*f[3    ][j] + 4.5*f[2    ][j] -  (4609.0/2520.0)*f[1    ][j] - 0.1*f[0    ][j];
		f_dx[N_x-2][j] = -(1.0/90.0)*f[N_x-11][j] +  0.125*f[N_x-10][j] -  (9.0/14.0)*f[N_x-9][j] + 2.0*f[N_x-8][j] - 4.2*f[N_x-7][j] +  6.3*f[N_x-6][j] -  7.0*f[N_x-5][j] +  6.0*f[N_x-4][j] - 4.5*f[N_x-3][j] +  (4609.0/2520.0)*f[N_x-2][j] + 0.1*f[N_x-1][j];
		
		f_dx[2    ][j] = -(1.0/360.0)*f[10    ][j] + (2.0/63.0)*f[9     ][j] - (1.0/6.0)*f[8    ][j] + (8.0/15.0)*f[7    ][j] -  (7.0/6.0)*f[6    ][j] +  (28.0/15.0)*f[5    ][j] -  (7.0/3.0)*f[4    ][j] + (8.0/3.0)*f[3    ][j] - (341.0/280.0)*f[2    ][j] - (2.0/9.0)*f[1    ][j] + (1.0/90.0)*f[0    ][j];
		f_dx[N_x-3][j] =  (1.0/360.0)*f[N_x-11][j] - (2.0/63.0)*f[N_x-10][j] + (1.0/6.0)*f[N_x-9][j] - (8.0/15.0)*f[N_x-8][j] +  (7.0/6.0)*f[N_x-7][j] -  (28.0/15.0)*f[N_x-6][j] +  (7.0/3.0)*f[N_x-5][j] - (8.0/3.0)*f[N_x-4][j] + (341.0/280.0)*f[N_x-3][j] + (2.0/9.0)*f[N_x-2][j] - (1.0/90.0)*f[N_x-1][j];
		
		f_dx[3    ][j] =  (1.0/840.0)*f[10    ][j] - (1.0/72.0)*f[9     ][j] + 0.075*f[8    ][j] - 0.25*f[7    ][j] + (7.0/12.0)*f[6    ][j] - 1.05*f[5    ][j] + 1.75*f[4    ][j] - (319.0/420.0)*f[3    ][j] - 0.375*f[2    ][j] + (1.0/24.0)*f[1    ][j] - (1.0/360.0)*f[0    ][j];
		f_dx[N_x-4][j] = -(1.0/840.0)*f[N_x-11][j] + (1.0/72.0)*f[N_x-10][j] - 0.075*f[N_x-9][j] + 0.25*f[N_x-8][j] - (7.0/12.0)*f[N_x-7][j] + 1.05*f[N_x-6][j] - 1.75*f[N_x-5][j] + (319.0/420.0)*f[N_x-4][j] + 0.375*f[N_x-3][j] - (1.0/24.0)*f[N_x-2][j] + (1.0/360.0)*f[N_x-1][j];
		
		f_dx[4    ][j] = -(1.0/1260.0)*f[10    ][j] + (1.0/105.0)*f[9     ][j] - (3.0/56.0)*f[8    ][j] + (4.0/21.0)*f[7    ][j] - 0.5*f[6    ][j] + 1.2*f[5    ][j] - (11.0/30.0)*f[4    ][j] - (4.0/7.0)*f[3    ][j] + (3.0/28.0)*f[2    ][j] - (1.0/63.0)*f[1    ][j] + (1.0/840.0)*f[0    ][j];
		f_dx[N_x-5][j] =  (1.0/1260.0)*f[N_x-11][j] - (1.0/105.0)*f[N_x-10][j] + (3.0/56.0)*f[N_x-9][j] - (4.0/21.0)*f[N_x-8][j] + 0.5*f[N_x-7][j] - 1.2*f[N_x-6][j] + (11.0/30.0)*f[N_x-5][j] + (4.0/7.0)*f[N_x-4][j] - (3.0/28.0)*f[N_x-3][j] + (1.0/63.0)*f[N_x-2][j] - (1.0/840.0)*f[N_x-1][j];
		
		for( int i=5; i<=N_x-6; i++ ){
			f_dx[i][j] = (1.0/1260.0)*(f[i+5][j]-f[i-5][j]) - (5.0/504.0)*(f[i+4][j]-f[i-4][j]) + (5.0/84.0)*(f[i+3][j]-f[i-3][j]) - (5.0/21.0)*(f[i+2][j]-f[i-2][j]) + (5.0/6.0)*(f[i+1][j]-f[i-1][j]);
		}
		
		for( int i=0; i<N_x; i++ ){
			f_dx[i][j] *= one_over_h;
		}
	
	}
	
	return f_dx;
	
}

std::vector< std::vector<double> > f_dy_FD_function_order_10_decimal( std::vector< std::vector<double> > &f, double h ){
	// Replace the division of large integers by each other, by decimals, in case this causes numerical errors.
	
	int N_x = f   .size();
	int N_y = f[0].size();
	std::vector< std::vector<double> > f_dy ( N_x, std::vector<double> ( N_y ) );
	double one_over_h = 1.0 / h;
	
	for( int i=0; i<N_x; i++ ){
	
		f_dy[i][0    ] = -0.1*f[i][10    ] + (10.0/9.0)*f[i][9     ] - 5.625*f[i][8    ] + (120.0/7.0)*f[i][7    ] - 35.0*f[i][6    ] + 50.4*f[i][5    ] - 52.5*f[i][4    ] + 40.0*f[i][3    ] - 22.5*f[i][2    ] + 10.0*f[i][1    ] - (7381.0/2520.0)*f[i][0    ];
		f_dy[i][N_y-1] =  0.1*f[i][N_y-11] - (10.0/9.0)*f[i][N_y-10] + 5.625*f[i][N_y-9] - (120.0/7.0)*f[i][N_y-8] + 35.0*f[i][N_y-7] - 50.4*f[i][N_y-6] + 52.5*f[i][N_y-5] - 40.0*f[i][N_y-4] + 22.5*f[i][N_y-3] - 10.0*f[i][N_y-2] + (7381.0/2520.0)*f[i][N_y-1];
		
		f_dy[i][1    ] =  (1.0/90.0)*f[i][10    ] -  0.125*f[i][9     ] +  (9.0/14.0)*f[i][8    ] - 2.0*f[i][7    ] + 4.2*f[i][6    ] -  6.3*f[i][5    ] +  7.0*f[i][4    ] -  6.0*f[i][3    ] + 4.5*f[i][2    ] -  (4609.0/2520.0)*f[i][1    ] - 0.1*f[i][0    ];
		f_dy[i][N_y-2] = -(1.0/90.0)*f[i][N_y-11] +  0.125*f[i][N_y-10] -  (9.0/14.0)*f[i][N_y-9] + 2.0*f[i][N_y-8] - 4.2*f[i][N_y-7] +  6.3*f[i][N_y-6] -  7.0*f[i][N_y-5] +  6.0*f[i][N_y-4] - 4.5*f[i][N_y-3] +  (4609.0/2520.0)*f[i][N_y-2] + 0.1*f[i][N_y-1];
		
		f_dy[i][2    ] = -(1.0/360.0)*f[i][10    ] + (2.0/63.0)*f[i][9     ] - (1.0/6.0)*f[i][8    ] + (8.0/15.0)*f[i][7    ] -  (7.0/6.0)*f[i][6    ] +  (28.0/15.0)*f[i][5    ] -  (7.0/3.0)*f[i][4    ] + (8.0/3.0)*f[i][3    ] - (341.0/280.0)*f[i][2    ] - (2.0/9.0)*f[i][1    ] + (1.0/90.0)*f[i][0    ];
		f_dy[i][N_y-3] =  (1.0/360.0)*f[i][N_y-11] - (2.0/63.0)*f[i][N_y-10] + (1.0/6.0)*f[i][N_y-9] - (8.0/15.0)*f[i][N_y-8] +  (7.0/6.0)*f[i][N_y-7] -  (28.0/15.0)*f[i][N_y-6] +  (7.0/3.0)*f[i][N_y-5] - (8.0/3.0)*f[i][N_y-4] + (341.0/280.0)*f[i][N_y-3] + (2.0/9.0)*f[i][N_y-2] - (1.0/90.0)*f[i][N_y-1];
		
		f_dy[i][3    ] =  (1.0/840.0)*f[i][10    ] - (1.0/72.0)*f[i][9     ] + 0.075*f[i][8    ] - 0.25*f[i][7    ] + (7.0/12.0)*f[i][6    ] - 1.05*f[i][5    ] + 1.75*f[i][4    ] - (319.0/420.0)*f[i][3    ] - 0.375*f[i][2    ] + (1.0/24.0)*f[i][1    ] - (1.0/360.0)*f[i][0    ];
		f_dy[i][N_y-4] = -(1.0/840.0)*f[i][N_y-11] + (1.0/72.0)*f[i][N_y-10] - 0.075*f[i][N_y-9] + 0.25*f[i][N_y-8] - (7.0/12.0)*f[i][N_y-7] + 1.05*f[i][N_y-6] - 1.75*f[i][N_y-5] + (319.0/420.0)*f[i][N_y-4] + 0.375*f[i][N_y-3] - (1.0/24.0)*f[i][N_y-2] + (1.0/360.0)*f[i][N_y-1];
		
		f_dy[i][4    ] = -(1.0/1260.0)*f[i][10    ] + (1.0/105.0)*f[i][9     ] - (3.0/56.0)*f[i][8    ] + (4.0/21.0)*f[i][7    ] - 0.5*f[i][6    ] + 1.2*f[i][5    ] - (11.0/30.0)*f[i][4    ] - (4.0/7.0)*f[i][3    ] + (3.0/28.0)*f[i][2    ] - (1.0/63.0)*f[i][1    ] + (1.0/840.0)*f[i][0    ];
		f_dy[i][N_y-5] =  (1.0/1260.0)*f[i][N_y-11] - (1.0/105.0)*f[i][N_y-10] + (3.0/56.0)*f[i][N_y-9] - (4.0/21.0)*f[i][N_y-8] + 0.5*f[i][N_y-7] - 1.2*f[i][N_y-6] + (11.0/30.0)*f[i][N_y-5] + (4.0/7.0)*f[i][N_y-4] - (3.0/28.0)*f[i][N_y-3] + (1.0/63.0)*f[i][N_y-2] - (1.0/840.0)*f[i][N_y-1];
		
		for( int j=5; j<=N_y-6; j++ ){
			f_dy[i][j] = (1.0/1260.0)*(f[i][j+5]-f[i][j-5]) - (5.0/504.0)*(f[i][j+4]-f[i][j-4]) + (5.0/84.0)*(f[i][j+3]-f[i][j-3]) - (5.0/21.0)*(f[i][j+2]-f[i][j-2]) + (5.0/6.0)*(f[i][j+1]-f[i][j-1]);
		}
		
		for( int j=0; j<N_y; j++ ){
			f_dy[i][j] *= one_over_h;
		}
	
	}
	
	return f_dy;
	
}