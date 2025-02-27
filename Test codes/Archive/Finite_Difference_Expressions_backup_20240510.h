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
