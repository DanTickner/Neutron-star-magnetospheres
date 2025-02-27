/*
Trigonometric_Integrals.h
Evaluate the sine integral and the cosine integral,
Si(x) = int_0^x sin(t)/t dt
Ci(x) = int_0^x cos(t)/t dt
using a recursive power series up to n_max terms. Set a default value for n_max.
Because we calculated recursive relations for the terms, there is minimal rist of values blowing up for large n.
Hence, good accuracy "can be counted on" even for large n. The series converge for all x in C, but does so slowly, so many terms are required.
See Book 14, p21-22.
*/


double Si( double x, int n_max=100 ){
	
	double prev = x;
	double ret  = prev;
	
	for( int n=1; n<=n_max; n++ ){
		double ratio = - x*x * ( 2.0*(double)n-1.0 ) / ( 2.0 * n * pow( 2.0*(double)n+1.0, 2 ) );
		double next = prev * ratio;
		ret += next;
		prev = next;
	}
	
	return ret;
}


double Ci( double x, int n_max=100 ){
	
	double prev = -0.25 * x*x;
	double ret  = prev;
	
	for( int n=2; n<=n_max; n++ ){
		double ratio = - x*x * ( n - 1.0 ) / ( 2.0 * n*n * ( 2.0*(double)n-1.0 ) );
		double next = prev * ratio;
		ret += next;
		prev = next;
	}

	double euler_mascheroni = 0.5772156649015328606065120900824024310421;
	return euler_mascheroni + log(x) + ret;
}