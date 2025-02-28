// cd OneDrive\PhD\Codes\Test mathematical expressions
// g++ sum_of_chebyshevs_and_spherical_harmonics.cpp -o sum_of_chebyshevs_and_spherical_harmonics

// Automated function to give a sum of Chebyshev polynomials and spherical harmonics with slowly incrementing multiplying factors.
// This guarantees a final function that depends on all of the individual basis functions, so that the coefficient-fitting code can be tested without fear of missing an error.

#include <iostream>
#include <math.h>
#include "../Generate_Basis_Functions.h"
#include "../Hardcoded_Chebyshev_Polynomials.h"
#include "../Hardcoded_Spherical_Harmonics.h"


//----- Global variables and function declarations -----
double LHS( int n );
double RHS( int n );

std::complex<double> f( double r, double t, double p, int n_max, int ell_max );
std::complex<double> f_hardcoded( double r, double t, double p );


int main(){
	
	//----- Define variables -----
	int n_max   = 3;	// Maximum order of Chebyshev polynomials.
	int ell_max = 3;	// Maximum order of spherical harmonics.
	
	double r = 0.3;
	double t = 0.32;
	double p = 0.53;
	
	
	//----- Generate basis functions -----
	generate_coeffs_chebyshev( n_max );
	generate_coeffs_spherical_harmonic( ell_max );
	
	
	
	
	//----- Output results -----
	std::cout << f( r, t, p, n_max, ell_max ) <<"\t"<< f_hardcoded( r, t, p ) << std::endl;
	
	
	return 0;
}




std::complex<double> f( double r, double t, double p, int n_max, int ell_max ){
	std::complex<double> ret = 0;
	double a  = 1.0;	// Real part of coefficient.
	double da = 0.1;	// Increment to real part of coefficient.
	double b  = 9.9;	// Imaginary part of coefficient.
	double db = -0.1;	// Increment to imaginary part of coefficient.
	double a_tot = 0;
	double b_tot = 0;
	for( int n=0; n<=n_max; n++ ){
		for( int ell=0; ell<=ell_max; ell++ ){
			for( int m=-ell; m<=ell; m++ ){
				std::cout <<"("<< a <<"+"<< b <<"i) * T_"<< n <<" * Y_"<< ell <<"^"<< m << std::endl;
				ret += std::complex<double>{a,b} * polynomial_chebyshev_1( r, n ) * ylm( t, p, ell, m );
				a_tot += a;
				b_tot += b;
				a += da;
				b += db;
			}
		}
	}
	std::cout << "Coeff total:\t" << std::complex<double>{a_tot,b_tot} << std::endl;
	return ret;
}




std::complex<double> f_hardcoded( double r, double t, double p ){
	
	//----- Case (0,0) -----
	//return std::complex<double>{1.0,9.9}*t0(r)*y00 (t,p);
	
	//----- Case (1,1)-----
	/*return std::complex<double>{1.0,9.9}*t0(r)*y00 (t,p)
	     + std::complex<double>{1.1,9.8}*t0(r)*y1m1(t,p) + std::complex<double>{1.2,9.7}*t0(r)*y10 (t,p) + std::complex<double>{1.3,9.6}*t0(r)*y11 (t,p)
		 + std::complex<double>{1.4,9.5}*t1(r)*y00 (t,p)
		 + std::complex<double>{1.5,9.4}*t1(r)*y1m1(t,p) + std::complex<double>{1.6,9.3}*t1(r)*y10 (t,p) + std::complex<double>{1.7,9.2}*t1(r)*y11 (t,p);*/
	
	//----- Case (1,3) -----
	/*return std::complex<double>{1.0,9.9}*t0(r)*y00 (t,p)
	     + std::complex<double>{1.1,9.8}*t0(r)*y1m1(t,p) + std::complex<double>{1.2,9.7}*t0(r)*y10 (t,p) + std::complex<double>{1.3,9.6}*t0(r)*y11 (t,p)
		 + std::complex<double>{1.4,9.5}*t0(r)*y2m2(t,p) + std::complex<double>{1.5,9.4}*t0(r)*y2m1(t,p) + std::complex<double>{1.6,9.3}*t0(r)*y20 (t,p) + std::complex<double>{1.7,9.2}*t0(r)*y21(t,p) + std::complex<double>{1.8,9.1}*t0(r)*y22(t,p)
		 + std::complex<double>{1.9,9.0}*t0(r)*y3m3(t,p) + std::complex<double>{2.0,8.9}*t0(r)*y3m2(t,p) + std::complex<double>{2.1,8.8}*t0(r)*y3m1(t,p) + std::complex<double>{2.2,8.7}*t0(r)*y30(t,p) + std::complex<double>{2.3,8.6}*t0(r)*y31(t,p) + std::complex<double>{2.4,8.5}*t0(r)*y32(t,p) + std::complex<double>{2.5,8.4}*t0(r)*y33(t,p)
		 + std::complex<double>{2.6,8.3}*t1(r)*y00 (t,p)
		 + std::complex<double>{2.7,8.2}*t1(r)*y1m1(t,p) + std::complex<double>{2.8,8.1}*t1(r)*y10 (t,p) + std::complex<double>{2.9,8.0}*t1(r)*y11 (t,p)
		 + std::complex<double>{3.0,7.9}*t1(r)*y2m2(t,p) + std::complex<double>{3.1,7.8}*t1(r)*y2m1(t,p) + std::complex<double>{3.2,7.7}*t1(r)*y20 (t,p) + std::complex<double>{3.3,7.6}*t1(r)*y21(t,p) + std::complex<double>{3.4,7.5}*t1(r)*y22(t,p)
		 + std::complex<double>{3.5,7.4}*t1(r)*y3m3(t,p) + std::complex<double>{3.6,7.3}*t1(r)*y3m2(t,p) + std::complex<double>{3.7,7.2}*t1(r)*y3m1(t,p) + std::complex<double>{3.8,7.1}*t1(r)*y30(t,p) + std::complex<double>{3.9,7.0}*t1(r)*y31(t,p) + std::complex<double>{4.0,6.9}*t1(r)*y32(t,p) + std::complex<double>{4.1,6.8}*t1(r)*y33(t,p);*/
		 
	//----- Case (2,3) -----
	/*return std::complex<double>{1.0,9.9}*t0(r)*y00 (t,p)
	     + std::complex<double>{1.1,9.8}*t0(r)*y1m1(t,p) + std::complex<double>{1.2,9.7}*t0(r)*y10 (t,p) + std::complex<double>{1.3,9.6}*t0(r)*y11 (t,p)
		 + std::complex<double>{1.4,9.5}*t0(r)*y2m2(t,p) + std::complex<double>{1.5,9.4}*t0(r)*y2m1(t,p) + std::complex<double>{1.6,9.3}*t0(r)*y20 (t,p) + std::complex<double>{1.7,9.2}*t0(r)*y21(t,p) + std::complex<double>{1.8,9.1}*t0(r)*y22(t,p)
		 + std::complex<double>{1.9,9.0}*t0(r)*y3m3(t,p) + std::complex<double>{2.0,8.9}*t0(r)*y3m2(t,p) + std::complex<double>{2.1,8.8}*t0(r)*y3m1(t,p) + std::complex<double>{2.2,8.7}*t0(r)*y30(t,p) + std::complex<double>{2.3,8.6}*t0(r)*y31(t,p) + std::complex<double>{2.4,8.5}*t0(r)*y32(t,p) + std::complex<double>{2.5,8.4}*t0(r)*y33(t,p)
		 + std::complex<double>{2.6,8.3}*t1(r)*y00 (t,p)
		 + std::complex<double>{2.7,8.2}*t1(r)*y1m1(t,p) + std::complex<double>{2.8,8.1}*t1(r)*y10 (t,p) + std::complex<double>{2.9,8.0}*t1(r)*y11 (t,p)
		 + std::complex<double>{3.0,7.9}*t1(r)*y2m2(t,p) + std::complex<double>{3.1,7.8}*t1(r)*y2m1(t,p) + std::complex<double>{3.2,7.7}*t1(r)*y20 (t,p) + std::complex<double>{3.3,7.6}*t1(r)*y21(t,p) + std::complex<double>{3.4,7.5}*t1(r)*y22(t,p)
		 + std::complex<double>{3.5,7.4}*t1(r)*y3m3(t,p) + std::complex<double>{3.6,7.3}*t1(r)*y3m2(t,p) + std::complex<double>{3.7,7.2}*t1(r)*y3m1(t,p) + std::complex<double>{3.8,7.1}*t1(r)*y30(t,p) + std::complex<double>{3.9,7.0}*t1(r)*y31(t,p) + std::complex<double>{4.0,6.9}*t1(r)*y32(t,p) + std::complex<double>{4.1,6.8}*t1(r)*y33(t,p)
		 + std::complex<double>{4.2,6.7}*t2(r)*y00 (t,p)
		 + std::complex<double>{4.3,6.6}*t2(r)*y1m1(t,p) + std::complex<double>{4.4,6.5}*t2(r)*y10 (t,p) + std::complex<double>{4.5,6.4}*t2(r)*y11 (t,p)
		 + std::complex<double>{4.6,6.3}*t2(r)*y2m2(t,p) + std::complex<double>{4.7,6.2}*t2(r)*y2m1(t,p) + std::complex<double>{4.8,6.1}*t2(r)*y20 (t,p) + std::complex<double>{4.9,6.0}*t2(r)*y21(t,p) + std::complex<double>{5.0,5.9}*t2(r)*y22(t,p)
		 + std::complex<double>{5.1,5.8}*t2(r)*y3m3(t,p) + std::complex<double>{5.2,5.7}*t2(r)*y3m2(t,p) + std::complex<double>{5.3,5.6}*t2(r)*y3m1(t,p) + std::complex<double>{5.4,5.5}*t2(r)*y30(t,p) + std::complex<double>{5.5,5.4}*t2(r)*y31(t,p) + std::complex<double>{5.6,5.3}*t2(r)*y32(t,p) + std::complex<double>{5.7,5.2}*t2(r)*y33(t,p);*/
	
	//----- Case (3,3) -----
	//----- Case (2,3) -----
	return std::complex<double>{1.0,9.9}*t0(r)*y00 (t,p)
	     + std::complex<double>{1.1,9.8}*t0(r)*y1m1(t,p) + std::complex<double>{1.2,9.7}*t0(r)*y10 (t,p) + std::complex<double>{1.3,9.6}*t0(r)*y11 (t,p)
		 + std::complex<double>{1.4,9.5}*t0(r)*y2m2(t,p) + std::complex<double>{1.5,9.4}*t0(r)*y2m1(t,p) + std::complex<double>{1.6,9.3}*t0(r)*y20 (t,p) + std::complex<double>{1.7,9.2}*t0(r)*y21(t,p) + std::complex<double>{1.8,9.1}*t0(r)*y22(t,p)
		 + std::complex<double>{1.9,9.0}*t0(r)*y3m3(t,p) + std::complex<double>{2.0,8.9}*t0(r)*y3m2(t,p) + std::complex<double>{2.1,8.8}*t0(r)*y3m1(t,p) + std::complex<double>{2.2,8.7}*t0(r)*y30(t,p) + std::complex<double>{2.3,8.6}*t0(r)*y31(t,p) + std::complex<double>{2.4,8.5}*t0(r)*y32(t,p) + std::complex<double>{2.5,8.4}*t0(r)*y33(t,p)
		 + std::complex<double>{2.6,8.3}*t1(r)*y00 (t,p)
		 + std::complex<double>{2.7,8.2}*t1(r)*y1m1(t,p) + std::complex<double>{2.8,8.1}*t1(r)*y10 (t,p) + std::complex<double>{2.9,8.0}*t1(r)*y11 (t,p)
		 + std::complex<double>{3.0,7.9}*t1(r)*y2m2(t,p) + std::complex<double>{3.1,7.8}*t1(r)*y2m1(t,p) + std::complex<double>{3.2,7.7}*t1(r)*y20 (t,p) + std::complex<double>{3.3,7.6}*t1(r)*y21(t,p) + std::complex<double>{3.4,7.5}*t1(r)*y22(t,p)
		 + std::complex<double>{3.5,7.4}*t1(r)*y3m3(t,p) + std::complex<double>{3.6,7.3}*t1(r)*y3m2(t,p) + std::complex<double>{3.7,7.2}*t1(r)*y3m1(t,p) + std::complex<double>{3.8,7.1}*t1(r)*y30(t,p) + std::complex<double>{3.9,7.0}*t1(r)*y31(t,p) + std::complex<double>{4.0,6.9}*t1(r)*y32(t,p) + std::complex<double>{4.1,6.8}*t1(r)*y33(t,p)
		 + std::complex<double>{4.2,6.7}*t2(r)*y00 (t,p)
		 + std::complex<double>{4.3,6.6}*t2(r)*y1m1(t,p) + std::complex<double>{4.4,6.5}*t2(r)*y10 (t,p) + std::complex<double>{4.5,6.4}*t2(r)*y11 (t,p)
		 + std::complex<double>{4.6,6.3}*t2(r)*y2m2(t,p) + std::complex<double>{4.7,6.2}*t2(r)*y2m1(t,p) + std::complex<double>{4.8,6.1}*t2(r)*y20 (t,p) + std::complex<double>{4.9,6.0}*t2(r)*y21(t,p) + std::complex<double>{5.0,5.9}*t2(r)*y22(t,p)
		 + std::complex<double>{5.1,5.8}*t2(r)*y3m3(t,p) + std::complex<double>{5.2,5.7}*t2(r)*y3m2(t,p) + std::complex<double>{5.3,5.6}*t2(r)*y3m1(t,p) + std::complex<double>{5.4,5.5}*t2(r)*y30(t,p) + std::complex<double>{5.5,5.4}*t2(r)*y31(t,p) + std::complex<double>{5.6,5.3}*t2(r)*y32(t,p) + std::complex<double>{5.7,5.2}*t2(r)*y33(t,p)
		 + std::complex<double>{5.8,5.1}*t3(r)*y00 (t,p)
		 + std::complex<double>{5.9,5.0}*t3(r)*y1m1(t,p) + std::complex<double>{6.0,4.9}*t3(r)*y10 (t,p) + std::complex<double>{6.1,4.8}*t3(r)*y11 (t,p)
		 + std::complex<double>{6.2,4.7}*t3(r)*y2m2(t,p) + std::complex<double>{6.3,4.6}*t3(r)*y2m1(t,p) + std::complex<double>{6.4,4.5}*t3(r)*y20 (t,p) + std::complex<double>{6.5,4.4}*t3(r)*y21(t,p) + std::complex<double>{6.6,4.3}*t3(r)*y22(t,p)
		 + std::complex<double>{6.7,4.2}*t3(r)*y3m3(t,p) + std::complex<double>{6.8,4.1}*t3(r)*y3m2(t,p) + std::complex<double>{6.9,4.0}*t3(r)*y3m1(t,p) + std::complex<double>{7.0,3.9}*t3(r)*y30(t,p) + std::complex<double>{7.1,3.8}*t3(r)*y31(t,p) + std::complex<double>{7.2,3.7}*t3(r)*y32(t,p) + std::complex<double>{7.3,3.6}*t3(r)*y33(t,p);
}