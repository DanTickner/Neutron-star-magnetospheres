#include <math.h>

//----- Chebyshev polynomials of the first kind -----
double t0( double x ) {
	return 1;
}
double t1( double x ) {
	return x;
}
double t2( double x ) {
	return 2*x*x - 1;
}
double t3( double x ) {
	return 4*x*x*x - 3*x;
}
double t4( double x ) {
	return 8*pow(x,4) - 8*x*x + 1;
}
double t5( double x ) {
	return 16*pow(x,5) - 20*x*x*x + 5*x;
}
double t6( double x ) {
	return 32*pow(x,6) - 48*pow(x,4) + 18*x*x - 1;
}
double t7( double x ) {
	return 64*pow(x,7) - 112*pow(x,5) + 56*x*x*x - 7*x;
}
double t8( double x ) {
	return 128*pow(x,8) - 256*pow(x,6) + 160*pow(x,4) - 32*x*x + 1;
}
double t9( double x ) {
	return 256*pow(x,9) - 576*pow(x,7) + 432*pow(x,5) - 120*x*x*x + 9*x;
}
double t10( double x ) {
	return 512*pow(x,10) - 1280*pow(x,8) + 1120*pow(x,6) - 400*pow(x,4) + 50*x*x - 1;
}




//----- Chebyshev polynomials of the second kind -----
double u0( double x ) {
	return 1;
}
double u1( double x ) {
	return 2*x;
}
double u2( double x ) {
	return 4*x*x - 1;
}
double u3( double x ) {
	return 8*x*x*x - 4*x;
}
double u4( double x ) {
	return 16*pow(x,4) - 12*x*x + 1;
}
double u5( double x ) {
	return 32*pow(x,5) - 32*x*x*x + 6*x;
}
double u6( double x ) {
	return 64*pow(x,6) - 80*pow(x,4) + 24*x*x - 1;
}
double u7( double x ) {
	return 128*pow(x,7) - 192*pow(x,5) + 80*x*x*x - 8*x;
}
double u8( double x ) {
	return 256*pow(x,8) - 448*pow(x,6) + 240*pow(x,4) - 40*x*x + 1;
}
double u9( double x ) {
	return 512*pow(x,9) - 1024*pow(x,7) + 672*pow(x,5) - 160*pow(x,3) + 10*x;
}
double u10( double x ) {
	return 1024*pow(x,10) - 2304*pow(x,8) + 1792*pow(x,6) - 560*pow(x,4) + 60*x*x - 1;
}




//----- First derivatives of Chebyshev polynomials of the first kind -----
double t0_dx( double x ) {
	return 0;
}
double t1_dx( double x ) {
	return 1;
}
double t2_dx( double x ) {
	return 4*x;
}
double t3_dx( double x ) {
	return 12*x*x - 3;
}
double t4_dx( double x ) {
	return 32*x*x*x - 16*x;
}
double t5_dx( double x ) {
	return 80*x*x*x*x - 60*x*x + 5;
}
double t6_dx( double x ) {
	return 192*pow(x,5) - 384*pow(x,3) + 36*x;
}
double t7_dx( double x ) {
	return 448*pow(x,6) - 560*pow(x,4) + 168*x*x - 7;
}
double t8_dx( double x ) {
	return 1024*pow(x,7) - 1536*pow(x,5) + 640*pow(x,3) - 64*x;
}
double t9_dx( double x ) {
	return 2304*pow(x,8) - 4032*pow(x,6) + 2160*pow(x,4) - 360*x*x + 9;
}
double t10_dx( double x ) {
	return 5120*pow(x,9) - 10240*pow(x,7) + 6720*pow(x,5) - 1600*pow(x,3) + 100*x;
}




//----- Second derivatives of Chebyshev polynomials of the first kind -----
double t0_dx2( double x ) {
	return 0;
}
double t1_dx2( double x ) {
	return 0;
}
double t2_dx2( double x ) {
	return 4;
}
double t3_dx2( double x ) {
	return 24*x;
}
double t4_dx2( double x ) {
	return 96*x*x - 16;
}
double t5_dx2( double x ) {
	return 320*x*x*x - 120*x;
}
double t6_dx2( double x ) {
	return 960*pow(x,4) - 1152*pow(x,2) + 36;
}
double t7_dx2( double x ) {
	return 2688*pow(x,5) - 2240*pow(x,3) + 336*x;
}
double t8_dx2( double x ) {
	return 7168*pow(x,6) - 7680*pow(x,4) + 1920*pow(x,2) - 64;
}
double t9_dx2( double x ) {
	return 18432*pow(x,7) - 24192*pow(x,5) + 8640*pow(x,3) - 720*x;
}
double t10_dx2( double x ) {
	return 46080*pow(x,8) - 71680*pow(x,6) + 33600*pow(x,4) - 4800*pow(x,2) + 100;
}




//----- First derivatives of Chebyshev polynomials of the second kind -----
double u0_dx( double x ) {
	return 0;
}
double u1_dx( double x ) {
	return 2;
}
double u2_dx( double x ) {
	return 8*x;
}
double u3_dx( double x ) {
	return 24*x*x - 4;
}
double u4_dx( double x ) {
	return 64*x*x*x - 24*x;
}
double u5_dx( double x ) {
	return 160*pow(x,4) - 96*x*x + 6;
}
double u6_dx( double x ) {
	return 384*pow(x,5) - 320*pow(x,3) + 48*x;
}
double u7_dx( double x ) {
	return 896*pow(x,6) - 960*pow(x,4) + 240*x*x - 8;
}
double u8_dx( double x ) {
	return 2048*pow(x,7) - 2688*pow(x,5) + 960*pow(x,3) - 80*x;
}
double u9_dx( double x ) {
	return 4608*pow(x,8) - 7168*pow(x,6) + 3360*pow(x,4) - 480*pow(x,2) + 10;
}
double u10_dx( double x ) {
	return 10240*pow(x,9) - 18432*pow(x,7) + 10752*pow(x,5) - 2240*pow(x,3) + 120*x;
}




//----- Second derivatives of Chebyshev polynomials of the second kind -----
double u0_dx2( double x ) {
	return 0;
}
double u1_dx2( double x ) {
	return 0;
}
double u2_dx2( double x ) {
	return 8;
}
double u3_dx2( double x ) {
	return 48*x;
}
double u4_dx2( double x ) {
	return 192*x*x - 24;
}
double u5_dx2( double x ) {
	return 640*pow(x,3) - 192*x;
}
double u6_dx2( double x ) {
	return 1920*pow(x,4) - 960*pow(x,2) + 48;
}
double u7_dx2( double x ) {
	return 5376*pow(x,5) - 3840*pow(x,3) + 480*x;
}
double u8_dx2( double x ) {
	return 14336*pow(x,6) - 13440*pow(x,4) + 2880*pow(x,2) - 80;
}
double u9_dx2( double x ) {
	return 36864*pow(x,7) - 43008*pow(x,5) + 13440*pow(x,3) - 960*x;
}
double u10_dx2( double x ) {
	return 92160*pow(x,8) - 129024*pow(x,6) + 53760*pow(x,4) - 6720*pow(x,2) + 120;
}