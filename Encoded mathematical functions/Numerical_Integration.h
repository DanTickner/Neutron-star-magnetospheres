#include <vector>
#include <math.h>
#include <complex>

//-------------------- ONE-DIMENSIONAL INTEGRALS --------------------
double rk4( double (*f)(double), double x_min, double x_max, int n_steps ) {
	// Fourth-order Runge-Kutta method for function only of x.
	double h   = ( x_max - x_min ) / n_steps;
	double x   = x_min+h;
	double ret = 0.5 * ( f(x_min) + 4*f(x_min+0.5*h) + f(x_max) );
	
	for( int i=1; i<n_steps; i++ ){
		ret += f(x) + 2*f(x+0.5*h);
		x   += h;
	}
	
	return ret * h/3.0;
}
std::complex<double> rk4( std::complex<double> (*f)(double), double x_min, double x_max, int n_steps ) {
	// Fourth-order Runge-Kutta method for complex-valued function only of x, where x is real.
	double               h   = ( x_max - x_min ) / n_steps;
	double               x   = x_min+h;
	std::complex<double> ret = 0.5 * ( f(x_min) + 4.0*f(x_min+0.5*h) + f(x_max) );
	
	for( int i=1; i<n_steps; i++ ){
		ret += f(x) + 2.0*f(x+0.5*h);
		x   += h;
	}
	
	return ret * h/3.0;
}
double rk4( double (*f)(double,double), double x_min, double y_min, double x_max, int n_steps ) {
	// Fourth-order Runge-Kutta method for function of x and y.
	double h = ( x_max - x_min ) / n_steps;
	double x = x_min;
	double y = y_min;
	
	for( int i=0; i<n_steps; i++ ){
		double k1 = h * f( x        , y          );
		double k2 = h * f( x + h*0.5, y + k1*0.5 );
		double k3 = h * f( x + h*0.5, y + k2*0.5 );
		double k4 = h * f( x + h    , y + k3     );
		y += ( k1 + k2*2 + k3*2 + k4 ) /6.0;
		x += h;
	}
	
	return y;
}
double rk4( double (*f)(double,int), int n, double x_min, double x_max, int n_steps ) {
	// Fourth-order Runge-Kutta method for function only of x, with extra integer argument n which may act as an index.
	double h   = ( x_max - x_min ) / n_steps;
	double x   = x_min+h;
	double ret = 0.5 * ( f(x_min,n) + 4*f(x_min+0.5*h,n) + f(x_max,n) );
	
	for( int i=1; i<n_steps; i++ ){
		ret += f(x,n) + 2*f(x+0.5*h,n);
		x   += h;
	}
	
	return ret * h/3.0;
}
std::complex<double> rk4( std::complex<double> (*f)(std::complex<double>,int), int n, std::complex<double> x_min, std::complex<double> x_max, int n_steps ) {
	// Fourth-order Runge-Kutta method for function only of x, with extra integer argument n which may act as an index.
	std::complex<double> h   = ( x_max - x_min ) / ( (double) n_steps );
	std::complex<double> x   = x_min+h;
	std::complex<double> ret = 0.5 * ( f(x_min,n) + 4.0*f(x_min+0.5*h,n) + f(x_max,n) );
	
	for( int i=1; i<n_steps; i++ ){
		ret += f(x,n) + 2.0*f(x+0.5*h,n);
		x   += h;
	}
	
	return ret * h/3.0;
}
std::complex<double> rk4( std::complex<double> (*f)(double,int), int n, double x_min, double x_max, int n_steps ) {
	// Fourth-order Runge-Kutta method for complex-valued function only of x with x real, with extra integer argument n which may act as an index.
	double h   = ( x_max - x_min ) / n_steps;
	double x   = x_min+h;
	std::complex<double> ret = 0.5 * ( f(x_min,n) + 4.0*f(x_min+0.5*h,n) + f(x_max,n) );
	
	for( int i=1; i<n_steps; i++ ){
		ret += f(x,n) + 2.0*f(x+0.5*h,n);
		x   += h;
	}
	
	return ret * h/3.0;
}
double rk4_extra_argument( double (*f)(double,double), double p, double x_min, double x_max, int n_steps ) {
	// Fourth-order Runge-Kutta method for function only of x, with extra double argument p which may act as an index. A unique function name is required for this.
	double h   = ( x_max - x_min ) / n_steps;
	double x   = x_min+h;
	double ret = 0.5 * ( f(x_min,p) + 4*f(x_min+0.5*h,p) + f(x_max,p) );
	
	for( int i=1; i<n_steps; i++ ){
		ret += f(x,p) + 2*f(x+0.5*h,p);
		x   += h;
	}
	
	return ret * h/3.0;
}
double rk4( double (*f)(double,int,double,double), int n, int a, int b, double x_min, double x_max, int n_steps ) {
	// Fourth-order Runge-Kutta method for function only of x, with three extra arguments n,a,b.
	double h   = ( x_max - x_min ) / n_steps;
	double x   = x_min+h;
	double ret = 0.5 * ( f(x_min,n,a,b) + 4*f(x_min+0.5*h,n,a,b) + f(x_max,n,a,b) );
	
	for( int i=1; i<n_steps; i++ ){
		ret += f(x,n,a,b) + 2*f(x+0.5*h,n,a,b);
		x   += h;
	}
	
	return ret * h/3.0;
}
std::complex<double> rk4_extra_argument( std::complex<double> (*f)(double,double), double p, double x_min, double x_max, int n_steps ) {
	// Fourth-order Runge-Kutta method for complex-valuedfunction only of x, with extra double argument p which may act as an index. A unique function name is required for this.
	double h   = ( x_max - x_min ) / n_steps;
	double x   = x_min+h;
	std::complex<double> ret = 0.5 * ( f(x_min,p) + 4.0*f(x_min+0.5*h,p) + f(x_max,p) );
	
	for( int i=1; i<n_steps; i++ ){
		ret += f(x,p) + 2.0*f(x+0.5*h,p);
		x   += h;
	}
	
	return ret * h/3.0;
}
double rk4( double (*f)(double,int,int), int n1, int n2, double x_min, double x_max, int n_steps ) {
	// Fourth-order Runge-Kutta method for function only of x, with three two integer arguments n_i which may act as indices.
	double h   = ( x_max - x_min ) / n_steps;
	double x   = x_min+h;
	double ret = 0.5 * ( f(x_min,n1,n2) + 4*f(x_min+0.5*h,n1,n2) + f(x_max,n1,n2) );
	
	for( int i=1; i<n_steps; i++ ){
		ret += f(x,n1,n2) + 2*f(x+0.5*h,n1,n2);
		x   += h;
	}
	
	return ret * h/3.0;
}
double rk4( double (*f)(double,int,int,int), int n1, int n2, int n3, double x_min, double x_max, int n_steps ) {
	// Fourth-order Runge-Kutta method for function only of x, with three extra integer arguments n_i which may act as indices.
	double h   = ( x_max - x_min ) / n_steps;
	double x   = x_min+h;
	double ret = 0.5 * ( f(x_min,n1,n2,n3) + 4*f(x_min+0.5*h,n1,n2,n3) + f(x_max,n1,n2,n3) );
	
	for( int i=1; i<n_steps; i++ ){
		ret += f(x,n1,n2,n3) + 2*f(x+0.5*h,n1,n2,n3);
		x   += h;
	}
	
	return ret * h/3.0;
}
double rk4( double (*f)(double,int,int,double,double), int n1, int n2, double a3, double a4, double x_min, double x_max, int n_steps ) {
	// Fourth-order Runge-Kutta method for function only of x, with two extra integer arguments n_i and two extra double arguments a_i.
	double h   = ( x_max - x_min ) / n_steps;
	double x   = x_min+h;
	double ret = 0.5 * ( f(x_min,n1,n2,a3,a4) + 4*f(x_min+0.5*h,n1,n2,a3,a4) + f(x_max,n1,n2,a3,a4) );
	
	for( int i=1; i<n_steps; i++ ){
		ret += f(x,n1,n2,a3,a4) + 2*f(x+0.5*h,n1,n2,a3,a4);
		x   += h;
	}
	
	return ret * h/3.0;
}





std::vector<double> rk4_with_error( double (*f)(double), double x_min, double x_max, int n_steps ) {
	// Fourth-order Runge-Kutta method for function only of x.
	double h                 = ( x_max - x_min ) / n_steps;
	double x                 = x_min;
	double integral          = 0; 
	double integral_halfstep = 0;
	double error             = 0;
	
	for( int i=0; i<n_steps; i++ ){
		integral          += ( f(x) + 4*f(x+0.5*h) + f(x+h) ) * h/6.0;
		integral_halfstep += ( f(x) + 4*f(x+0.25*h) + 2*f(x+0.5*h) + 4*f(x+0.75*h) + f(x+h) ) * h/12.0;
		error             += pow( integral_halfstep - integral, 2 ) * 256/225.0;
		x                 += h;
	}
	
	return std::vector<double> { integral, sqrt( error ) };
}
std::vector<double> rk4_with_error( double (*f)(double, double), double x_min, double y_min, double x_max, int n_steps ) {
	// Fourth-order Runge-Kutta method for function of x and y.
	double h                 = ( x_max - x_min ) / n_steps;
	double x                 = x_min;
	double integral          = y_min;
	double integral_halfstep = y_min;
	double error             = 0;
	
	for( int i=0; i<n_steps; i++ ){
		
		double k1a = h * f( x        , integral           );
		double k2a = h * f( x + h*0.5, integral + k1a*0.5 );
		double k3a = h * f( x + h*0.5, integral + k2a*0.5 );
		double k4a = h * f( x + h    , integral + k3a     );
		integral += ( k1a + k2a*2 + k3a*2 + k4a ) /6.0;
		
		double k1b = h*0.5 * f( x         , integral_halfstep           );
		double k2b = h*0.5 * f( x + h*0.25, integral_halfstep + k1b*0.5 );
		double k3b = h*0.5 * f( x + h*0.25, integral_halfstep + k2b*0.5 );
		double k4b = h*0.5 * f( x + h*0.5 , integral_halfstep + k3b     );
		integral_halfstep += ( k1b + k2b*2 + k3b*2 + k4b ) /6.0;
		
		double k1c = h*0.5 * f( x + h*0.5 , integral_halfstep           );
		double k2c = h*0.5 * f( x + h*0.75, integral_halfstep + k1c*0.5 );
		double k3c = h*0.5 * f( x + h*0.75, integral_halfstep + k2c*0.5 );
		double k4c = h*0.5 * f( x + h     , integral_halfstep + k3c     );
		integral_halfstep += ( k1c + k2c*2 + k3c*2 + k4c ) /6.0;
		
		error             += pow( integral_halfstep - integral, 2 ) * 256/225.0;
		x                 += h;
	}
	
	return std::vector<double> { integral, sqrt( error ) };
}
std::vector<double> rk4_with_error( double (*f)(double,int), int n, double x_min, double x_max, int n_steps ) {
	// Fourth-order Runge-Kutta method for function only of x, with extra integer argument n.
	double h                 = ( x_max - x_min ) / n_steps;
	double x                 = x_min;
	double integral          = 0; 
	double integral_halfstep = 0;
	double error             = 0;
	
	for( int i=0; i<n_steps; i++ ){
		integral          += ( f(x,n) + 4*f(x+0.5*h,n) + f(x+h,n) ) * h/6.0;
		integral_halfstep += ( f(x,n) + 4*f(x+0.25*h,n) + 2*f(x+0.5*h,n) + 4*f(x+0.75*h,n) + f(x+h,n) ) * h/12.0;
		error             += pow( integral_halfstep - integral, 2 ) * 256/225.0;
		x                 += h;
	}
	
	return std::vector<double> { integral, sqrt( error ) };
}
std::vector<double> rk4_with_error_extra_argument( double (*f)(double,double), double p, double x_min, double x_max, int n_steps ) {
	// Fourth-order Runge-Kutta method for function only of x, with extra double argument p. A unique function name is required for this.
	double h                 = ( x_max - x_min ) / n_steps;
	double x                 = x_min;
	double integral          = 0; 
	double integral_halfstep = 0;
	double error             = 0;
	
	for( int i=0; i<n_steps; i++ ){
		integral          += ( f(x,p) + 4*f(x+0.5*h,p) + f(x+h,p) ) * h/6.0;
		integral_halfstep += ( f(x,p) + 4*f(x+0.25*h,p) + 2*f(x+0.5*h,p) + 4*f(x+0.75*h,p) + f(x+h,p) ) * h/12.0;
		error             += pow( integral_halfstep - integral, 2 ) * 256/225.0;
		x                 += h;
	}
	
	return std::vector<double> { integral, sqrt( error ) };
}




std::vector<double> rkf( double (*f)(double), double x_min, double x_max, int n_steps ) {
	// Runge-Kutta-Felhberg method for function of x only.
	double h     = ( x_max - x_min ) / n_steps;
	double x     = x_min;
	double y     = 0;
	double error = 0;
	
	for( int i=0; i<n_steps; i++ ){
		double k1 = h * f( x             );
		double k3 = h * f( x + h*0.375   );	// k2 not needed.
		double k4 = h * f( x + h*12/13.0 );
		double k5 = h * f( x + h         );
		double k6 = h * f( x + h*0.5     );
		y     += k1*16/135.0 + k3*6656/12825.0 + k4*28561/56430.0 - k5*0.18 + k6*2/55.0;
		error += pow( k1/360.0 - k3*128/4375.0 - k4*2197/75240.0 + k5*0.02 + k6*2/55.0, 2 );
		x     += h;
	}
	
	return std::vector<double> { y, sqrt( error ) };
}
std::vector<double> rkf( double (*f)(double,double), double x_min, double y_min, double x_max, int n_steps ) {
	// Runge-Kutta-Felhberg method for function of x and y.
	double h      = ( x_max - x_min ) / n_steps;
	double x      = x_min;
	double y      = y_min;
	double error  = 0;
	
	for( int i=0; i<n_steps; i++ ){
		double k1 = h * f( x            , y                                                                 );
		double k2 = h * f( x + h*0.25   , y + k1*0.25                                                       );
		double k3 = h * f( x + h*0.375  , y + ( k1*3 + k2*9 ) / 32.0                                        );
		double k4 = h * f( x + h*12/13.0, y + ( k1*1932 - k2*7200 + k3*7296 ) / 2197.0                      );
		double k5 = h * f( x + h        , y + k1*439/216.0 - k2*8 + k3*3680/513.0 - k4*845/4104.0           );
		double k6 = h * f( x + h*0.5    , y - k1*8/27.0 + k2*2 - k3*3544/2565.0 + k4*1859/4104.0 - k5*0.275 );
		y     += k1*16/135.0 + k3*6656/12825.0 + k4*28561/56430.0 - k5*0.18 + k6*2/55.0;
		error += pow( k1/360.0 - k3*128/4375.0 - k4*2197/75240.0 + k5*0.02 + k6*2/55.0, 2 );
		x     += h;
	}
	
	return std::vector<double> { y, sqrt( error ) };
}
std::vector<double> rkf( double (*f)(double,int), int n, double x_min, double x_max, int n_steps ) {
	// Runge-Kutta-Felhberg method for function of x only.
	double h     = ( x_max - x_min ) / n_steps;
	double x     = x_min;
	double y     = 0;
	double error = 0;
	
	for( int i=0; i<n_steps; i++ ){
		double k1 = h * f( x            , n );
		double k3 = h * f( x + h*0.375  , n );
		double k4 = h * f( x + h*12/13.0, n );
		double k5 = h * f( x + h        , n );
		double k6 = h * f( x + h*0.5    , n );
		y     += k1*16/135.0 + k3*6656/12825.0 + k4*28561/56430.0 - k5*0.18 + k6*2/55.0;
		error += pow( k1/360.0 - k3*128/4375.0 - k4*2197/75240.0 + k5*0.02 + k6*2/55.0, 2 );
		x     += h;
	}
	
	return std::vector<double> { y, sqrt( error ) };
}
std::vector<double> rkf_extra_argument( double (*f)(double,double), double p, double x_min, double x_max, int n_steps ) {
	// Runge-Kutta-Felhberg method for function of x only, with extra double argument p. A unique function name is required for this.
	double h     = ( x_max - x_min ) / n_steps;
	double x     = x_min;
	double y     = 0;
	double error = 0;
	
	for( int i=0; i<n_steps; i++ ){
		double k1 = h * f( x            , p );
		double k3 = h * f( x + h*0.375  , p );
		double k4 = h * f( x + h*12/13.0, p );
		double k5 = h * f( x + h        , p );
		double k6 = h * f( x + h*0.5    , p );
		y     += k1*16/135.0 + k3*6656/12825.0 + k4*28561/56430.0 - k5*0.18 + k6*2/55.0;
		error += pow( k1/360.0 - k3*128/4375.0 - k4*2197/75240.0 + k5*0.02 + k6*2/55.0, 2 );
		x     += h;
	}
	
	return std::vector<double> { y, sqrt( error ) };
}




double trapezium( double (*f)(double), double x_min, double x_max, int n_steps ) {
	// Trapezium rule for function only of x.
	double h   = ( x_max - x_min ) / n_steps;
	double ret = 0.5 * ( f(x_min) + f(x_max) );
	
	for( int i=1; i<=n_steps-1; i++ ){
		ret += f( x_min + i * h );
	}
	
	return ret * h;
}
double trapezium( double (*f)(double,int), int n, double x_min, double x_max, int n_steps ) {
	// Trapezium rule for function only of x, with extra integer argument n.
	double h   = ( x_max - x_min ) / n_steps;
	double ret = 0.5 * ( f(x_min,n) + f(x_max,n) );
	
	for( int i=1; i<=n_steps-1; i++ ){
		ret += f( x_min + i * h, n );
	}
	
	return ret * h;
}
double trapezium_extra_argument( double (*f)(double,double), double p, double x_min, double x_max, int n_steps ) {
	// Trapezium rule for function only of x, with extra double argument p. A unique function name is required for this.
	double h   = ( x_max - x_min ) / n_steps;
	double ret = 0.5 * ( f(x_min,p) + f(x_max,p) );
	
	for( int i=1; i<=n_steps-1; i++ ){
		ret += f( x_min + i * h, p );
	}
	
	return ret * h;
}

double euler( double (*f)(double), double x_min, double x_max, int n_steps ) {
	// Euler method for function only of x.
	double h   = ( x_max - x_min ) / n_steps;
	double ret = 0;
	
	for( int i=0; i<n_steps; i++ ){
		ret += f( x_min + i * h );
	}
	return ret * h;
}
double trapezium( std::vector<double> f, std::vector<double> x ) {
	// Trapezium rule for function only of x.
	double ret = f.back()*x.back() - f[0]*x[0];
	
	for( int i=0; i<x.size()-1; i++ ){
		ret += f[i]*x[i+1] - f[i+1]*x[i];
	}
	
	return ret * 0.5;
}
double trapezium_minus( std::vector<double> f, std::vector<double> x ) {
	// Trapezium rule for function only of x. This will be used when a substitution introduces a minus sign.
	double ret = f[0]*x[0] - f.back()*x.back();
	
	for( int i=0; i<x.size()-1; i++ ){
		ret += f[i+1]*x[i] - f[i]*x[i+1];
	}
	
	return ret * 0.5;
}





std::vector<double> euler_with_error( double (*f)(double), double x_min, double x_max, int n_steps ) {
	// Fourth-order Runge-Kutta method for function only of x.
	double h                 = ( x_max - x_min ) / n_steps;
	double x                 = x_min;
	double integral          = 0; 
	double integral_halfstep = 0;
	double error             = 0;
	
	for( int i=0; i<n_steps; i++ ){
		integral          += f(x) * h;
		integral_halfstep += ( f(x) + f(x+0.5*h) ) * 0.5*h;
		error             += pow( integral_halfstep - integral, 2 ) * 16/9.0;
		x                 += h;
	}
	
	return std::vector<double> { integral, sqrt( error ) };
}




//-------------------- TWO-DIMENSIONAL INTEGRALS --------------------
double integral_2d( std::vector< std::vector<double> > f, std::vector< std::vector<double> > du ){
	// Integral of a 2D function given as a list, with the area elements already calculated.
	// Uses Notes/Numerical Integration Proposition 4.6.
	int    N1  = du.size();
	int    N2  = du[0].size();
	double ret = 0;
	
	for( int i=0; i<N1-1; i++ ){
		for( int j=0; j<N2-1; j++ ){
			ret += ( f[i][j] + f[i][j+1] + f[i+1][j] + f[i+1][j+1] ) * du[i][j];
		}
	}
	return ret * 0.25;
}



double integral_2d( std::vector< std::vector<double> > f, std::vector<double> du1, std::vector<double> du2 ){
	// Integral of a 2D function given as a list, with the area elements already calculated.
	// More computationally efficient than, but equivalent to, integral_2d( std::vector< std::vector<double> > f, std::vector< std::vector<double> > du ).
	int    N1  = du1.size();
	int    N2  = du2.size();
	double ret = 0;
	
	for( int i=1; i<N1-1; i++ ){
		ret += ( f[i][0] * du2[0] + f[i].back() * du2[N2-2] ) * ( du1[i-1] + du1[i] );
		for( int j=1; j<N2-1; j++ ){
			ret += f[i][j] * ( ( du1[i-1] + du1[i] ) * ( du2[j-1] + du2[j] ) );
		}
	}
	for( int j=1; j<N2-1; j++ ){
		ret += ( f[0][j] * du1[0] + f.back()[j] * du1[N1-2] ) * ( du2[j-1] + du2[j] );
	}
	ret += ( f[0][0] * du2[0] + f[0].back() * du2[N2-2] ) * du1[0] + ( f.back()[0] * du2[0] + f.back().back() * du2[N2-2] ) * du1[N1-2];
	
	return ret * 0.25;
}
std::complex<double> integral_2d( std::vector< std::vector< std::complex<double> > > f, std::vector<double> du1, std::vector<double> du2 ){
	// Integral of a 2D function given as a list, with the area elements already calculated.
	int    N1  = du1.size();
	int    N2  = du2.size();
	std::complex<double> ret = 0;
	
	for( int i=1; i<N1-1; i++ ){
		ret += ( f[i][0] * du2[0] + f[i].back() * du2[N2-2] ) * ( du1[i-1] + du1[i] );
		for( int j=1; j<N2-1; j++ ){
			ret += f[i][j] * ( ( du1[i-1] + du1[i] ) * ( du2[j-1] + du2[j] ) );
		}
	}
	for( int j=1; j<N2-1; j++ ){
		ret += ( f[0][j] * du1[0] + f.back()[j] * du1[N1-2] ) * ( du2[j-1] + du2[j] );
	}
	ret += ( f[0][0] * du2[0] + f[0].back() * du2[N2-2] ) * du1[0] + ( f.back()[0] * du2[0] + f.back().back() * du2[N2-2] ) * du1[N1-2];
	
	return ret * 0.25;
}




double integral_cartesian( double (*f)(double,double), double x_min, double x_max, double y_min, double y_max, int nx, int ny ){
	// Multi-step finite-area method for surface integral in Cartesian coordinates (x,y). Function is real-valued.
	
	//----- Calculate subintervals in each dimension -----
	double hx  = ( x_max - x_min ) / ( (double) nx );
	double hy  = ( y_max - y_min ) / ( (double) ny );
	
	//----- Perform the integral -----
	double ret = 0;
	double x   = x_min + hx;
	double y   = y_min + hy;
	
	//--- Double sum, and single sum over x ---
	for( int i=1; i<nx; i++ ) {
		y = y_min + hy;
		ret += 0.5 * ( f( x, y_min ) + f( x, y_max ) );	// Single sum over x
		for( int j=1; j<ny; j++ ) {
			ret += f( x, y );							// Double sum
			y   += hy;
		}
		x += hx;
	}
	
	//--- Single sum over y ---
	y = y_min + hy;
	for( int i=1; i<ny; i++ ) {
		ret += 0.5 * ( f( x_min, y ) + f( x_max, y ) );
		y   += hy;
	}
	
	//--- Vertices ---
	ret += 0.25 * ( f( x_min, y_min ) + f( x_min, y_max ) + f( x_max, y_min ) + f( x_max, y_max ) );
	
	return ret * hx * hy;
}
double integral_cartesian( double (*f)(double,double,int,int), int a, int b, double x_min, double x_max, double y_min, double y_max, int nx, int ny ){
	// Multi-step finite-area method for surface integral in Cartesian coordinates (x,y). Function is real-valued and can take two integer parameters a and b.
	
	//----- Calculate subintervals in each dimension -----
	double hx  = ( x_max - x_min ) / ( (double) nx );
	double hy  = ( y_max - y_min ) / ( (double) ny );
	
	//----- Perform the integral -----
	double ret = 0;
	double x   = x_min + hx;
	double y   = y_min + hy;
	
	//--- Double sum, and single sum over x ---
	for( int i=1; i<nx; i++ ) {
		y = y_min + hy;
		ret += 0.5 * ( f( x, y_min, a, b ) + f( x, y_max, a, b ) );	// Single sum over x
		for( int j=1; j<ny; j++ ) {
			ret += f( x, y, a, b );									// Double sum
			y   += hy;
		}
		x += hx;
	}
	
	//--- Single sum over y ---
	y = y_min + hy;
	for( int i=1; i<ny; i++ ) {
		ret += 0.5 * ( f( x_min, y, a, b ) + f( x_max, y, a, b ) );
		y   += hy;
	}
	
	//--- Vertices ---
	ret += 0.25 * ( f( x_min, y_min, a, b ) + f( x_min, y_max, a, b ) + f( x_max, y_min, a, b ) + f( x_max, y_max, a, b ) );
	
	return ret * hx * hy;
}




double integral_spherical_surface( std::vector< std::vector<double> > f, std::vector<double> t, std::vector<double> p ){
	// In this version, f is a list, not a function, and theta,phi are lists with arbitrary spacing rules, not simply linear.
	const double pi = 3.14159265358979323846;
	int    N1  = t.size();
	int    N2  = p.size();
	double ret = 0;			// Final value to be returned.
	
	for( int i=1; i<N1-1; i++ ){
		for( int j=1; j<N2-1; j++ ){
			ret += f[i][j] * (cos(t[i])-cos(t[i+1]))*(p[j+1]-p[j]);
		}
	}
	for( int i=1; i<N1-1; i++ ){
		ret += 0.5*(cos(t[i])-cos(t[i+1])) * ( f[i][0] * p[1] + f[i].back() * (2*pi-p[N2-2]) );
	}
	for( int j=1; j<N2-1; j++ ){
		ret += 0.5*(p[j+1]-p[j]) * ( f[0][j] * (1-cos(t[1])) + f.back()[j] * (cos(t[N1-2])+1) );
	}
	ret += 0.25 * ( (1-cos(t[1])) * ( f[0][0] * p[1] + f[0].back() * (2*pi-p[N2-2]) ) + (cos(t[N1-2])+1) * ( f.back()[0] * p[1] + f.back().back() * (2*pi-p[N2-2]) ) );
	
	return ret;
}
std::complex<double> integral_spherical_surface( std::vector< std::vector< std::complex<double> > > f, std::vector<double> t, std::vector<double> p ){
	// In this version, f is a list, not a function, and theta,phi are lists with arbitrary spacing rules, not simply linear.
	// The range of theta MUST be from 0 to pi, and the range of phi MUST be from 0 to 2pi.
	const double pi = 3.14159265358979323846;
	int    N1  = t.size();
	int    N2  = p.size();
	std::complex<double> ret = 0;			// Final value to be returned.
	
	for( int i=1; i<N1-1; i++ ){
		for( int j=1; j<N2-1; j++ ){
			ret += f[i][j] * (cos(t[i])-cos(t[i+1]))*(p[j+1]-p[j]);
		}
	}
	for( int i=1; i<N1-1; i++ ){
		ret += 0.5*(cos(t[i])-cos(t[i+1])) * ( f[i][0] * p[1] + f[i].back() * (2*pi-p[N2-2]) );
	}
	for( int j=1; j<N2-1; j++ ){
		ret += 0.5*(p[j+1]-p[j]) * ( f[0][j] * (1-cos(t[1])) + f.back()[j] * (cos(t[N1-2])+1) );
	}
	ret += 0.25 * ( (1-cos(t[1])) * ( f[0][0] * p[1] + f[0].back() * (2*pi-p[N2-2]) ) + (cos(t[N1-2])+1) * ( f.back()[0] * p[1] + f.back().back() * (2*pi-p[N2-2]) ) );
	
	return ret;
}
std::complex<double> integral_spherical_surface_using_references( std::vector< std::vector< std::complex<double> > >& f, std::vector<double>& t, std::vector<double>& p ){
	// In this version, f is a list, not a function, and theta,phi are lists with arbitrary spacing rules, not simply linear.
	// The range of theta MUST be from 0 to pi, and the range of phi MUST be from 0 to 2pi.
	const double pi = 3.14159265358979323846;
	int    N1  = t.size();
	int    N2  = p.size();
	std::complex<double> ret = 0;			// Final value to be returned.
	
	for( int i=1; i<N1-1; i++ ){
		for( int j=1; j<N2-1; j++ ){
			ret += f[i][j] * (cos(t[i])-cos(t[i+1]))*(p[j+1]-p[j]);
		}
	}
	for( int i=1; i<N1-1; i++ ){
		ret += 0.5*(cos(t[i])-cos(t[i+1])) * ( f[i][0] * p[1] + f[i].back() * (2*pi-p[N2-2]) );
	}
	for( int j=1; j<N2-1; j++ ){
		ret += 0.5*(p[j+1]-p[j]) * ( f[0][j] * (1-cos(t[1])) + f.back()[j] * (cos(t[N1-2])+1) );
	}
	ret += 0.25 * ( (1-cos(t[1])) * ( f[0][0] * p[1] + f[0].back() * (2*pi-p[N2-2]) ) + (cos(t[N1-2])+1) * ( f.back()[0] * p[1] + f.back().back() * (2*pi-p[N2-2]) ) );
	
	return ret;
}
		

double integral_spherical_surface( double (*f)(double,double), double t_min, double t_max, double p_min, double p_max, double n_steps_t, double n_steps_p ) {
	// THESE ARE OUTDATED AND NEED TO BE REPLACED BY THE METHOD ABOVE, MODIFIED TO DEAL WITH FUNCTIONS INSTEAD OF LISTS.
	// Multi-step finite-area method for surface integral in spherical polar coordinates (t,p)=(theta,phi). Function is real-valued.
	double h_t = ( t_max - t_min ) / n_steps_t;
	double h_p = ( p_max - p_min ) / n_steps_p;
	double t   = t_min;
	double p   = p_min;
	double dA  = 0;
	double ret = 0;
	
	for( int i=0; i<n_steps_t; i++ ) {
		p  = p_min;
		dA = h_p * ( cos(t) * ( 1 - cos(h_t) ) + sin(t)*sin(h_t) );		// Constant for given phi.
		for( int j=0; j<n_steps_p; j++ ) {
			ret += 0.25 * dA * ( f( t, p ) + f( t+h_t, p ) + f( t, p+h_p ) + f( t+h_t, p+h_p ) );
			p   += h_p;
		}
		t += h_t;
	}
	
	return ret;
}
std::complex<double> integral_spherical_surface( std::complex<double> (*f)(double,double), double t_min, double t_max, double p_min, double p_max, double n_steps_t, double n_steps_p ) {
	// Multi-step finite-area method for surface integral in spherical polar coordinates (t,p)=(theta,phi). Function is complex.
	double h_t               = ( t_max - t_min ) / n_steps_t;
	double h_p               = ( p_max - p_min ) / n_steps_p;
	double t                 = t_min;
	double p                 = p_min;
	double dA                = 0;
	std::complex<double> ret = 0;
	
	for( int i=0; i<n_steps_t; i++ ) {
		p  = p_min;
		dA = h_p * ( cos(t) * ( 1 - cos(h_t) ) + sin(t)*sin(h_t) );		// Constant for given phi.
		for( int j=0; j<n_steps_p; j++ ) {
			ret += 0.25 * dA * ( f( t, p ) + f( t+h_t, p ) + f( t, p+h_p ) + f( t+h_t, p+h_p ) );
			p   += h_p;
		}
		t += h_t;
	}
	
	return ret;
}
std::complex<double> integral_spherical_surface( std::complex<double> (*f)(double,double,int,int), double t_min, double t_max, double p_min, double p_max, double n_steps_t, double n_steps_p, int ell, int m ) {
	// Multi-step finite-area method for surface integral in spherical polar coordinates (t,p)=(theta,phi). Function is complex and may take two extra integer arguments. These will most commonly be the spherical harmonic indices ell, m and are labelled as such.
	double h_t               = ( t_max - t_min ) / n_steps_t;
	double h_p               = ( p_max - p_min ) / n_steps_p;
	double t                 = t_min;
	double p                 = p_min;
	double dA                = 0;
	std::complex<double> ret = 0;
	
	for( int i=0; i<n_steps_t; i++ ) {
		p  = p_min;
		dA = h_p * ( cos(t) * ( 1 - cos(h_t) ) + sin(t)*sin(h_t) );		// Constant for given phi.
		for( int j=0; j<n_steps_p; j++ ) {
			ret += 0.25 * dA * ( f( t, p, ell, m ) + f( t+h_t, p, ell, m ) + f( t, p+h_p, ell, m ) + f( t+h_t, p+h_p, ell, m ) );
			p += h_p;
		}
		t += h_t;
	}
	
	return ret;
}		



//-------------------- THREE-DIMENSIONAL INTEGRALS --------------------
double integral_3d_cartesian( double (*f)(double,double,double), double x_min, double x_max, double y_min, double y_max, double z_min, double z_max,
                              int nx, int ny, int nz ) {
	// Multi-step finite-volume method for volume integral in Cartesian coordinates (x,y,z).
	
	//----- Calculate subintervals in each dimension -----
	double hx  = ( x_max - x_min ) / ( (double) nx );
	double hy  = ( y_max - y_min ) / ( (double) ny );
	double hz  = ( z_max - z_min ) / ( (double) nz );
	
	//----- Perform the integral -----
	double ret = 0;
	double x   = x_min + hx;
	double y   = y_min + hy;
	double z   = z_min + hz;
	
	//--- Triple sum; Double sum over x,y; Double sum over x,z; Single sum over x ---
	for( int i=1; i<nx; i++ ) {
		y = y_min + hy;
		for( int j=1; j<ny; j++ ) {
			z = z_min + hz;
			ret += 0.5 * ( f( x, y, z_min ) + f( x, y, z_max ) );														// Double sum over x,y
			for( int k=1; k<nz; k++ ) {
				ret += f( x, y , z );																					// Triple sum
				z   += hz;
			}
			y += hy;
		}
		z = z_min + hz;
		for( int k=1; k<nz; k++ ) {
			ret += 0.5 * ( f( x, y_min, z ) + f( x, y_max, z ) );														// Double sum over x,z
			z   += hz;
		}
		ret += 0.25 * ( f( x, y_min, z_min ) + f( x, y_min, z_max ) + f( x, y_max, z_min ) + f( x, y_max, z_max ) );	// Single sum over x
		x += hx;
	}
	
	//--- Double sum over y,z; Single sum over z ---	
	y = y_min + hy;
	for( int j=1; j<ny; j++ ) {
		z = z_min + hz;
		for( int k=1; k<nz; k++ ) {
			ret += 0.5 * ( f( x_min, y, z ) + f( x_max, y, z ) );														// Double sum over y,z
			z   += hz;
		}
		ret += 0.25 * ( f( x_min, y, z_min ) + f( x_min, y, z_max ) + f( x_max, y, z_min ) + f( x_max, y, z_max ) );	// Single sum over y
		y += hy;
	}
	
	//--- Single sum over z ---
	z = z_min + hz;
	for( int k=1; k<nz; k++ ) {
		ret += 0.25 * ( f( x_min, y_min, z ) + f( x_min, y_max, z ) + f( x_max, y_min, z ) + f( x_max, y_max, z ) );
		z   += hz;
	}
	
	//--- Vertices ---
	ret += 0.125 * (   f( x_min, y_min, z_min ) + f( x_min, y_min, z_max ) + f( x_min, y_max, z_min ) + f( x_min, y_max, z_max )
	                 + f( x_max, y_min, z_min ) + f( x_max, y_min, z_max ) + f( x_max, y_max, z_min ) + f( x_max, y_max, z_max ) );
	
	return ret * hx * hy * hz;
}





double integral_3d_spherical( double (*f)(double,double,double), double r_min, double r_max, double t_min, double t_max, double p_min, double p_max, int nr, int nt, int np ){
	// Multi-step finite-volume method for volume integral in spherical polar coordinates (r,t,p)=(r,theta,phi).
	
	//----- Calculate subintervals in each dimension, and commonly used volume elements -----
	double hr          = ( r_max - r_min ) / ( (double) nr );
	double ht          = ( t_max - t_min ) / ( (double) nt );
	double hp          = ( p_max - p_min ) / ( (double) np );
	double r_max_m1    = r_min + (nr-1) * hr;
	double t_max_m1    = t_min + (nt-1) * ht;
	double dV_r_min    = r_min   *r_min   *hr + r_min   *hr*hr + hr*hr*hr/3.0;
	double dV_r_max_m1 = r_max_m1*r_max_m1*hr + r_max_m1*hr*hr + hr*hr*hr/3.0;
	double dV_t_min    = cos(t_min   ) - cos(t_min+ht   );
	double dV_t_max_m1 = cos(t_max_m1) - cos(t_max_m1+ht);
	
	//----- Perform the integral -----
	double ret = 0;
	double r   = r_min + hr;
	double t   = t_min + ht;
	double p   = p_min + hp;
	
	//--- Triple sum; Double sum over r,t; Double sum over r,p; Single sum over r ---
	for( int i=1; i<nr; i++ ){
		t = t_min + ht;
		double dV_r = 2 * hr * ( r*r + hr*hr/3.0 );
		for( int j=1; j<nt; j++ ){
			p = p_min + hp;
			double dV = 2 * dV_r * sin(t) * sin(ht);
			ret += ( f(r,t,p_min) + f(r,t,p_max) ) * dV;																				// Double sum over i,j
			for( int k=1; k<np; k++ ){
				ret += f(r,t,p) * 2 * dV;																								// Triple sum
				p   += hp;
			}
			t += ht;
		}
		p = p_min + hp;
		for( int k=1; k<np; k++ ){
			ret += ( f(r,t_min,p) * dV_t_min + f(r,t_max,p) * dV_t_max_m1 ) * 2 * dV_r;													// Double sum over i,k
			p   += hp;
		}
		ret += ( ( f(r,t_min,p_min) + f(r,t_min,p_max) ) * dV_t_min + ( f(r,t_max,p_min) + f(r,t_max,p_max) ) * dV_t_max_m1 ) * dV_r;	// Single sum over i
		r += hr;
	}
	
	//--- Double sum over t,p; single sum over t ---
	t = t_min + ht;
	for( int j=1; j<nt; j++ ){
		p = p_min + hp;
		double dV_t = 2 * sin(t) * sin(ht);
		for( int k=1; k<np; k++ ){
			ret += ( f(r_min,t,p) * dV_r_min + f(r_max,t,p) * dV_r_max_m1 ) * 2 * dV_t;													// Double sum over j,k
			p   += hp;
		}
		ret += ( ( f(r_min,t,p_min) + f(r_min,t,p_max) ) * dV_r_min + ( f(r_max,t,p_min) + f(r_max,t,p_max) ) * dV_r_max_m1 ) * dV_t;	// Single sum over j
		t   += ht;
	}
	
	//--- Single sum over p ---
	p = p_min + hp;
	for( int k=1; k<np; k++ ){
		ret += 2 * ( f(r_min,t_min,p)*dV_t_min + f(r_min,t_max,p)*dV_t_max_m1 ) * dV_r_min;
		ret += 2 * ( f(r_max,t_min,p)*dV_t_min + f(r_max,t_max,p)*dV_t_max_m1 ) * dV_r_max_m1;
		p   += hp;
	}
	
	//--- Coordinate endpoints ---
	ret += ( ( f(r_min,t_min,p_min) + f(r_min,t_min,p_max) )*dV_t_min + ( f(r_min,t_max,p_min) + f(r_min,t_max,p_max) )*dV_t_max_m1 ) * dV_r_min;
	ret += ( ( f(r_max,t_min,p_min) + f(r_max,t_min,p_max) )*dV_t_min + ( f(r_max,t_max,p_min) + f(r_max,t_max,p_max) )*dV_t_max_m1 ) * dV_r_max_m1;
	
	return ret * hp * 0.125;
}
double integral_3d_spherical( std::vector< std::vector< std::vector<double> > > f, std::vector<double> r, std::vector<double> t, std::vector<double> p ){
	// Multi-step finite-volume method for volume integral in spherical polar coordinates (r,t,p)=(r,theta,phi).
	// In this version, f is a real-valued 3D list instead of a function, and r,t,p are lists with arbitrary and variable stepsize.
	int    N1 = r.size();
	int    N2 = t.size();
	int    N3 = p.size();
	double ret    = 0;		// Final value to be returned.
	
	//----- Triple sum -----
	for( int i=1; i<N1-1; i++ ){
		for( int j=1; j<N2-1; j++ ){
			for( int k=1; k<N3-1; k++ ){
				ret += f[i][j][k] * (pow(r[i+1],3)-pow(r[i],3))*(cos(t[j])-cos(t[j+1]))*(p[k+1]-p[k])/3.0;
			}
		}
	}
	
	//----- Double sums -----
	for( int i=1; i<N1-1; i++ ){
		for( int j=1; j<N2-1; j++ ){
			ret += 1.0/6 * (pow(r[i+1],3)-pow(r[i],3))*(cos(t[j])-cos(t[j+1])) * ( f[i][j][0] * (p[1]-p[0]) + f[i][j].back() * (p.back()-p[N3-2]) );
		}
	}
	
	for( int i=1; i<N1-1; i++ ){
		for( int k=1; k<N3-1; k++ ){
			ret += 1.0/6 * (pow(r[i+1],3)-pow(r[i],3))*(p[k+1]-p[k]) * ( f[i][0][k] * (cos(t[0])-cos(t[1])) + f[i].back()[k] * (cos(t[N2-2])-cos(t.back())) );
		}
	}
	
	for( int j=1; j<N2-1; j++ ){
		for( int k=1; k<N3-1; k++ ){
			ret += 1.0/6 * (cos(t[j])-cos(t[j+1]))*(p[k+1]-p[k]) * ( f[0][j][k] * (pow(r[1],3)-pow(r[0],3)) + f.back()[j][k] * (pow(r.back(),3)-pow(r[N1-2],3)) );
		}
	}
	
	//----- Single sums -----
	for( int i=1; i<N1-1; i++ ){
		ret += 1.0/12 * (pow(r[i+1],3)-pow(r[i],3)) * ( (cos(t[0])-cos(t[1])) * ( f[i][0][0]*(p[1]-p[0]) + f[i][0].back()*(p.back()-p[N3-2]) ) + (cos(t[N2-2])-cos(t.back())) * ( f[i].back()[0]*(p[1]-p[0]) + f[i].back().back()*(p.back()-p[N3-2]) ) );
	}
	
	for( int j=1; j<N2-1; j++ ){
		ret += 1.0/12 * (cos(t[j])-cos(t[j+1])) * ( (pow(r[1],3)-pow(r[0],3)) * ( f[0][j][0]*(p[1]-p[0]) + f[0][j].back()*(p.back()-p[N3-2]) ) + (pow(r.back(),3)-pow(r[N1-2],3)) * ( f.back()[j][0]*(p[1]-p[0]) + f.back()[j].back()*(p.back()-p[N3-2]) ) );
	}
	
	for( int k=1; k<N3-1; k++ ){
		ret += 1.0/12 * (p[k+1]-p[k]) * ( (pow(r[1],3)-pow(r[0],3)) * ( f[0][0][k]*(cos(t[0])-cos(t[1])) + f[0].back()[k]*(cos(t[N2-2])-cos(t.back())) ) + (pow(r.back(),3)-pow(r[N1-2],3)) * ( f.back()[0][k]*(cos(t[0])-cos(t[1])) + f.back().back()[k]*(cos(t[N2-2])-cos(t.back())) ) );
	}
	
	//----- Vertices -----
	ret += 1.0/24 * ( (pow(r[1]    ,3)-pow(r[0]   ,3)) * ( (cos(t[0])-cos(t[1])) * ( f[0    ][0][0]*(p[1]-p[0]) + f[0    ][0].back()*(p.back()-p[N3-2]) ) + (cos(t[N2-2])-cos(t.back())) * ( f[0    ].back()[0]*(p[1]-p[0]) + f[0    ].back().back()*(p.back()-p[N3-2]) ) )
	                + (pow(r.back(),3)-pow(r[N1-2],3)) * ( (cos(t[0])-cos(t[1])) * ( f.back()[0][0]*(p[1]-p[0]) + f.back()[0].back()*(p.back()-p[N3-2]) ) + (cos(t[N2-2])-cos(t.back())) * ( f.back().back()[0]*(p[1]-p[0]) + f.back().back().back()*(p.back()-p[N3-2]) ) ) );
	
	return ret;
}
std::complex<double> integral_3d_spherical( std::vector< std::vector< std::vector< std::complex<double> > > > f, std::vector<double> r, std::vector<double> t, std::vector<double> p ){
	// Multi-step finite-volume method for volume integral in spherical polar coordinates (r,t,p)=(r,theta,phi).
	// In this version, f is a complex-valued 3D list instead of a function, and r,t,p are lists with arbitrary and variable stepsize.
	int    N1 = r.size();
	int    N2 = t.size();
	int    N3 = p.size();
	std::complex<double> ret    = 0;		// Final value to be returned.
	
	//----- Triple sum -----
	for( int i=1; i<N1-1; i++ ){
		for( int j=1; j<N2-1; j++ ){
			for( int k=1; k<N3-1; k++ ){
				ret += f[i][j][k] * (pow(r[i+1],3)-pow(r[i],3))*(cos(t[j])-cos(t[j+1]))*(p[k+1]-p[k])/3.0;
			}
		}
	}
	
	//----- Double sums -----
	for( int i=1; i<N1-1; i++ ){
		for( int j=1; j<N2-1; j++ ){
			ret += 1.0/6 * (pow(r[i+1],3)-pow(r[i],3))*(cos(t[j])-cos(t[j+1])) * ( f[i][j][0] * (p[1]-p[0]) + f[i][j].back() * (p.back()-p[N3-2]) );
		}
	}
	
	for( int i=1; i<N1-1; i++ ){
		for( int k=1; k<N3-1; k++ ){
			ret += 1.0/6 * (pow(r[i+1],3)-pow(r[i],3))*(p[k+1]-p[k]) * ( f[i][0][k] * (cos(t[0])-cos(t[1])) + f[i].back()[k] * (cos(t[N2-2])-cos(t.back())) );
		}
	}
	
	for( int j=1; j<N2-1; j++ ){
		for( int k=1; k<N3-1; k++ ){
			ret += 1.0/6 * (cos(t[j])-cos(t[j+1]))*(p[k+1]-p[k]) * ( f[0][j][k] * (pow(r[1],3)-pow(r[0],3)) + f.back()[j][k] * (pow(r.back(),3)-pow(r[N1-2],3)) );
		}
	}
	
	//----- Single sums -----
	for( int i=1; i<N1-1; i++ ){
		ret += 1.0/12 * (pow(r[i+1],3)-pow(r[i],3)) * ( (cos(t[0])-cos(t[1])) * ( f[i][0][0]*(p[1]-p[0]) + f[i][0].back()*(p.back()-p[N3-2]) ) + (cos(t[N2-2])-cos(t.back())) * ( f[i].back()[0]*(p[1]-p[0]) + f[i].back().back()*(p.back()-p[N3-2]) ) );
	}
	
	for( int j=1; j<N2-1; j++ ){
		ret += 1.0/12 * (cos(t[j])-cos(t[j+1])) * ( (pow(r[1],3)-pow(r[0],3)) * ( f[0][j][0]*(p[1]-p[0]) + f[0][j].back()*(p.back()-p[N3-2]) ) + (pow(r.back(),3)-pow(r[N1-2],3)) * ( f.back()[j][0]*(p[1]-p[0]) + f.back()[j].back()*(p.back()-p[N3-2]) ) );
	}
	
	for( int k=1; k<N3-1; k++ ){
		ret += 1.0/12 * (p[k+1]-p[k]) * ( (pow(r[1],3)-pow(r[0],3)) * ( f[0][0][k]*(cos(t[0])-cos(t[1])) + f[0].back()[k]*(cos(t[N2-2])-cos(t.back())) ) + (pow(r.back(),3)-pow(r[N1-2],3)) * ( f.back()[0][k]*(cos(t[0])-cos(t[1])) + f.back().back()[k]*(cos(t[N2-2])-cos(t.back())) ) );
	}
	
	//----- Vertices -----
	ret += 1.0/24 * ( (pow(r[1]    ,3)-pow(r[0]   ,3)) * ( (cos(t[0])-cos(t[1])) * ( f[0    ][0][0]*(p[1]-p[0]) + f[0    ][0].back()*(p.back()-p[N3-2]) ) + (cos(t[N2-2])-cos(t.back())) * ( f[0    ].back()[0]*(p[1]-p[0]) + f[0    ].back().back()*(p.back()-p[N3-2]) ) )
	                + (pow(r.back(),3)-pow(r[N1-2],3)) * ( (cos(t[0])-cos(t[1])) * ( f.back()[0][0]*(p[1]-p[0]) + f.back()[0].back()*(p.back()-p[N3-2]) ) + (cos(t[N2-2])-cos(t.back())) * ( f.back().back()[0]*(p[1]-p[0]) + f.back().back().back()*(p.back()-p[N3-2]) ) ) );
	
	return ret;
}