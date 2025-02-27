#include <complex>
#include <vector>

//----- First derivative -----
double dx( double (*f)(double), double x, double h=1e-10 ){
	// First derivative of a univariate function f(x) at the point x, given by the traditional definition of the derivative.
	return ( f( x + h ) - f( x ) ) / h;
}
double dx( double (*f)(double,double), double x, double y, double h=1e-10 ){
	// First partial x-derivative of a bivariate function f(x,y) at the point (x,y), given by the traditional definition of the derivative.
	return ( f( x + h, y ) - f( x, y ) ) / h;
}
double dx( double (*f)(double,int), double x, int n, double h=1e-10 ){
	// First derivative of a univariate function f(x) with optional integer paramter n, at the point x, given by the traditional definition of the derivative.
	return ( f( x + h, n ) - f( x, n ) ) / h;
}
double dx( double (*f)(double,int,int), double x, int n, int m, double h=1e-10 ){
	// First derivative of a univariate function f(x) with two optional integer paramters n,m, at the point x, given by the traditional definition of the derivative.
	return ( f( x + h, n, m ) - f( x, n, m ) ) / h;
}
std::complex<double> dx( std::complex<double> (*f)(double), double x, double h=1e-10 ){
	// First derivative of a univariate complex function f(x) at the point x, given by the traditional definition of the derivative.
	return ( f( x + h ) - f( x ) ) / h;
}
std::complex<double> dx( std::complex<double> (*f)(double,double), double x, double n, double h=1e-10 ){
	// First derivative of a univariate complex function f(x) with extra parameter n, at the point x, given by the traditional definition of the derivative.
	return ( f( x + h, n ) - f( x, n ) ) / h;
}
std::complex<double> dx( std::complex<double> (*f)(double,double,int,int), double x, double y, int n, int m, double h=1e-10 ){
	// First partial x-derivative of a bivariate function f(x,y) with two optional integer paramters n,m, at the point (x,y), given by the traditional definition of the derivative.
	return ( f( x + h, y, n, m ) - f( x, y, n, m ) ) / h;
}


double dx_symmetric( double (*f)(double), double x, double h=1e-10 ){
	// First derivative of a univariate function f(x) at the point x, given by the definition of the symmetric derivative.
	return ( f( x + h ) - f( x - h ) ) / ( 2 * h );
}
double dx_symmetric( double (*f)(double,double), double x, double y, double h=1e-10 ){
	// First partial x-derivative of a bivariate function f(x,y) at the point (x,y), given by the definition of the symmetric derivative.
	return ( f( x + h, y ) - f( x - h, y ) ) / ( 2 * h );
}
double dx_symmetric( double (*f)(double,int), double x, int n, double h=1e-10 ){
	// First derivative of a univariate function f(x) with optional integer parameter n, at the point x, given by the definition of the symmetric derivative.
	return ( f( x + h, n ) - f( x - h, n ) ) / ( 2 * h );
}
double dx_symmetric( double (*f)(double,int,int), double x, int n, int m, double h=1e-10 ){
	// First derivative of a univariate function f(x) with two optional integer parameters n,m, at the point x, given by the definition of the symmetric derivative.
	return ( f( x + h, n, m ) - f( x - h, n, m ) ) / ( 2 * h );
}
std::complex<double> dx_symmetric( std::complex<double> (*f)(double), double x, double h=1e-10 ){
	// First derivative of a complex univariate function f(x) at the point x, given by the definition of the symmetric derivative.
	return ( f( x + h ) - f( x - h ) ) / ( 2 * h );
}
std::complex<double> dx_symmetric( std::complex<double> (*f)(double,double), double x, double n, double h=1e-10 ){
	// First derivative of a complex univariate function f(x) with extra parameter n, at the point x, given by the definition of the symmetric derivative.
	return ( f( x + h, n ) - f( x - h, n ) ) / ( 2 * h );
}
std::complex<double> dx_symmetric( std::complex<double> (*f)(double,double,int,int), double x, double y, int n, int m, double h=1e-10 ){
	// First partial x-derivative of a bivariate function f(x,y) with two optional integer parameters n,m, at the point (x,y), given by the definition of the symmetric derivative.
	return ( f( x + h, y, n, m ) - f( x - h, y, n, m ) ) / ( 2 * h );
}


double partial_x1( double (*f)(double,double), double x1, double x2, double h=1e-10 ){
	// First partial x1-derivative of a bivariate function f(x1,x2) at the point (x1,x2) given by the traditional definition of the derivative.
	return ( f( x1 + h, x2     ) - f( x1, x2 ) ) / h;
}
double partial_x2( double (*f)(double,double), double x1, double x2, double h=1e-10 ){
	// First partial x2-derivative of a bivariate function f(x1,x2) at the point (x1,x2) given by the traditional definition of the derivative.
	return ( f( x1    , x2 + h ) - f( x1, x2 ) ) / h;
}

double partial_x1_symmetric( double (*f)(double,double), double x1, double x2, double h=1e-10 ){
	// First partial x1-derivative of a bivariate function f(x1,x2) at the point (x1,x2), given by the definition of the symmetric derivative.
	return ( f( x1 + h , x2     ) - f( x1 - h, x2    ) ) / ( 2 * h );
}
double partial_x2_symmetric( double (*f)(double,double), double x1, double x2, double h=1e-10 ){
	// First partial x1-derivative of a bivariate function f(x1,x2) at the point (x1,x2), given by the definition of the symmetric derivative.
	return ( f( x1     , x2 + h ) - f( x1    , x2 - h ) ) / ( 2 * h );
}



double partial_x1( double (*f)(double,double,double), double x1, double x2, double x3, double h=1e-10 ){
	// First partial x1-derivative of a trivariate function f(x1,x2,x3) at the point (x1,x2,x3) given by the traditional definition of the derivative.
	return ( f( x1 + h, x2    , x3     ) - f( x1, x2, x3 ) ) / h;
}
double partial_x2( double (*f)(double,double,double), double x1, double x2, double x3, double h=1e-10 ){
	// First partial x2-derivative of a trivariate function f(x,y,z) at the point (x,y,z) given by the traditional definition of the derivative.
	return ( f( x1    , x2 + h, x3     ) - f( x1, x2, x3 ) ) / h;
}
double partial_x3( double (*f)(double,double,double), double x1, double x2, double x3, double h=1e-10 ){
	// First partial x3-derivative of a trivariate function f(x,y,z) at the point (x,y,z) given by the traditional definition of the derivative.
	return ( f( x1    , x2    , x3 + h ) - f( x1, x2, x3 ) ) / h;
}

double partial_x1_symmetric( double (*f)(double,double,double), double x1, double x2, double x3, double h=1e-10 ){
	// First partial x1-derivative of a trivariate function f(x1,x2,x3) at the point (x1,x2,x3), given by the definition of the symmetric derivative.
	return ( f( x1 + h, x2    , x3     ) - f( x1 - h, x2    , x3     ) ) / ( 2 * h );
}
double partial_x2_symmetric( double (*f)(double,double,double), double x1, double x2, double x3, double h=1e-10 ){
	// First partial x2-derivative of a trivariate function f(x1,x2,x3) at the point (x1,x2,x3), given by the definition of the symmetric derivative.
	return ( f( x1    , x2 + h, x3     ) - f( x1    , x2 - h, x3     ) ) / ( 2 * h );
}
double partial_x3_symmetric( double (*f)(double,double,double), double x1, double x2, double x3, double h=1e-10 ){
	// First partial x3-derivative of a trivariate function f(x1,x2,x3) at the point (x1,x2,x3), given by the definition of the symmetric derivative.
	return ( f( x1    , x2    , x3 + h ) - f( x1    , x2    , x3 - h ) ) / ( 2 * h );
}


//----- Second derivatives (note that the factor 1/h^2 means h can only be around sqrt() the size used for a first derivative, before numerical error dominates) -----
double dx2( double (*f)(double), double x, double h=1e-5 ){
	// Seond derivative of a univariate function f(x) at the point x, given by the traditional definition of the derivative.
	return ( f( x + 2 * h ) - 2 * f( x + h ) + f( x ) ) / ( h * h );
}
double dx2( double (*f)(double,double), double x, double y, double h=1e-5 ){
	// Seond partial x-derivative of a bivariate function f(x,y) at the point (x,y), given by the traditional definition of the derivative.
	return ( f( x + 2 * h, y ) - 2 * f( x + h, y ) + f( x, y ) ) / ( h * h );
}
double dx2( double (*f)(double,int), double x, int n, double h=1e-5 ){
	// Seond derivative of a univariate function f(x) with optional integer parameter n, at the point x, given by the traditional definition of the derivative.
	return ( f( x + 2 * h, n ) - 2 * f( x + h, n ) + f( x, n ) ) / ( h * h );
}
double dx2( double (*f)(double,int,int), double x, int n, int m, double h=1e-5 ){
	// Seond derivative of a univariate function f(x) with two optional integer parameters n,m, at the point x, given by the traditional definition of the derivative.
	return ( f( x + 2 * h, n, m ) - 2 * f( x + h, n, m ) + f( x, n, m ) ) / ( h * h );
}
std::complex<double> dx2( std::complex<double> (*f)(double), double x, double h=1e-5 ){
	// Seond derivative of a univariate complex function f(x) at the point x, given by the traditional definition of the derivative.
	return ( f( x + 2 * h ) - 2.0 * f( x + h ) + f( x ) ) / ( h * h );
}
std::complex<double> dx2( std::complex<double> (*f)(double,double), double x, double n, double h=1e-5 ){
	// Seond derivative of a univariate complex function f(x) with extra parameter n, at the point x, given by the traditional definition of the derivative.
	return ( f( x + 2 * h, n ) - 2.0 * f( x + h, n ) + f( x, n ) ) / ( h * h );
}
std::complex<double> dx2( std::complex<double> (*f)(double,double,int,int), double x, double y, int n, int m, double h=1e-5 ){
	// Seond partial x-derivative of a bivariate function f(x,y) with two optional integer parameters n,m, at the point (x,y), given by the traditional definition of the derivative.
	return ( f( x + 2 * h, y, n, m ) - 2.0 * f( x + h, y, n, m ) + f( x, y, n, m ) ) / ( h * h );
}


double dx2_symmetric( double (*f)(double), double x, double h=1e-5 ){
	// Second derivative of a univariate function f(x) at the point x, given by the definition of the symmetric derivative.
	return ( f( x + h ) - 2 * f( x ) + f( x - h ) ) / ( h * h );
}
double dx2_symmetric( double (*f)(double,double), double x, double y, double h=1e-5 ){
	// Second partial x-derivative of a bivariate function f(x,y) at the point (x,y), given by the definition of the symmetric derivative.
	return ( f( x + h, y ) - 2 * f( x, y ) + f( x - h, y ) ) / ( h * h );
}
double dx2_symmetric( double (*f)(double,int), double x, int n, double h=1e-5 ){
	// Second derivative of a univariate function f(x) with optional integer parameter n, at the point x, given by the definition of the symmetric derivative.
	return ( f( x + h, n ) - 2 * f( x, n ) + f( x - h, n ) ) / ( h * h );
}
double dx2_symmetric( double (*f)(double,int,int), double x, int n, int m, double h=1e-5 ){
	// Second derivative of a univariate function f(x) with two optional integer parameters n,m at the point x, given by the definition of the symmetric derivative.
	return ( f( x + h, n, m ) - 2 * f( x, n, m ) + f( x - h, n, m ) ) / ( h * h );
}
std::complex<double> dx2_symmetric( std::complex<double> (*f)(double), double x, double h=1e-5 ){
	// Second derivative of a univariate complex function f(x) at the point x, given by the definition of the symmetric derivative.
	return ( f( x + h ) - 2.0 * f( x ) + f( x - h ) ) / ( h * h );
}
std::complex<double> dx2_symmetric( std::complex<double> (*f)(double,double), double x, double n, double h=1e-5 ){
	// Second derivative of a univariate complex function f(x) with extra parameter n, at the point x, given by the definition of the symmetric derivative.
	return ( f( x + h, n ) - 2.0 * f( x, n ) + f( x - h, n ) ) / ( h * h );
}
std::complex<double> dx2_symmetric( std::complex<double> (*f)(double,double,int,int), double x, double y, int n, int m, double h=1e-5 ){
	// Second partial x-derivative of a bivariate function f(x,y) with two optional integer parameters n,m, at the point (x,y), given by the definition of the symmetric derivative.
	return ( f( x + h, y, n, m ) - 2.0 * f( x, y, n, m ) + f( x - h, y, n, m ) ) / ( h * h );
}

double partial2_x1( double (*f)(double,double,double), double x1, double x2, double x3, double h=1e-5 ){
	// Second partial x1-derivative of a trivariate function f(x1,x2,x3) at the point (x1,x2,x3) given by the traditional definition of the derivative.
	return ( f( x1 + 2*h, x2      , x3       ) - 2 * f( x1 + h, x2    , x3     ) + f( x1, x2, x3 ) ) / ( h * h );
}
double partial2_x2( double (*f)(double,double,double), double x1, double x2, double x3, double h=1e-5 ){
	// Second partial x2-derivative of a trivariate function f(x1,x2,x3) at the point (x1,x2,x3) given by the traditional definition of the derivative.
	return ( f( x1      , x2 + 2*h, x3       ) - 2 * f( x1    , x2 + h, x3     ) + f( x1, x2, x3 ) ) / ( h * h );
}
double partial2_x3( double (*f)(double,double,double), double x1, double x2, double x3, double h=1e-5 ){
	// Second partial x3-derivative of a trivariate function f(x1,x2,x3) at the point (x1,x2,x3) given by the traditional definition of the derivative.
	return ( f( x1      , x2      , x3 + 2*h ) - 2 * f( x1    , x2    , x3 + h ) + f( x1, x2, x3 ) ) / ( h * h );
}

double partial2_x1_symmetric( double (*f)(double,double,double), double x1, double x2, double x3, double h=1e-5 ){
	// Second partial x1-derivative of a trivariate function f(x1,x2,x3) at the point (x1,x2,x3), given by the definition of the symmetric derivative.
	return ( f( x1 + h, x2    , x3     ) - 2 * f( x1, x2, x3 ) + f( x1 - h, x2    , x3     ) ) / ( h * h );
}
double partial2_x2_symmetric( double (*f)(double,double,double), double x1, double x2, double x3, double h=1e-5 ){
	// Second partial x1-derivative of a trivariate function f(x1,x2,x3) at the point (x1,x2,x3), given by the definition of the symmetric derivative.
	return ( f( x1    , x2 + h, x3     ) - 2 * f( x1, x2, x3 ) + f( x1    , x2 - h, x3     ) ) / ( h * h );
}
double partial2_x3_symmetric( double (*f)(double,double,double), double x1, double x2, double x3, double h=1e-5 ){
	// Second partial x1-derivative of a trivariate function f(x1,x2,x3) at the point (x1,x2,x3), given by the definition of the symmetric derivative.
	return ( f( x1    , x2    , x3 + h ) - 2 * f( x1, x2, x3 ) + f( x1    , x2    , x3 - h ) ) / ( h * h );
}


//----- Mixed second partial derivatives -----
double partial2_x2_x1( double (*f)(double,double,double), double x1, double x2, double x3, double h=1e-5 ){
	// Mixed second partial x2-x1-derivative of a trivariate function f(x1,x2,x3) at the point (x1,x2,x3) given by the traditional definition of the derivative.
	return ( f( x1 + h, x2 + h, x3     ) - f( x1 + h, x2    , x3     ) - f( x1    , x2 + h, x3     ) + f( x1, x2, x3 ) ) / ( h * h );
}
double partial2_x3_x1( double (*f)(double,double,double), double x1, double x2, double x3, double h=1e-5 ){
	// Mixed second partial x3-x1-derivative of a trivariate function f(x1,x2,x3) at the point (x1,x2,x3) given by the traditional definition of the derivative.
	return ( f( x1 + h, x2    , x3 + h ) - f( x1 + h, x2    , x3     ) - f( x1    , x2    , x3 + h ) + f( x1, x2, x3 ) ) / ( h * h );
}
double partial2_x3_x2( double (*f)(double,double,double), double x1, double x2, double x3, double h=1e-5 ){
	// Mixed second partial x3-x2-derivative of a trivariate function f(x1,x2,x3) at the point (x1,x2,x3) given by the traditional definition of the derivative.
	return ( f( x1    , x2 + h, x3 + h ) - f( x1    , x2 + h, x3     ) - f( x1    , x2    , x3 + h ) + f( x1, x2, x3 ) ) / ( h * h );
}
double partial2_x1_x2( double (*f)(double,double,double), double x1, double x2, double x3, double h=1e-5 ){ return partial2_x2_x1( f, x1, x2, x3, h ); }
double partial2_x1_x3( double (*f)(double,double,double), double x1, double x2, double x3, double h=1e-5 ){ return partial2_x3_x1( f, x1, x2, x3, h ); }
double partial2_x2_x3( double (*f)(double,double,double), double x1, double x2, double x3, double h=1e-5 ){ return partial2_x3_x2( f, x1, x2, x3, h ); }

double partial2_x2_x1_symmetric( double (*f)(double,double,double), double x1, double x2, double x3, double h=1e-5 ){
	// Mixed second partial x2-x1-derivative of a trivariate function f(x1,x2,x3) at the point (x1,x2,x3) given by the symmetric definition of the derivative.
	return ( f( x1 + h, x2 + h, x3     ) - f( x1 - h, x2 + h, x3     ) - f( x1 + h, x2 - h, x3     ) + f( x1 - h, x2 - h, x3     ) ) / ( 4 * h * h );
}
double partial2_x3_x1_symmetric( double (*f)(double,double,double), double x1, double x2, double x3, double h=1e-5 ){
	// Mixed second partial x3-x1-derivative of a trivariate function f(x1,x2,x3) at the point (x1,x2,x3) given by the symmetric definition of the derivative.
	return ( f( x1 + h, x2    , x3 + h ) - f( x1 - h, x2    , x3 + h ) - f( x1 + h, x2    , x3 - h ) + f( x1 - h, x2    , x3 - h ) ) / ( 4 * h * h );
}
double partial2_x3_x2_symmetric( double (*f)(double,double,double), double x1, double x2, double x3, double h=1e-5 ){
	// Mixed second partial x3-x2-derivative of a trivariate function f(x1,x2,x3) at the point (x1,x2,x3) given by the symmetric definition of the derivative.
	return ( f( x1    , x2 + h, x3 + h ) - f( x1    , x2 - h, x3 + h ) - f( x1    , x2 + h, x3 - h ) + f( x1    , x2 - h, x3 - h ) ) / ( 4 * h * h );
}
double partial2_x1_x2_symmetric( double (*f)(double,double,double), double x1, double x2, double x3, double h=1e-5 ){ return partial2_x2_x1_symmetric( f, x1, x2, x3, h ); }
double partial2_x1_x3_symmetric( double (*f)(double,double,double), double x1, double x2, double x3, double h=1e-5 ){ return partial2_x3_x1_symmetric( f, x1, x2, x3, h ); }
double partial2_x2_x3_symmetric( double (*f)(double,double,double), double x1, double x2, double x3, double h=1e-5 ){ return partial2_x3_x2_symmetric( f, x1, x2, x3, h ); }




//----- Divergence of a vector -----
double div_cartesian( double (*A_x)(double,double,double), double (*A_y)(double,double,double), double (*A_z)(double,double,double), double x, double y, double z ){
	return partial_x1(A_x,x,y,z) + partial_x2(A_y,x,y,z) + partial_x3(A_z,x,y,z);
}

double div_cartesian_cout( double (*A_x)(double,double,double), double (*A_y)(double,double,double), double (*A_z)(double,double,double), double x, double y, double z ){
	std::cout << "d/dx:\t" << partial_x1(A_x,x,y,z) << std::endl;
	std::cout << "d/dy:\t" << partial_x2(A_y,x,y,z) << std::endl;
	std::cout << "d/dz:\t" << partial_x3(A_z,x,y,z) << std::endl;
	return partial_x1(A_x,x,y,z) + partial_x2(A_y,x,y,z) + partial_x3(A_z,x,y,z);
}

double div_cylindrical( double (*A_varpi)(double,double,double), double (*A_phi)(double,double,double), double (*A_z)(double,double,double), double varpi, double phi, double z ){
	return partial_x1(A_varpi,varpi,phi,z) + pow(varpi,-1) * ( A_varpi(varpi,phi,z) + partial_x2(A_phi,varpi,phi,z) ) + partial_x3(A_z,varpi,phi,z);
}

double div_spherical( double (*A_r)(double,double,double), double (*A_t)(double,double,double), double (*A_p)(double,double,double), double r, double t, double p ){
	// Divergence of a vector in spherical polar coordinates, whose components are each given by separate functions named above.
	//--- V1 ---
	/*double ret = 0;
	ret += partial_x1(A_r,r,t,p);
	ret += 2.0 * pow(r,-1) * A_r(r,t,p);
	ret += pow(r,-1) * partial_x2(A_t,r,t,p);
	ret += cos(t) * pow( r*sin(t), -1 ) * A_t(r,t,p);
	ret += pow( r*sin(t), -1 ) * partial_x3(A_p,r,t,p);
	return ret;*/
	
	//--- V2 ---
	return partial_x1(A_r,r,t,p) + pow(r,-1) * ( 2.0*A_r(r,t,p) + partial_x2(A_t,r,t,p) + pow(sin(t),-1) * ( cos(t)*A_t(r,t,p) + partial_x3(A_p,r,t,p) ) );
	
}

double div_spherical( double (*A_r)(double,double), double (*A_t)(double,double), double (*A_p)(double,double), double r, double t ){
	// Divergence of an axisymmetric (phi-independent) vector in spherical polar coordinates, whose components are each given by separate functions named above.
	return partial_x1(A_r,r,t) + pow(r,-1) * ( 2.0*A_r(r,t) + partial_x2(A_t,r,t) + cos(t)*pow(sin(t),-1)*A_t(r,t) );
	
}




//----- Curl of a vector -----
std::vector<double> curl_cartesian( double (*A_x)(double,double,double), double (*A_y)(double,double,double), double (*A_z)(double,double,double), double x, double y, double z ){
	std::vector<double> ret ( 3 );
	ret[0] = partial_x2(A_z,x,y,z) - partial_x3(A_y,x,y,z);
	ret[1] = partial_x3(A_x,x,y,z) - partial_x1(A_z,x,y,z);
	ret[2] = partial_x1(A_y,x,y,z) - partial_x2(A_x,x,y,z);
	return ret;
}

std::vector<double> curl_cylindrical( double (*A_varpi)(double,double,double), double (*A_phi)(double,double,double), double (*A_z)(double,double,double), double varpi, double phi, double z ){
	std::vector<double> ret ( 3 );
	ret[0] = pow(varpi,-1) * partial_x2(A_z,varpi,phi,z) - partial_x3(A_phi,varpi,phi,z);
	ret[1] = partial_x3(A_varpi,varpi,phi,z) - partial_x1(A_z,varpi,phi,z);
	ret[2] = partial_x1(A_phi,varpi,phi,z) + pow(varpi,-1) * ( A_phi(varpi,phi,z) - partial_x2(A_varpi,varpi,phi,z) );
	return ret;
}

std::vector<double> curl_spherical( double (*A_r)(double,double,double), double (*A_t)(double,double,double), double (*A_p)(double,double,double), double r, double t, double p ){
	// Curl of a vector in spherical polar coordinates, whose components are each given by separate functions named above.
	// In most cases, #include<vector> will already be used within the parent file, so it shouldn't be needed here. But be aware that it could be the cause of any issues.
	std::vector<double> ret ( 3 );
	ret[0] = pow(r,-1) * ( partial_x2(A_p,r,t,p) + pow(sin(t),-1) * ( cos(t)*A_p(r,t,p) - partial_x3(A_t,r,t,p) ) );
	ret[1] = pow(r*sin(t),-1) * partial_x3(A_r,r,t,p) - partial_x1(A_p,r,t,p) - pow(r,-1)*A_p(r,t,p);
	ret[2] = partial_x1(A_t,r,t,p) + pow(r,-1) * ( A_t(r,t,p) - partial_x2(A_r,r,t,p) );
	return ret;
}

std::vector<double> curl_spherical( double (*A_r)(double,double), double (*A_t)(double,double), double (*A_p)(double,double), double r, double t ){
	// Curl of an axisymmetric (phi-independent) vector in spherical polar coordinates, whose components are each given by separate functions named above.
	std::vector<double> ret ( 3 );
	ret[0] = pow(r,-1) * ( partial_x2(A_p,r,t) + pow(sin(t),-1) * ( cos(t)*A_p(r,t) ) );
	ret[1] = - partial_x1(A_p,r,t) - pow(r,-1)*A_p(r,t);
	ret[2] = partial_x1(A_t,r,t) + pow(r,-1) * ( A_t(r,t) - partial_x2(A_r,r,t) );
	return ret;
}

std::vector<double> curl_spherical_cout( double (*A_r)(double,double), double (*A_t)(double,double), double (*A_p)(double,double), double r, double t ){
	// Curl of an axisymmetric (phi-independent) vector in spherical polar coordinates, whose components are each given by separate functions named above.
	// Cout each of the summed values so that mathematical errors can be found.
	std::vector<double> ret ( 3 );
	double r1 = pow(r,-1) * partial_x2(A_p,r,t);
	double r2 = pow(r*sin(t),-1) * ( cos(t)*A_p(r,t) );
	double t1 = - partial_x1(A_p,r,t);
	double t2 = - pow(r,-1) * A_p(r,t);
	double p1 = partial_x1(A_t,r,t);
	double p2 = pow(r,-1) * A_t(r,t);
	double p3 = pow(r,-1) * - partial_x2(A_r,r,t);
	ret[0] = r1 + r2;
	ret[1] = t1 + t2;
	ret[2] = p1 + p2 + p3;
	std::cout << "r:\t" << r1 <<"\t"<< r2 << std::endl;
	std::cout << "t:\t" << t1 <<"\t"<< t2 << std::endl;
	std::cout << "p:\t" << p1 <<"\t"<< p2 <<"\t"<< p3 << std::endl;
	return ret;
}