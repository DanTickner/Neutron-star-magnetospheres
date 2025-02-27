#include <vector>
#include <math.h>
#include <complex>
#include <string>
#include <iomanip>

std::vector<double> cart_to_spher( std::vector<double> v ){
	// Converts a vector with three Cartesian components to one with three spherical polar components.
	double vr = pow( v[0]*v[0] + v[1]*v[1] + v[2]*v[2], 0.5 );
	return std::vector<double> { vr, acos( v[2] / vr ), atan2( v[1], v[0] ) };
}




std::vector<double> spher_to_cart( std::vector<double> v ){
	// Converts a vector with three spherical polar components to one with three Cartesian components.
	return std::vector<double> { v[0]*sin(v[1])*cos(v[2]), v[0]*sin(v[1])*sin(v[2]), v[0]*cos(v[1]) };
}



void cout_vector( std::vector<int> v ){
	// Prints an n-dimensional vector to the screen.
	for( int i=0; i<v.size(); i++ ){
		std::cout << v[i] << "\t";
	}
	std::cout << std::endl;
}
void cout_vector( std::vector<double> v ){
	// Prints an n-dimensional vector to the screen.
	for( int i=0; i<v.size(); i++ ){
		std::cout << v[i] << "\t";
	}
	std::cout << std::endl;
}
void cout_vector( std::vector< std::complex<double> > v ){
	// Prints an n-dimensional vector to the screen.
	for( int i=0; i<v.size(); i++ ){
		std::cout << v[i] << "\t";
	}
	std::cout << std::endl;
}
void cout_vector( std::vector<std::string> v ){
	// Prints an n-dimensional vector to the screen.
	for( int i=0; i<v.size(); i++ ){
		std::cout << v[i] << "\t";
	}
	std::cout << std::endl;
}
/*void cout_vector( std::vector<double> v ){
	// Prints a three-element vector to the screen.
	std::cout << v[0] << "\t" << v[1] << "\t" << v[2] << std::endl;
}*/
/*void cout_vector( std::vector<std::complex<double>> v ){
	// Prints a three-element vector with complex elements to the screen.
	std::cout << v[0] << "\t" << v[1] << "\t" << v[2] << std::endl;
}*/
void cout_vector( std::vector<double> v, std::string preamble ){
	// Prints an n-dimensional vector to the screen.
	std::cout << preamble << ":\t";
	for( int i=0; i<v.size(); i++ ){
		std::cout << v[i] << "\t";
	}
	std::cout << std::endl;
}
void cout_vector( std::vector< std::complex<double> > v, std::string preamble ){
	// Prints an n-dimensional vector to the screen.
	std::cout << preamble << ":\t";
	for( int i=0; i<v.size(); i++ ){
		std::cout << v[i] << "\t";
	}
	std::cout << std::endl;
}
void cout_vector( std::vector<double> v, std::string preamble, std::string postscript ){
	// Prints an n-dimensional vector to the screen.
	std::cout << preamble << ":\t";
	for( int i=0; i<v.size(); i++ ){
		std::cout << v[i] << "\t";
	}
	std::cout << postscript << std::endl;
}
void cout_vector( std::vector< std::complex<double> > v, std::string preamble, std::string postscript ){
	// Prints an n-dimensional vector to the screen.
	std::cout << preamble << ":\t";
	for( int i=0; i<v.size(); i++ ){
		std::cout << v[i] << "\t";
	}
	std::cout << postscript << std::endl;
}


void cout_vector_setw( std::vector<double> v, int w, std::string preamble ){
	// Prints an n-dimensional vector to the screen.
	std::cout << preamble << ":\t";
	for( int i=0; i<v.size(); i++ ){
		std::cout <<std::left<<std::setw(w)<< v[i];
	}
	std::cout << std::endl;
}
void cout_vector_setw( std::vector< std::complex<double> > v, int w, std::string preamble ){
	// Prints an n-dimensional vector to the screen.
	std::cout << preamble << ":\t";
	for( int i=0; i<v.size(); i++ ){
		std::cout <<std::left<<std::setw(w)<< v[i];
	}
	std::cout << std::endl;
}
void cout_vector_setw( std::vector< std::complex<double> > v, int w, std::string preamble, std::string postscript ){
	// Prints an n-dimensional vector to the screen.
	std::cout << preamble << ":\t";
	for( int i=0; i<v.size(); i++ ){
		std::cout <<std::left<<std::setw(w)<< v[i];
	}
	std::cout << postscript << std::endl;
}




std::vector<double> vector_sum_cart( std::vector<double> a, std::vector<double> b ){
	// Sum of two vectors which are given in Cartesian coordinates. Result is in Cartesian coordinates.
	return std::vector<double> { a[0]+b[0], a[1]+b[1], a[2]+b[2] };
}
std::vector< std::complex<double> > vector_sum_cart( std::vector< std::complex<double> > a, std::vector< std::complex<double> > b ){
	// Sum of two vectors which are given in Cartesian coordinates. Result is in Cartesian coordinates.
	return std::vector< std::complex<double> > { a[0]+b[0], a[1]+b[1], a[2]+b[2] };
}
std::vector<double> vector_sum( std::vector<double> a, std::vector<double> b ){
	// Sum of two vectors which are given in Cartesian coordinates. Result is in Cartesian coordinates.
	return std::vector<double> { a[0]+b[0], a[1]+b[1], a[2]+b[2] };
}
std::vector< std::complex<double> > vector_sum( std::vector< std::complex<double> > a, std::vector< std::complex<double> > b ){
	// Sum of two vectors which are given in Cartesian coordinates. Result is in Cartesian coordinates.
	return std::vector< std::complex<double> > { a[0]+b[0], a[1]+b[1], a[2]+b[2] };
}
std::vector< std::complex<double> > vector_sum( std::vector<double> a, std::vector< std::complex<double> > b ){
	// Sum of two vectors which are given in Cartesian coordinates. Result is in Cartesian coordinates.
	return std::vector< std::complex<double> > { a[0]+b[0], a[1]+b[1], a[2]+b[2] };
}
std::vector< std::complex<double> > vector_sum( std::vector< std::complex<double> > a, std::vector<double> b ){
	// Sum of two vectors which are given in Cartesian coordinates. Result is in Cartesian coordinates.
	return std::vector< std::complex<double> > { a[0]+b[0], a[1]+b[1], a[2]+b[2] };
}




std::vector<double> vector_diff_cart( std::vector<double> a, std::vector<double> b ){
	// Difference of two vectors which are given in Cartesian coordinates. Result is in Cartesian coordinates.
	return std::vector<double> { a[0]-b[0], a[1]-b[1], a[2]-b[2] };
}




std::vector<double> vector_sum_spher( std::vector<double> a, std::vector<double> b ){
	// Sum of two vectors which are given in spherical coordinates. Result is in spherical coordinates.
	double r     = pow( a[0]*a[0] + b[0]*b[0] + 2*a[0]*b[0]*( sin(a[1]) * sin(b[1]) * cos(a[2]-b[2]) + cos(a[1]) * cos(b[1]) ), 0.5 );
	double theta = acos( ( a[0] * cos(a[1]) + b[0] * cos(b[1]) ) / r );
	double phi   = atan2( a[0] * sin(a[1]) * sin(a[2]) + b[0] * sin(b[1]) * sin(b[2]), a[0] * sin(a[1]) * cos(a[2]) + b[0] * sin(b[1]) * cos(b[2]) );
	return std::vector<double> { r, theta, phi };
}
std::vector< std::complex<double> > vector_sum_spher( std::vector< std::complex<double> > a, std::vector< std::complex<double> > b ){
	// Sum of two vectors which are given in spherical coordinates. Result is in spherical coordinates.
	// NOT TESTED MATHEMATICALLY.
	std::complex<double> r     = pow( a[0]*a[0] + b[0]*b[0] + 2.0*a[0]*b[0]*( sin(a[1]) * sin(b[1]) * cos(a[2]-b[2]) + cos(a[1]) * cos(b[1]) ), 0.5 );
	std::complex<double> theta = acos( ( a[0] * cos(a[1]) + b[0] * cos(b[1]) ) / r );
	
	std::complex<double> phi_y = a[0] * sin(a[1]) * sin(a[2]) + b[0] * sin(b[1]) * sin(b[2]);
	std::complex<double> phi_x = a[0] * sin(a[1]) * cos(a[2]) + b[0] * sin(b[1]) * cos(b[2]);
	//std::cout << phi_y <<"\t"<< phi_x << std::endl;
	//std::cout << a[0] * sin(a[1]) * sin(a[2]) <<"\t"<< sin(a[2]) + b[0] * sin(b[1]) * sin(b[2]) <<"\t\t"<< a[0] * sin(a[1]) * cos(a[2]) <<"\t"<< cos(a[2]) + b[0] * sin(b[1]) * cos(b[2]) << std::endl;
	std::complex<double> phi   = 0;
	if( not( ( abs(phi_y)<1e-12 ) && ( abs(phi_x)<1e-12 ) ) ){
		//std::cout << "updating phi" << std::endl;
		phi = atan( phi_y / phi_x );
	}
	return std::vector< std::complex<double> > { r, theta, phi };
}




std::vector<double> vector_diff_spher( std::vector<double> a, std::vector<double> b ){
	// Difference of two vectors which are given in spherical coordinates. Result is in spherical coordinates.
	double r = pow( a[0]*a[0] + b[0]*b[0] - 2*a[0]*b[0]*( sin(a[1]) * sin(b[1]) * cos(a[2]-b[2]) + cos(a[1]) * cos(b[1]) ), 0.5 );
	double theta = acos( ( a[0] * cos(a[1]) - b[0] * cos(b[1]) ) / r );
	double phi = atan2( a[0] * sin(a[1]) * sin(a[2]) - b[0] * sin(b[1]) * sin(b[2]), a[0] * sin(a[1]) * cos(a[2]) - b[0] * sin(b[1]) * cos(b[2]) );
	return std::vector<double> { r, theta, phi };
}




std::vector<double> scalar_times_vector_cart( double A, std::vector<double> v ){
	// Scalar multiplication of a vector in Cartesian coordinates. Result is in cartesian coordinates.
	return std::vector<double> { A*v[0], A*v[1], A*v[2] };
}
std::vector< std::complex<double> > scalar_times_vector_cart( std::complex<double> A, std::vector< std::complex<double> > v ){
	// Scalar multiplication of a vector in Cartesian coordinates. Result is in cartesian coordinates.
	return std::vector< std::complex<double> > { A*v[0], A*v[1], A*v[2] };
}

std::vector<double> scalar_times_vector( double A, std::vector<double> v ){
	// General function which doesn't specify the coordinate system. Useful for, e.g., Cartesian, or for spherical polar vectors both defined at the same point.
	return std::vector<double> { A*v[0], A*v[1], A*v[2] };
}
std::vector< std::complex<double> > scalar_times_vector( std::complex<double> A, std::vector< std::complex<double> > v ){
	// General function which doesn't specify the coordinate system. Useful for, e.g., Cartesian, or for spherical polar vectors both defined at the same point.
	return std::vector< std::complex<double> > { A*v[0], A*v[1], A*v[2] };
}



std::vector<double> scalar_times_vector_spher( double A, std::vector<double> v ){
	// Scalar multiplication of a vector in spherical coordinates. Result is in spherical coordinates.
	if( A > 0 ){
		return std::vector<double> { A*v[0], v[1], v[2] };
	} else if( A < 0 ){
		return std::vector<double> { abs(A)*v[0], 3.14159265358979323846 - v[1], atan2( -sin(v[1])*sin(v[2]), -sin(v[1])*cos(v[2]) ) };
	}
	else {
		return std::vector<double> { 0, 0, 0 };
	}
}
std::vector< std::complex<double> > scalar_times_vector_spher( std::complex<double> A, std::vector< std::complex<double> > v ){
	// Scalar multiplication of a vector in spherical coordinates. Result is in spherical coordinates.
	// NOT PROPERLY CHECKED MATHEMATICALLY. WHAT HAPPENS IF re(A)<0?
	return std::vector< std::complex<double> > { A*v[0], v[1], v[2] };
}




std::vector<double> cross_product_cart( std::vector<double> a, std::vector<double> b ){
	return std::vector<double> { a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0] };
}

std::vector<double> cross_product( std::vector<double> a, std::vector<double> b ){
	// General function which doesn't specify the coordinate system. Useful for, e.g., Cartesian, or for spherical polar vectors both defined at the same point.
	return std::vector<double> { a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0] };
}
std::vector< std::complex<double> > cross_product( std::vector< std::complex<double> > a, std::vector< std::complex<double> > b ){
	// General function which doesn't specify the coordinate system. Useful for, e.g., Cartesian, or for spherical polar vectors both defined at the same point.
	return std::vector< std::complex<double> > { a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0] };
}

double dot_product_cart( std::vector<double> a, std::vector<double> b ){
	return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}
std::complex<double> dot_product_cart( std::vector<double> a, std::vector< std::complex<double> > b ){
	return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}
std::complex<double> dot_product_cart( std::vector< std::complex<double> > a, std::vector<double> b ){
	return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}
std::complex<double> dot_product_cart( std::vector< std::complex<double> > a, std::vector< std::complex<double> > b ){
	return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

double dot_product( std::vector<double> a, std::vector<double> b ){
	return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}
std::complex<double> dot_product( std::vector<double> a, std::vector< std::complex<double> > b ){
	return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}
std::complex<double> dot_product( std::vector< std::complex<double> > a, std::vector<double> b ){
	return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}
std::complex<double> dot_product( std::vector< std::complex<double> > a, std::vector< std::complex<double> > b ){
	return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

std::vector< std::complex<double> > vector_conj( std::vector< std::complex<double> > A ){
	// Return the complex conjugate of a vector A, defined as the vector of the complex conjugates of the elements of A.
	std::vector< std::complex<double> > ret ( A.size() );
	for( int i=0; i<A.size(); i++ ){
		ret[i] = std::conj( A[i] );
	}
	return ret;
}




double vector_magnitude( std::vector<double> v ){
	if( v.size() != 3 ){
		return 999;
	}
	return sqrt( pow( v[0], 2 ) + pow( v[1], 2 ) + pow( v[2], 2 ) );
}