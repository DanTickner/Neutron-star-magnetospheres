#include <math.h>
#include <vector>

//----- Global variables -----
const double pi = 3.14159265358979323846;
std::vector<std::vector<int   >> coeffs_chebyshev_1;							//Store precalculated polynomial coefficients for faster evaluation of Chebyshev polynomials of the first kind.
std::vector<std::vector<int   >> coeffs_chebyshev_2;							//Store precalculated polynomial coefficients for faster evaluation of Chebyshev polynomials of the second kind.
std::vector<std::vector<int   >> coeffs_chebyshev_1_dr;							//Store precalculated polynomial coefficients for faster evaluation of first derivatives of Chebyshev polynomials of the first kind.
std::vector<std::vector<int   >> coeffs_chebyshev_2_dr;							//Store precalculated polynomial coefficients for faster evaluation of first derivatives of Chebyshev polynomials of the second kind.
std::vector<std::vector<int   >> coeffs_chebyshev_1_dr2;						//Store precalculated polynomial coefficients for faster evaluation of second derivatives of Chebyshev polynomials of the first kind.
std::vector<std::vector<int   >> coeffs_chebyshev_2_dr2;						//Store precalculated polynomial coefficients for faster evaluation of second derivatives of Chebyshev polynomials of the second kind.


//----- Functions for Chebyshev polynomials -----

void generate_coeffs_chebyshev( int order ){
	order = std::max( order, 1 );	// Entering order=0 causes nothing to happen.
	//----- Build arrays of Chebyshev polynomial coefficients, only if not enough already exist -----
	if( coeffs_chebyshev_1.size() < order ){
		coeffs_chebyshev_1     = { std::vector<int>{1}, std::vector<int>{0,1} };
		coeffs_chebyshev_2     = { std::vector<int>{1}, std::vector<int>{0,2} };
		coeffs_chebyshev_1_dr  = { std::vector<int>{0}                        };		// In this and the following, we wust specify first vector because the for-loop will have 1-1=0 terms.
		coeffs_chebyshev_2_dr  = { std::vector<int>{0}                        };
		coeffs_chebyshev_1_dr2 = { std::vector<int>{0}                        };
		coeffs_chebyshev_2_dr2 = { std::vector<int>{0}                        };
		
		//----- Build arrays of Chebyshev polynomial coefficients (first and second kind) -----
		for( int n=0; n<order+1; n++ ){
			int index = coeffs_chebyshev_1.size() - 2;
			coeffs_chebyshev_1[index].push_back( 0 );
			coeffs_chebyshev_2[index].push_back( 0 );
			std::vector<int> coeffs_new_1 = { -coeffs_chebyshev_1[index][0] };
			std::vector<int> coeffs_new_2 = { -coeffs_chebyshev_2[index][0] };
			for( int i=0; i<coeffs_chebyshev_1[index].size()-1; i++ ){
				coeffs_new_1.push_back( 2*coeffs_chebyshev_1.back()[i] - coeffs_chebyshev_1[index][i+1] );
				coeffs_new_2.push_back( 2*coeffs_chebyshev_2.back()[i] - coeffs_chebyshev_2[index][i+1] );
			}
			coeffs_new_1.push_back( 2 * coeffs_chebyshev_1.back().back() );
			coeffs_chebyshev_1[index].pop_back();
			coeffs_chebyshev_1.push_back( coeffs_new_1 );
			coeffs_new_2.push_back( 2 * coeffs_chebyshev_2.back().back() );
			coeffs_chebyshev_2[index].pop_back();
			coeffs_chebyshev_2.push_back( coeffs_new_2 );
		}
		
		//----- Build array of derivatives of Chebyshev polynomial coefficients (first and second kind) -----
		for( int n=1; n<order+1; n++ ){
			std::vector<int> coeffs_dr;
			for( int i=0; i<coeffs_chebyshev_1[n].size()-1; i++ ){
				coeffs_dr.push_back( ( i + 1 ) * coeffs_chebyshev_1[n][i+1] );
			}
			coeffs_chebyshev_1_dr.push_back( coeffs_dr );
		}
		for( int n=1; n<order+1; n++ ){
			std::vector<int> coeffs_dr2;
			for( int i=0; i<coeffs_chebyshev_1_dr[n].size()-1; i++ ){
				coeffs_dr2.push_back( ( i + 1 ) * coeffs_chebyshev_1_dr[n][i+1] );
			}
			coeffs_chebyshev_1_dr2.push_back( coeffs_dr2 );
		}
		for( int n=1; n<order+1; n++ ){
			std::vector<int> coeffs_dr;
			for( int i=0; i<coeffs_chebyshev_2[n].size()-1; i++ ){
				coeffs_dr.push_back( ( i + 1 ) * coeffs_chebyshev_2[n][i+1] );
			}
			coeffs_chebyshev_2_dr.push_back( coeffs_dr );
		}
		for( int n=1; n<order+1; n++ ){
			std::vector<int> coeffs_dr2;
			for( int i=0; i<coeffs_chebyshev_2_dr[n].size()-1; i++ ){
				coeffs_dr2.push_back( ( i + 1 ) * coeffs_chebyshev_2_dr[n][i+1] );
			}
			coeffs_chebyshev_2_dr2.push_back( coeffs_dr2 );
		}
	}
}




double polynomial_chebyshev_1( double r, int order ){
	// Evaluate the Chebyshev polynomial of the first kind to a given order and a given value of r.
	double ret = 0;
	double rn  = pow( r, order%2 );	// r^n
	for( int i=order%2; i<=coeffs_chebyshev_1[order].size(); i+=2 ){
		ret += coeffs_chebyshev_1[order][i] * rn;
		rn  *= r*r;
	}
	return ret;
}




double polynomial_chebyshev_2( double r, int order ){
	// Evaluate the Chebyshev polynomial of the second kind to a given order and a given value of r.
	double ret = 0;
	double rn  = pow( r, order%2 );	// r^n
	for( int i=order%2; i<=coeffs_chebyshev_2[order].size(); i+=2 ){
		ret += coeffs_chebyshev_2[order][i] * rn;
		rn  *= r*r;
	}
	return ret;
}




double polynomial_chebyshev_1_dr( double r, int order ){
	// Evaluate the first derivative of the Chebyshev polynomial of the first kind to a given order and a given value of r.
	double ret = 0;
	double rn  = pow( r, (order-1)%2 );	// r^n
	for( int i=(order-1)%2; i<=coeffs_chebyshev_1_dr[order].size(); i+=2 ){
		ret += coeffs_chebyshev_1_dr[order][i] * rn;
		rn  *= r*r;
	}
	return ret;
}




double polynomial_chebyshev_1_dr2( double r, int order ){
	// Evaluate the second derivative of the Chebyshev polynomial of the first kind to a given order and a given value of r.
	double ret = 0;
	double rn  = pow( r, order%2 );	// r^n
	for( int i=order%2; i<=coeffs_chebyshev_1_dr2[order].size(); i+=2 ){	// Really (order-2)%2 but this is equivalent.
		ret += coeffs_chebyshev_1_dr2[order][i] * rn;
		rn  *= r*r;
	}
	return ret;
}




double polynomial_chebyshev_2_dr( double r, int order ){
	// Evaluate the first derivative of the Chebyshev polynomial of the second kind to a given order and a given value of r.
	double ret = 0;
	double rn  = pow( r, order%2 );	// r^n
	for( int i=order%2; i<=coeffs_chebyshev_2_dr[order].size(); i+=2 ){
		ret += coeffs_chebyshev_2_dr[order][i] * rn;
		rn  *= r*r;
	}
	return ret;
}




double polynomial_chebyshev_2_dr2( double r, int order ){
	// Evaluate the second derivative of the Chebyshev polynomial of the second kind to a given order and a given value of r.
	double ret = 0;
	double rn  = pow( r, order%2 );	// r^n
	for( int i=order%2; i<=coeffs_chebyshev_2_dr2[order].size(); i+=2 ){
		ret += coeffs_chebyshev_2_dr2[order][i] * rn;
		rn  *= r*r;
	}
	return ret;
}




std::vector<double> generate_series_coeffs_chebyshev( double (*f)(double), int n_max, int n_steps=1000 ){
	std::vector<double> coeffs;
	double h   = pi / ( (double) n_steps );
	for( int n=0; n<=n_max; n++ ){
		double coeff = 0.5 * ( f(1) + 4*f(cos(0.5*h))*cos(n*0.5*h) + f(-1)*pow(-1,n) );
		double t     = h;
		for( int i=1; i<n_steps; i++ ){
			coeff += f(cos(t))*cos(n*t) + 2*f(cos(t+0.5*h))*cos(n*(t+0.5*h));
			t     += h;
		}
		coeffs.push_back( 2.0/pi * coeff * h/3.0 );
	}
	return coeffs;
}




double series_chebyshev( double x, std::vector<double> coeffs ){
	double ret = 0.5 * coeffs[0];
	for( int n=1; n<coeffs.size(); n++ ){
		ret += coeffs[n] * polynomial_chebyshev_1(x,n);
	}
	return ret;
}