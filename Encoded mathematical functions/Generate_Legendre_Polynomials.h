// polynomial_legendre_dx() and polynomial_legendre_dx2() require #include "Generate_Associated_Legendre_Functions". They are commented-out unless needed.

#include <math.h>
#include <vector>

//----- Global variables -----
//const double pi = 3.14159265358979323846;
std::vector<std::vector<double>> coeffs_legendre;								// Store precalculated polynomial coefficients for faster evaluation of Legendre polynomials.


//----- Functions for Legendre polynomials -----

void generate_coeffs_legendre( int order ){
	// Build arrays of Legendre polynomial coefficients, only if not enough already exist.
	if( coeffs_legendre.size() < order ){
		coeffs_legendre.clear();
		for( int ell=0; ell<=order; ell++ ){
			std::vector<double> coeffs(ell+1);
			for( int k=ell%2; k<=ell; k+=2 ){
				coeffs[k] = pow(2,ell) * tgamma(0.5*(k+1+ell)) / ( tgamma(k+1) * tgamma(ell-k+1) * tgamma(0.5*(k+1-ell) ) );
			}
			coeffs_legendre.push_back( coeffs );
		}
	}
}

double polynomial_legendre( double x, int ell ){
	if( ell < 0 ){
		return 0;
	}
	double ret = 0;
	for( int k=ell%2; k<=ell; k+=2 ){
		ret += pow(x,k) * coeffs_legendre[ell][k];
	}
	return ret;
}
/*double polynomial_legendre_dx( double x, int ell ){
	// Evaluate the associated Legendre function with indices ell,m at x.
	double ret = 0;
	for( int k=(ell+1)%2; k<=ell-1; k+=2 ){
		ret -= coeffs_associated_legendre[ell][ell+1][k] * pow(x,k);
	}
	return ret;
}
double polynomial_legendre_dx2( double x, int ell ){
	// Evaluate the associated Legendre function with indices ell,m at x.
	double ret = 0;
	for( int k=(ell+2)%2; k<=ell-2; k+=2 ){
		ret += coeffs_associated_legendre[ell][ell+2][k] * pow(x,k);
	}
	return ret;
}*/

double series_legendre( double x, std::vector<double> coeffs ){
	// Evaluate the legendre series of a univariate function of x, given a list of coefficients.
	double ret = 0;
	for( int ell=0; ell<coeffs.size(); ell++ ){
		ret += polynomial_legendre(x,ell) * coeffs[ell];
	}
	return ret;
}
double series_legendre( double x, double y, std::vector< std::vector<double> > coeffs ){
	// Evaluate the legendre series of a bivariate function of x and y, given a 2D list of coefficients.
	double ret = 0;
	for( int ell=0; ell<coeffs.size(); ell++ ){
		for( int m=0; m<coeffs[ell].size(); m++ ){
			ret += coeffs[ell][m] * polynomial_legendre( x, ell ) * polynomial_legendre( y, m );
		}
	}
	return ret;
}

std::vector<double> generate_series_coeffs_legendre( double (*f)(double), int ell_max, int n_steps=1000 ){
	std::vector<double> coeffs;
	double h   = 2 / ( (double) n_steps );
	for( int ell=0; ell<=ell_max; ell++ ){
		double x     = -1+h;
		double coeff = 0.5 * ( f(-1)*pow(-1,ell) + 4*f(-1+0.5*h)*polynomial_legendre(-1+0.5*h,ell) + f(1) );
		for( int i=1; i<n_steps; i++ ){
			coeff += f(x)*polynomial_legendre(x,ell) + 2*f(x+0.5*h)*polynomial_legendre(x+0.5*h,ell);
			x     += h;
		}
		coeffs.push_back( (ell+0.5) * coeff * h/3.0 );
	}
	return coeffs;
}

std::vector<double> generate_series_coeffs_legendre( std::vector<double> f, std::vector<double> x, int ell_max ){
	// Generate the Legendre series coefficients up to some ell_max for f(x), where both f and x are given as lists.
	std::vector<double> coeffs;
	for( int ell=0; ell<=ell_max; ell++ ){
		double coeff = 0;
		for( int i=0; i<f.size()-1; i++ ){
			coeff += ( f[i]*polynomial_legendre(x[i],ell) + f[i+1]*polynomial_legendre(x[i+1],ell) ) * ( x[i+1] - x[i] );
		}
		coeffs.push_back( (ell+0.5) * 0.5 * coeff );
	}
	return coeffs;
}





std::vector<double> generate_series_coeffs_legendre_constant_stepsize( std::vector<double> f, std::vector<double> x, int ell_max ){
	// Generate the Legendre series coefficients up to some ell_max for f(x), where both f and x are given as lists.
	// This version is more efficient, but requires the stepsize x[i+1]-x[i] to be constant.
	std::vector<double> coeffs;
	for( int ell=0; ell<=ell_max; ell++ ){
		double coeff = 0.5 * ( f[0] * pow(-1,ell) + f.back() );
		for( int i=1; i<f.size()-1; i++ ){
			coeff += f[i] * polynomial_legendre(x[i],ell);
		}
		coeffs.push_back( (ell+0.5) * ( x[1] - x[0] ) * coeff );
	}
	return coeffs;
}




/*std::vector< std::vector<double> > generate_series_coeffs_legendre( double (*f)(double,double), int ell_max, int m_max, int n_steps_x, int n_steps_y ){
	// Generate the coefficients of the Legendre series of a continuous bivariate function f(x,y) up to maximum indices ell_max and m_max, where ell refers to the indices for x and m to y. Return the coefficients as a 2D vector.
	std::vector< std::vector<double> > coeffs( ell_max+1, std::vector<double> ( m_max+1 ) );
	double hx  = 2.0 / ( (double) n_steps_x );
	double hy  = 2.0 / ( (double) n_steps_y );
	
	for( int ell=0; ell<=ell_max; ell++ ){
		double P_ell_m1 = pow( -1, ell );
		for( int m=0; m<=m_max; m++ ){
			
			double P_m_m1 = pow( -1, m );
			double coeff  = 0;
			double x      = -1 + hx;
			double y      = -1 + hy;
			
			//--- Double sum, and single sum over x ---
			for( int i=1; i<n_steps_x; i++ ){
				y = -1 + hy;
				coeff += 0.5 * ( ( f(x,-1)*P_m_m1 + f(x,1) ) * polynomial_legendre(x,ell) );	// Single sum over x
				for( int j=1; j<n_steps_y; j++ ) {
					coeff += f(x,y)*polynomial_legendre(x,ell)*polynomial_legendre(y,m);		// Double sum
					y     += hy;
				}
				x += hx;
			}
			
			//--- Single sum over y ---
			y = -1 + hy;
			for( int i=1; i<n_steps_y; i++ ){
				coeff += 0.5 * ( ( f(-1,y)*P_ell_m1 + f(1,y) ) * polynomial_legendre(y,m) );
				y     += hy;
			}
			
			//--- Vertices ---
			coeff += 0.25 * ( f(-1,-1)*P_ell_m1*P_m_m1 + f(-1,1)*P_ell_m1 + f(1,-1)*P_m_m1 + f(1,1) );
			
			coeffs[ell][m] = (ell+0.5) * (m+0.5) * hx * hy * coeff;
		}
	}
	return coeffs;
}*/