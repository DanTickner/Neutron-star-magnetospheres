#include <math.h>
#include <vector>
#include <complex>

//----- Global variables -----
const double pi = acos(-1);										// Need this only if we're calling this header file on its own. Otherwise pi defined elsewhere.
std::vector< std::vector< std::vector<double> > > coeffs_sh;					//Store precalculated coefficients for faster spherical harmonic evaluation.


//----- Functions for spherical harmonics -----

void generate_coeffs_spherical_harmonic( int ell_max ){
	coeffs_sh.clear();
	for( int ell=0; ell<=ell_max; ell++ ){
		std::vector< std::vector< double > > coeffs_sh_for_this_ell;
		for( int m=-ell; m<=ell; m++ ){
			int M = abs( m );
			std::vector<double> coeffs_sh_for_this_ell_and_m( ell-M+1, 0 );
			double factor = sqrt( (2*ell+1) * tgamma(ell-M+1) / ( 4*pi * tgamma(ell+M+1) )) * pow(2,ell);
			for( int k=(ell+M)%2; k<=ell-M; k+=2 ){
				coeffs_sh_for_this_ell_and_m[k] =  factor * tgamma(0.5*(k+M+1+ell)) / ( tgamma(0.5*(k+M+1-ell)) * tgamma(ell-k-M+1) * tgamma(k+1) );
				if( m > 0 ){
					coeffs_sh_for_this_ell_and_m[k] *= pow( -1, m );
				}
			}
			coeffs_sh_for_this_ell.push_back( coeffs_sh_for_this_ell_and_m );
		}
		coeffs_sh.push_back( coeffs_sh_for_this_ell );
	}
}




std::complex<double> ylm( double theta, double phi, int ell, int m ){
	double ret = 0;
	double x   = cos( theta );
	for( int k=(ell+m)%2; k<=ell-abs(m); k+=2 ){
		ret += coeffs_sh[ell][ell+m][k] * pow(x,k);
	}
	return ret * exp(std::complex<double>{0,m*phi}) * pow(sin(theta),abs(m));
}




std::complex<double> ylm_over_sintheta( double theta, double phi, int ell, int m ){
	// The cases m=0 and m=+-1 are not yet covered.
	double ret = 0;
	double x   = cos( theta );
	for( int k=(ell+m)%2; k<=ell-abs(m); k+=2 ){
		ret += coeffs_sh[ell][ell+m][k] * pow(x,k);
	}
	return ret * exp(std::complex<double>{0,m*phi}) * pow(sin(theta),abs(m)-1);
}




std::complex<double> ylm_dtheta( double theta, double phi, int ell, int m ){
	if( abs(m) != ell ){
		return 0.5*( sqrt((ell-m)*(ell+m+1)) * exp(std::complex<double>{0,-phi}) * ylm(theta,phi,ell,m+1) - sqrt((ell+m)*(ell-m+1)) * exp(std::complex<double>{0,phi}) * ylm(theta,phi,ell,m-1) );
	}
	if( ell == 0 ){
		return 0;
	}
	if( m == ell ){
		return - exp(std::complex<double>{0,phi}) * sqrt(0.5*ell) * ylm(theta,phi,ell,ell-1);
	}
	return exp(std::complex<double>{0,-phi}) * sqrt(0.5*ell) * ylm(theta,phi,ell,-ell+1); // Only remaining possibility is m = -ell.
}




std::complex<double> ylm_dtheta2( double theta, double phi, int ell, int m ){
	if( ( abs(m) != ell ) && ( abs(m) != ell-1 ) ){
		return std::complex<double>{cos(2*phi),-sin(2*phi)} * 0.25 * sqrt((ell+m+2)*(ell+m+1)*(ell-m)*(ell-m-1)) * ylm(theta,phi,ell,m+2) + exp(std::complex<double>{0,2*phi}) * 0.25 * sqrt((ell+m)*(ell+m-1)*(ell-m+2)*(ell-m+1)) * ylm(theta,phi,ell,m-2) + 0.5*(m*m-ell*(ell+1)) * ylm(theta,phi,ell,m);
	}
	if( ell == 1 ){
		return - ylm(theta,phi,ell,m);
	}
	if( ell == 0 ){
		return 0;
	}	
	if( m == ell-1 ){
		if( ell == 2 ){
			return exp(std::complex<double>{0,phi}) * sqrt(7.5/pi) * sin(2*theta);
		}
		return -((double)ell*ell) * ylm(theta,phi,ell,ell-1) + exp(std::complex<double>{0,2*phi}) * 0.5 * sqrt((2*ell+1)*(2*ell-1)*(ell-1)*(ell-2)) * ylm(theta,phi,ell-2,ell-3);
	}
	if( m == -ell+1 ){
		if( ell == 2 ){
			return - exp(std::complex<double>{0,-phi}) * sqrt(7.5/pi) * sin(2*theta);
		}
		return std::complex<double>{cos(2*(ell-1)*phi),-sin(2*(ell-1)*phi)} * pow(-1,ell) * ((double)ell*ell) * ylm(theta,phi,ell,ell-1) + std::complex<double>{cos(2*(ell-2)*phi),-sin(2*(ell-2)*phi)} * pow(-1,ell+1) * 0.5 * sqrt((2*ell+1)*(2*ell-1)*(ell-1)*(ell-2)) * ylm(theta,phi,ell-2,ell-3);
	}
	if( m == ell ){
		return -((double)ell*ell) * ylm(theta,phi,ell,ell) + exp(std::complex<double>{0,2*phi}) * 0.5 * sqrt(ell*(ell-1)*(2*ell+1)*(2*ell-1)) * ylm(theta,phi,ell-2,ell-2);
	}
	return std::complex<double>{cos(2*ell*phi),-sin(2*ell*phi)} * pow(-1,ell+1) * ((double)ell*ell) * ylm(theta,phi,ell,ell) + std::complex<double>{cos((2*ell-2)*phi),-sin((2*ell-2)*phi)} * pow(-1,ell) * 0.5 * sqrt(ell*(ell-1)*(2*ell+1)*(2*ell-1)) * ylm(theta,phi,ell-2,ell-2);	// Only remaining possibility is m = -ell.
}




std::complex<double> ylm_dphi( double theta, double phi, int ell, int m ){
	// Evaluate the first phi-derivative of the spherical harmonic of a given ell, m at a given theta, phi.
	if( m == 0 ) {
		return std::complex<double> (0,0);
	}
	return std::complex<double>(0,m) * ylm( theta, phi, ell, m );
}

std::complex<double> ylm_dphi2( double theta, double phi, int ell, int m ){
	// Evaluate the second phi-derivative of the spherical harmonic of a given ell, m at a given theta, phi. int * complex not defined.
	if( m == 0 ) {
		return std::complex<double> (0,0);
	}
	return std::complex<double>(-m*m,0) * ylm( theta, phi, ell, m );
}

std::complex<double> ylm_dphi_over_sintheta( double theta, double phi, int ell, int m ){
	// Evaluate 1/sin(theta) times the first phi-derivative of the spherical harmonic of a given ell, m at a given theta, phi.
	// Used extensively in vector spherical harmonics.
	if( m == 0 ) {
		return std::complex<double> (0,0);
	}
	std::complex<double> ret_a = sqrt( (ell+m+2.0)*(ell+m+1.0) ) * ylm( theta, phi, ell+1, m+1 ) * exp( std::complex<double>{0,-phi} );
	std::complex<double> ret_b = sqrt( (ell-m+2.0)*(ell-m+1.0) ) * ylm( theta, phi, ell+1, m-1 ) * exp( std::complex<double>{0, phi} );
	return -0.5 * std::complex<double>{0,1} * sqrt( (2.0*ell+1.0)/((double)(2.0*ell+3.0)) ) * ( ret_a + ret_b );
}




//----- Vector spherical harmonics -----
std::vector< std::complex<double> > VSH_Y( double theta, double phi, int ell, int m ){
	return std::vector< std::complex<double> > { ylm(theta,phi,ell,m), 0, 0 };
}

std::vector< std::complex<double> > VSH_Psi( double theta, double phi, int ell, int m ){
	return std::vector< std::complex<double> > { 0, ylm_dtheta(theta,phi,ell,m), ylm_dphi_over_sintheta(theta,phi,ell,m) };
}

std::vector< std::complex<double> > VSH_Phi( double theta, double phi, int ell, int m ){
	return std::vector< std::complex<double> > { 0, -ylm_dphi_over_sintheta(theta,phi,ell,m), ylm_dtheta(theta,phi,ell,m) };
}




//----- Wigner 3j symbols -----
double wigner_3j( int L1, int L2, int L3, int M1, int M2, int M3 ){
	// https://en.wikipedia.org/wiki/3-j_symbol#Explicit_expression
	
	bool condition_M = ( M1+M2+M3 != 0 );
	bool condition_triangle = ( ( L1+L2-L3 +1 < 0 ) or ( L1-L2+L3 < 0 ) or ( -L1+L2+L3 < 0 ) or ( L1+L2+L3+1 < 0 ) );
	bool condition_LM = ( ( L1-M1 < 0 ) or ( L1+M1 < 0 ) or ( L2-M2 < 0 ) or ( L2+M2 < 0 ) or ( L3-M3 < 0 ) or ( L3+M3 < 0 ) );
	
	if( condition_M or condition_triangle or condition_LM ){
		return 0;
	}
	
	
	double ret = 0;
	
	int k_min = std::max( std::max( 0, -L3+L2-M1 ), -L3+L1+M2 );
	int k_max = std::min( std::min( L1+L2-L3, L1-M1 ), L2+M2 );
	
	for( int k=k_min; k<=k_max; k++ ){
		ret += pow(-1,k) / ( tgamma( k +1 ) * tgamma( L1+L2-L3-k +1 ) * tgamma( L1-M1-k +1 ) * tgamma( L2+M2-k +1 ) * tgamma( L3-L2+M1+k +1 ) * tgamma( L3-L1-M2+k +1 ) );
	}
	
	ret *= sqrt( tgamma( L1-M1 +1 ) * tgamma( L1+M1 +1 ) * tgamma( L2-M2 +1 ) * tgamma( L2+M2 +1 ) * tgamma( L3-M3 +1 ) * tgamma( L3+M3 + 1 ) );
	ret *= sqrt( tgamma( L1+L2-L3 +1 ) * tgamma( L1-L2+L3 +1 ) * tgamma( -L1+L2+L3 +1 ) / tgamma( L1+L2+L3+1 + 1 ) );
	ret *= pow( -1, L1-L2-M3 );
	
	return ret;
}

double wigner_3j_000( int L1, int L2, int L3 ){
	// https://en.wikipedia.org/wiki/3-j_symbol#Explicit_expression
	
	bool condition_triangle = ( ( L1+L2-L3 +1 < 0 ) or ( L1-L2+L3 < 0 ) or ( -L1+L2+L3 < 0 ) or ( L1+L2+L3+1 < 0 ) );
	bool condition_L = ( ( L1 < 0 ) or ( L2 < 0 ) or ( L3 < 0 ) );
	
	if( condition_triangle or condition_L ){
		return 0;
	}
	
	
	double ret = 0;
	
	int k_min = std::max( std::max( 0, L2-L3 ), L1-L3 );
	int k_max = std::min( std::min( L1+L2-L3, L1 ), L2 );
	
	for( int k=k_min; k<=k_max; k++ ){
		ret += pow(-1,k) / ( tgamma( k +1 ) * tgamma( L1+L2-L3-k +1 ) * tgamma( L1-k +1 ) * tgamma( L2-k +1 ) * tgamma( L3-L2+k +1 ) * tgamma( L3-L1+k +1 ) );
	}
	
	ret *= tgamma( L1 +1 ) * tgamma( L2 +1 ) * tgamma( L3 +1 );
	ret *= sqrt( tgamma( L1+L2-L3 +1 ) * tgamma( L1-L2+L3 +1 ) * tgamma( -L1+L2+L3 +1 ) / tgamma( L1+L2+L3+1 + 1 ) );
	ret *= pow( -1, L1-L2 );
	
	return ret;
}