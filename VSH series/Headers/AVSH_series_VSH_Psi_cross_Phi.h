/*
AVSH_series_VSH_Psi_cross_Phi.h

*/

#include "../../Generate_Spherical_Harmonics.h"			// VSH and wigner_3j()


int L1 = 3;
int L2 = 4;
	
	
	
double A_r_function( double r, double theta ){
	generate_coeffs_spherical_harmonic( std::max( L1, L2 ) );
	return cross_product( VSH_Psi( theta, 0, L1, 0 ), VSH_Phi( theta, 0, L2, 0 ) )[0].real();
}


double A_theta_function( double r, double theta ){
	return 0;
}

double A_phi_function( double r, double theta ){
	return 0;
}




double I_011( int L1, int L2, int L3 ){
	// 1/(4sqrt(pi)) sqrt(2L1+1) sqrt(L2(L2+1)(2L2+1)) sqrt(L3(L3+1)(2L3+1)) int_0^pi P_{L1}^0 P_{L2}^1 P_{L3}^1 sin(theta) dtheta
	// ../Check integrals/20230807_P_L1_0_P_L2_1_P_L_1.cpp
	
	double ret = 0;
	
	for( int L12=std::max(abs(L1-L2),1); L12<=L1+L2; L12++ ){
		if( (L1+L2+L12)%2 == 0 ){
		
			double G12 = - ( 2*L12+1 ) * wigner_3j_000( L1, L2, L12 ) * wigner_3j( L1, L2, L12, 0, 1, -1 );
			
			if( (L12+L3)%2 == 0 ){
				for( int L123=std::max(abs(L12-L3),2); L123<=L12+L3; L123++ ){
					if( L123 %2 == 0 ){
						double G123          = ( 2*L123+1 ) * wigner_3j_000( L12, L3, L123 ) * wigner_3j( L12, L3, L123, 1, 1, -2 );
						double sqrt_fraction = pow( (L123+2) * (L123+1) * L123 * (L123-1), -0.5 );
						
						ret += G12 * G123 * sqrt_fraction;
					}
				}
			}
			
		}
	}
	
	return ret * sqrt( (2*L1+1) * L2*(L2+1)*(2*L2+1) * L3*(L3+1)*(2*L3+1) / pi );
}




double A_r_ell_guess( double r, int L ){
	
	/*
	//----- V1 -----
	if( L > L1+L2 ){
		return 0;
	}
	
	double ret = 0;
	
	for( int L12=std::max(abs(L1-L2),2); L12<=L1+L2; L12++ ){
		if( (L1+L2+L12)%2 == 0 ){
		
			double G12 = ( 2*L12+1 ) * wigner_3j_000( L1, L2, L12 ) * wigner_3j( L1, L2, L12, 1, 1, -2 );
			
			if( (L12+L)%2 == 0 ){
				for( int L123=std::max(abs(L12-L),2); L123<=L12+L; L123++ ){
					if( L123 %2 == 0 ){
						double G123          = ( 2*L123+1 ) * wigner_3j_000( L12, L, L123 ) * wigner_3j( L12, L, L123, 2, 0, -2 );
						double sqrt_fraction = pow( (L123+2) * (L123+1) * L123 * (L123-1), -0.5 );
						
						ret += G12 * G123 * sqrt_fraction;
					}
				}
			}
			
		}
	}
	
	return ret * sqrt( (2*L+1)*(2*L1+1)*(2*L2+1) * L1*(L1+1) * L2*(L2+1) / pi );
	*/
	
	//----- V2 -----
	if( L==0 ){
		return 0;
	}
	if( L > L1+L2 ){
		return 0;
	}
	return I_011( L, L1, L2 );
}

double A_1_ell_guess( double r, int ell ){
	return 0;
}

double A_2_ell_guess( double r, int ell ){
	return 0;
}