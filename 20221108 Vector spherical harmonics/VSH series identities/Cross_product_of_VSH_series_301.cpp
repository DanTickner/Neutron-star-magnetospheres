// cd OneDrive\PhD\Codes\20221108 Vector spherical harmonics\VSH series identities
// g++ Cross_product_of_VSH_series_301.cpp -o Cross_product_of_VSH_series_301

/*
Fit vector spherical harmonic series to two vectors A(theta,phi) and B(theta,phi) and calculate their cross product.
Compare this to an analytical expression derived on paper (Book 10, p58).
No longer force as many terms as possible into ell_1 and m_1. Hence, no need for k-indices.
*/

#include <iostream>
#include <iomanip>
#include <chrono>
#include "../../Generate_Spherical_Harmonics.h"
//#include "../../Vector_Operations.h"	Needed but called within Generate_Spherical_Harmonics.h
#include "../Headers_vector_function/Vector_function_05.h"
#include "../Headers_vector_function/Vector_function_06.h"
#include "../Headers_VSH_series/VSH_series_calculation_101.h"
#include "../Headers_VSH_series/VSH_series_evaluation_101.h"

int main(){
	
	auto time_start = std::chrono::high_resolution_clock::now();
	
	//----- Define variables -----
	int ell_max = 3;
	int N_t     = 30;
	int N_p     = 30;
	int i_test  = 18;
	int j_test  = 20;
	
	std::vector<double> t_list ( N_t );
	std::vector<double> p_list ( N_p );
	
	std::vector< std::vector< std::vector< std::complex<double> > > > A( 3, std::vector< std::vector< std::complex<double> > > ( N_t, std::vector< std::complex<double> > ( N_p ) ) );		// Precalculated values of the function that we wish to fit.
	std::vector< std::vector< std::vector< std::complex<double> > > > B( 3, std::vector< std::vector< std::complex<double> > > ( N_t, std::vector< std::complex<double> > ( N_p ) ) );		// Precalculated values of the function that we wish to fit.
	
	// Totals of absolute values of fitted coefficients in each direction. No physical interpretation but useful as a quick verification tool.
	std::vector< std::complex<double> > coeffs_total_A ( 3 );
	std::vector< std::complex<double> > coeffs_total_B ( 3 );
	
	
	//----- Fail if the test indices are greater than the maximum number of points -----
	if( i_test > N_t-1 ){
		std::cout << "Test theta index (" << i_test << ") exceeds number of datapoints (" << N_t << "). Stopping the code to prevent undefined behaviour." << std::endl;
		return 1;
	}
	if( j_test > N_p-1 ){
		std::cout << "Test phi index (" << j_test << ") exceeds number of datapoints (" << N_p << "). Stopping the code to prevent undefined behaviour." << std::endl;
		return 1;
	}
	
	
	//----- Build coordinate lists -----
	for( int i=0; i<N_t; i++ ){
		t_list[i] = i * pi / ( (double) N_t - 1 );
	}
	for( int j=0; j<N_p; j++ ){
		p_list[j] = j * 2*pi / ( (double) N_p - 1 );
	}
	
	
	//----- Generate spherical harmonic coefficients -----
	generate_coeffs_spherical_harmonic( ell_max );
	
	
	//----- Precalculate vector A to be fitted -----
	for( int i=0; i<N_t; i++ ){
		for( int j=0; j<N_p; j++ ){
			
			A[0][i][j] = A_r_function( t_list[i], p_list[j] );
			A[1][i][j] = A_t_function( t_list[i], p_list[j] );
			A[2][i][j] = A_p_function( t_list[i], p_list[j] );
			
			B[0][i][j] = B_r_function( t_list[i], p_list[j] );
			B[1][i][j] = B_t_function( t_list[i], p_list[j] );
			B[2][i][j] = B_p_function( t_list[i], p_list[j] );
		}
	}
	
	
	//----- Calculate VSH series -----
	std::vector< std::vector< std::vector< std::complex<double> > > > coeffs_A = generate_vsh_series( A, t_list, p_list, ell_max );
	std::vector< std::vector< std::vector< std::complex<double> > > > coeffs_B = generate_vsh_series( B, t_list, p_list, ell_max );
	
	
	//----- Output VSH series -----
	for( int ell=0; ell<=ell_max; ell++ ){
		for( int m=-ell; m<=ell; m++ ){
			for( int c=0; c<3; c++ ){
				coeffs_total_A[c] += std::complex<double> { abs( coeffs_A[c][ell][ell+m].real() ), abs( coeffs_A[c][ell][ell+m].imag() ) };
				coeffs_total_B[c] += std::complex<double> { abs( coeffs_B[c][ell][ell+m].real() ), abs( coeffs_B[c][ell][ell+m].imag() ) };
			}
		}
	}
	
	std::cout << "Coeffs total A:\t" << coeffs_total_A[0] <<"\t"<< coeffs_total_A[1] <<"\t"<< coeffs_total_A[2] <<"\t\t"<< coeffs_total_A[0]+coeffs_total_A[1]+coeffs_total_A[2] << std::endl;
	std::cout << "Coeffs total B:\t" << coeffs_total_B[0] <<"\t"<< coeffs_total_B[1] <<"\t"<< coeffs_total_B[2] <<"\t\t"<< coeffs_total_B[0]+coeffs_total_B[1]+coeffs_total_B[2] << std::endl;
	
	
	//----- Evaluate cross product at the test point -----
	double t = t_list[i_test];
	double p = p_list[j_test];
	
	std::cout << "\nTest fit for (theta,phi) = (" << t <<","<< p << "):" << std::endl;
	
	std::vector< std::complex<double> > A_test_calc  = evaluate_vsh_series( t, p, coeffs_A );
	std::vector< std::complex<double> > B_test_calc  = evaluate_vsh_series( t, p, coeffs_B );
	
	std::vector< std::complex<double> > A_test_exact{ A[0][i_test][j_test], A[1][i_test][j_test], A[2][i_test][j_test] };
	std::vector< std::complex<double> > B_test_exact{ B[0][i_test][j_test], B[1][i_test][j_test], B[2][i_test][j_test] };
	
	std::vector< std::complex<double> > cross_product_exact{ A_test_exact[1]*B_test_exact[2]-A_test_exact[2]*B_test_exact[1],
	                                                         A_test_exact[2]*B_test_exact[0]-A_test_exact[0]*B_test_exact[2],
															 A_test_exact[0]*B_test_exact[1]-A_test_exact[1]*B_test_exact[0] };
	
	std::vector< std::complex<double> > cross_product_from_VSH_series{ A_test_calc[1]*B_test_calc[2]-A_test_calc[2]*B_test_calc[1],
	                                                                   A_test_calc[2]*B_test_calc[0]-A_test_calc[0]*B_test_calc[2],
															           A_test_calc[0]*B_test_calc[1]-A_test_calc[1]*B_test_calc[0] };
	
	print_vector( cross_product_exact          , "Cross product exact          " );
	print_vector( cross_product_from_VSH_series, "Cross product from VSH series" );
	
	
	//----- Calculate cross product from VSH coefficients -----
	int k_max = pow(ell_max+1,2) - 1;
	std::vector< std::complex<double> > cross_product_from_VSH_coeffs   { 0, 0, 0 };
	std::vector< std::complex<double> > cross_product_from_VSH_coeffs_a { 0, 0, 0 };
	
	std::vector< std::complex<double> > cross_product_from_VSH_coeffs_01 { 0, 0, 0 };
	std::vector< std::complex<double> > cross_product_from_VSH_coeffs_02 { 0, 0, 0 };
	std::vector< std::complex<double> > cross_product_from_VSH_coeffs_10 { 0, 0, 0 };
	std::vector< std::complex<double> > cross_product_from_VSH_coeffs_11 { 0, 0, 0 };
	std::vector< std::complex<double> > cross_product_from_VSH_coeffs_12 { 0, 0, 0 };
	std::vector< std::complex<double> > cross_product_from_VSH_coeffs_20 { 0, 0, 0 };
	std::vector< std::complex<double> > cross_product_from_VSH_coeffs_21 { 0, 0, 0 };
	std::vector< std::complex<double> > cross_product_from_VSH_coeffs_22 { 0, 0, 0 };
	
	
	//----- Calculate cross product -----
	//--- V1 --- 
	/*for( int L1=0; L1<=ell_max; L1++ ){
		for( int M1=-L1; M1<=L1; M1++ ){
			
			std::vector< std::complex<double> > Ar12 { coeffs_A[0][L1][L1+M1], coeffs_A[1][L1][L1+M1], coeffs_A[2][L1][L1+M1] };
			std::vector< std::complex<double> > Psi1 = VSH_Psi(t,p,L1,M1);
			std::vector< std::complex<double> > Phi1 = VSH_Phi(t,p,L1,M1);
			
			for( int L2=0; L2<=ell_max; L2++ ){
				for( int M2=-L2; M2<=L2; M2++ ){
					
					std::vector< std::complex<double> > Br12 { coeffs_B[0][L2][L2+M2], coeffs_B[1][L2][L2+M2], coeffs_B[2][L2][L2+M2] };
					std::vector< std::complex<double> > Psi2 = VSH_Psi(t,p,L2,M2);
					std::vector< std::complex<double> > Phi2 = VSH_Phi(t,p,L2,M2);
					
					//--- Radial component ---
					std::complex<double> coeff_Psi = Ar12[1]*Br12[2] - Ar12[2]*Br12[1];
					std::complex<double> coeff_Phi = Ar12[1]*Br12[1] + Ar12[2]*Br12[2];
					std::vector< std::complex<double> > coeff { 0, coeff_Psi*Psi1[1]+coeff_Phi*Phi1[1], coeff_Psi*Psi1[2]+coeff_Phi*Phi1[2] };
					cross_product_from_VSH_coeffs_a[0] += coeff[1]*Psi2[1] + coeff[2]*Psi2[2];
					
					//--- Angular components ---
					cross_product_from_VSH_coeffs_a = vector_sum( cross_product_from_VSH_coeffs_a, scalar_times_vector(  Ar12[2]*Br12[0]*ylm(t,p,L2,M2), Psi1 ) );
					cross_product_from_VSH_coeffs_a = vector_sum( cross_product_from_VSH_coeffs_a, scalar_times_vector( -Ar12[1]*Br12[0]*ylm(t,p,L2,M2), Phi1 ) );
					cross_product_from_VSH_coeffs_a = vector_sum( cross_product_from_VSH_coeffs_a, scalar_times_vector( -Ar12[0]*Br12[2]*ylm(t,p,L1,M1), Psi2 ) );
					cross_product_from_VSH_coeffs_a = vector_sum( cross_product_from_VSH_coeffs_a, scalar_times_vector(  Ar12[0]*Br12[1]*ylm(t,p,L1,M1), Phi2 ) );
				}
			}
		}
	}*/
	
	//--- V2 ---
	/*for( int L1=0; L1<=ell_max; L1++ ){
		for( int M1=-L1; M1<=L1; M1++ ){
			
			std::vector< std::complex<double> > Ar12 { coeffs_A[0][L1][L1+M1], coeffs_A[1][L1][L1+M1], coeffs_A[2][L1][L1+M1] };
			std::vector< std::complex<double> > Psi1 = VSH_Psi(t,p,L1,M1);
			std::vector< std::complex<double> > Phi1 = VSH_Phi(t,p,L1,M1);
			
			std::vector< std::complex<double> > B_r_Y_vector { 0, Ar12[2]*Psi1[1]-Ar12[1]*Phi1[1], Ar12[2]*Psi1[2]-Ar12[1]*Phi1[2] };
			
			for( int L2=0; L2<=ell_max; L2++ ){
				for( int M2=-L2; M2<=L2; M2++ ){
					
					std::vector< std::complex<double> > Br12 { coeffs_B[0][L2][L2+M2], coeffs_B[1][L2][L2+M2], coeffs_B[2][L2][L2+M2] };
					std::vector< std::complex<double> > Psi2 = VSH_Psi(t,p,L2,M2);
					std::vector< std::complex<double> > Phi2 = VSH_Phi(t,p,L2,M2);
					
					//--- Radial component ---
					std::complex<double> coeff_Psi = Ar12[1]*Br12[2] - Ar12[2]*Br12[1];
					std::complex<double> coeff_Phi = Ar12[1]*Br12[1] + Ar12[2]*Br12[2];
					std::vector< std::complex<double> > coeff { 0, coeff_Psi*Psi1[1]+coeff_Phi*Phi1[1], coeff_Psi*Psi1[2]+coeff_Phi*Phi1[2] };
					cross_product_from_VSH_coeffs_a[0] += coeff[1]*Psi2[1] + coeff[2]*Psi2[2];
					
					//--- Angular components ---
					std::vector< std::complex<double> > A_r_Y_vector { 0, Br12[1]*Phi2[1]-Br12[2]*Psi2[1], Br12[1]*Phi2[2]-Br12[2]*Psi2[2] };
					for( int i=1; i<=2; i++ ){
						cross_product_from_VSH_coeffs_a[i] += Br12[0]*ylm(t,p,L2,M2) * B_r_Y_vector[i] + Ar12[0]*ylm(t,p,L1,M1) * A_r_Y_vector[i];
					}
					
				}
			}
		}
	}*/
	
	/*//--- V3 (Don't use function for VSHs; turn B_r_Y_vector and A_r_Y_vector into the final expression now that we know it's correct) ---
	for( int L1=0; L1<=ell_max; L1++ ){
		for( int M1=-L1; M1<=L1; M1++ ){
			
			//--- Precalculate VSHs for first indices ---
			std::vector< std::complex<double> > Ar12 { coeffs_A[0][L1][L1+M1], coeffs_A[1][L1][L1+M1], coeffs_A[2][L1][L1+M1] };
			std::vector< std::complex<double> > Psi1 { 0, ylm_dtheta(t,p,L1,M1), 0 };
			std::vector< std::complex<double> > Phi1 { 0, 0, Psi1[1] };
			if( M1 != 0 ){
				Psi1[2] =  std::complex<double>{0,(double)M1} * ylm_over_sintheta(t,p,L1,M1);
				Phi1[1] = -Psi1[2];
			}
			
			for( int L2=0; L2<=ell_max; L2++ ){
				for( int M2=-L2; M2<=L2; M2++ ){
					
					//--- Precalculate VSHs for second indices ---
					std::vector< std::complex<double> > Br12 { coeffs_B[0][L2][L2+M2], coeffs_B[1][L2][L2+M2], coeffs_B[2][L2][L2+M2] };
					std::vector< std::complex<double> > Psi2 { 0, ylm_dtheta(t,p,L2,M2), 0 };
					std::vector< std::complex<double> > Phi2 { 0, 0, Psi2[1] };
					if( M2 != 0 ){
						Psi2[2] =  std::complex<double>{0,(double)M2} * ylm_over_sintheta(t,p,L2,M2);
						Phi2[1] = -Psi2[2];
					}
					
					//--- Radial component ---
					std::complex<double> coeff_Psi = Ar12[1]*Br12[2] - Ar12[2]*Br12[1];
					std::complex<double> coeff_Phi = Ar12[1]*Br12[1] + Ar12[2]*Br12[2];
					cross_product_from_VSH_coeffs_a[0] += ( coeff_Psi*Psi1[1]+coeff_Phi*Phi1[1] )*Psi2[1] + ( coeff_Psi*Psi1[2]+coeff_Phi*Phi1[2] )*Psi2[2];
					
					//--- Angular components ---
					for( int i=1; i<=2; i++ ){
						cross_product_from_VSH_coeffs_a[i] += Br12[0]*ylm(t,p,L2,M2) * ( Ar12[2]*Psi1[i]-Ar12[1]*Phi1[i] ) + Ar12[0]*ylm(t,p,L1,M1) * ( Br12[1]*Phi2[i]-Br12[2]*Psi2[i] );
						
					}
				}
			}
		}
	}*/
	
	
	/*//--- V4 (Explicitly calculate ell_2 = 0 terms separately) ---
	for( int L1=0; L1<=ell_max; L1++ ){
		for( int M1=-L1; M1<=L1; M1++ ){
			
			//--- Precalculate VSHs for first indices ---
			std::vector< std::complex<double> > Ar12 { coeffs_A[0][L1][L1+M1], coeffs_A[1][L1][L1+M1], coeffs_A[2][L1][L1+M1] };
			std::vector< std::complex<double> > Psi1 { 0, ylm_dtheta(t,p,L1,M1), 0 };
			std::vector< std::complex<double> > Phi1 { 0, 0, Psi1[1] };
			if( M1 != 0 ){
				Psi1[2] =  std::complex<double>{0,(double)M1} * ylm_over_sintheta(t,p,L1,M1);
				Phi1[1] = -Psi1[2];
			}
			
			//--- ell_2 = 0 ---
			for( int i=0; i<3; i++ ){
				cross_product_from_VSH_coeffs_a[i] += coeffs_B[0][0][0] * pow(4.0*pi,-0.5) * ( Ar12[2] * Psi1[i] - Ar12[1] * Phi1[i] );
			}
			
			// --- ell_2 > 0 ---
			for( int L2=1; L2<=ell_max; L2++ ){
				for( int M2=-L2; M2<=L2; M2++ ){
					
					//--- Precalculate VSHs for second indices ---
					std::vector< std::complex<double> > Br12 { coeffs_B[0][L2][L2+M2], coeffs_B[1][L2][L2+M2], coeffs_B[2][L2][L2+M2] };
					std::vector< std::complex<double> > Psi2 { 0, ylm_dtheta(t,p,L2,M2), 0 };
					std::vector< std::complex<double> > Phi2 { 0, 0, Psi2[1] };
					if( M2 != 0 ){
						Psi2[2] =  std::complex<double>{0,(double)M2} * ylm_over_sintheta(t,p,L2,M2);
						Phi2[1] = -Psi2[2];
					}
					
					//--- Radial component ---
					std::complex<double> coeff_Psi = Ar12[1]*Br12[2] - Ar12[2]*Br12[1];
					std::complex<double> coeff_Phi = Ar12[1]*Br12[1] + Ar12[2]*Br12[2];
					cross_product_from_VSH_coeffs_a[0] += ( coeff_Psi*Psi1[1]+coeff_Phi*Phi1[1] )*Psi2[1] + ( coeff_Psi*Psi1[2]+coeff_Phi*Phi1[2] )*Psi2[2];
					
					//--- Angular components ---
					for( int i=1; i<=2; i++ ){
						cross_product_from_VSH_coeffs_a[i] += Br12[0]*ylm(t,p,L2,M2) * ( Ar12[2]*Psi1[i]-Ar12[1]*Phi1[i] ) + Ar12[0]*ylm(t,p,L1,M1) * ( Br12[1]*Phi2[i]-Br12[2]*Psi2[i] );
					}
				}
			}
		}
	}*/
	
	/*//--- V5 (Explicitly calculate ell_1 = 0 and ell_2 = 0 terms separately) ---
	//--- ell_1 = 0 ---	
	for( int L2=1; L2<=ell_max; L2++ ){
		for( int M2=-L2; M2<=L2; M2++ ){
			
			//--- Precalculate VSHs for second indices ---
			std::vector< std::complex<double> > Br12 { coeffs_B[0][L2][L2+M2], coeffs_B[1][L2][L2+M2], coeffs_B[2][L2][L2+M2] };
			std::vector< std::complex<double> > Psi2 { 0, ylm_dtheta(t,p,L2,M2), 0 };
			std::vector< std::complex<double> > Phi2 { 0, 0, Psi2[1] };
			if( M2 != 0 ){
				Psi2[2] =  std::complex<double>{0,(double)M2} * ylm_over_sintheta(t,p,L2,M2);
				Phi2[1] = -Psi2[2];
			}
			
			//--- Angular components ---
			for( int i=1; i<=2; i++ ){
				cross_product_from_VSH_coeffs_a[i] += coeffs_A[0][0][0]*pow(4.0*pi,-0.5) * ( Br12[1]*Phi2[i]-Br12[2]*Psi2[i] );
			}
		}
	}
	
	//--- ell_1 > 0 ---
	for( int L1=1; L1<=ell_max; L1++ ){
		for( int M1=-L1; M1<=L1; M1++ ){
			
			//--- Precalculate VSHs for first indices ---
			std::vector< std::complex<double> > Ar12 { coeffs_A[0][L1][L1+M1], coeffs_A[1][L1][L1+M1], coeffs_A[2][L1][L1+M1] };
			std::vector< std::complex<double> > Psi1 { 0, ylm_dtheta(t,p,L1,M1), 0 };
			std::vector< std::complex<double> > Phi1 { 0, 0, Psi1[1] };
			if( M1 != 0 ){
				Psi1[2] =  std::complex<double>{0,(double)M1} * ylm_over_sintheta(t,p,L1,M1);
				Phi1[1] = -Psi1[2];
			}
			
			//--- ell_2 = 0 ---
			for( int i=0; i<3; i++ ){
				cross_product_from_VSH_coeffs_a[i] += coeffs_B[0][0][0] * pow(4.0*pi,-0.5) * ( Ar12[2] * Psi1[i] - Ar12[1] * Phi1[i] );
			}
			
			// --- ell_2 > 0 ---
			for( int L2=1; L2<=ell_max; L2++ ){
				for( int M2=-L2; M2<=L2; M2++ ){
					
					//--- Precalculate VSHs for second indices ---
					std::vector< std::complex<double> > Br12 { coeffs_B[0][L2][L2+M2], coeffs_B[1][L2][L2+M2], coeffs_B[2][L2][L2+M2] };
					std::vector< std::complex<double> > Psi2 { 0, ylm_dtheta(t,p,L2,M2), 0 };
					std::vector< std::complex<double> > Phi2 { 0, 0, Psi2[1] };
					if( M2 != 0 ){
						Psi2[2] =  std::complex<double>{0,(double)M2} * ylm_over_sintheta(t,p,L2,M2);
						Phi2[1] = -Psi2[2];
					}
					
					//--- Radial component ---
					std::complex<double> coeff_Psi = Ar12[1]*Br12[2] - Ar12[2]*Br12[1];
					std::complex<double> coeff_Phi = Ar12[1]*Br12[1] + Ar12[2]*Br12[2];
					cross_product_from_VSH_coeffs_a[0] += ( coeff_Psi*Psi1[1]+coeff_Phi*Phi1[1] )*Psi2[1] + ( coeff_Psi*Psi1[2]+coeff_Phi*Phi1[2] )*Psi2[2];
					
					//--- Angular components ---
					for( int i=1; i<=2; i++ ){
						cross_product_from_VSH_coeffs_a[i] += Br12[0]*ylm(t,p,L2,M2) * ( Ar12[2]*Psi1[i]-Ar12[1]*Phi1[i] ) + Ar12[0]*ylm(t,p,L1,M1) * ( Br12[1]*Phi2[i]-Br12[2]*Psi2[i] );
					}
				}
			}
		}
	}*/
	
	
	//--- V6 (As a VSH series, Book 11, p56) ---
	// Calculate components of B
	std::vector< std::complex<double> > A_test  = { 0, 0, 0 };
	std::vector< std::complex<double> > B_test  = { 0, 0, 0 };
	// r-component
	for( int ell=0; ell<=ell_max; ell++ ){
		for( int m=-ell; m<=ell; m++ ){
			A_test[0] += coeffs_A[0][ell][ell+m] * ylm(t,p,ell,m);
			B_test[0] += coeffs_B[0][ell][ell+m] * ylm(t,p,ell,m);
		}
	}
	
	// theta-component and phi-component
	for( int ell=1; ell<=ell_max; ell++ ){
		for( int m=-ell; m<=ell; m++ ){
			A_test[1] += coeffs_A[1][ell][ell+m] * ylm_dtheta(t,p,ell,m) - coeffs_A[2][ell][ell+m] * ylm_dphi(t,p,ell,m) * pow(sin(t),-1);
			A_test[2] += coeffs_A[2][ell][ell+m] * ylm_dtheta(t,p,ell,m) + coeffs_A[1][ell][ell+m] * ylm_dphi(t,p,ell,m) * pow(sin(t),-1);
			B_test[1] += coeffs_B[1][ell][ell+m] * ylm_dtheta(t,p,ell,m) - coeffs_B[2][ell][ell+m] * ylm_dphi(t,p,ell,m) * pow(sin(t),-1);
			B_test[2] += coeffs_B[2][ell][ell+m] * ylm_dtheta(t,p,ell,m) + coeffs_B[1][ell][ell+m] * ylm_dphi(t,p,ell,m) * pow(sin(t),-1);
		}
	}
	print_vector( A_test, "A at test point" );
	print_vector( B_test, "B at test point" );
	
	// Now calculate cross product
	for( int ell=0; ell<=ell_max; ell++ ){
		for( int m=-ell; m<=ell; m++ ){
			std::complex<double> Y_term   = ( coeffs_A[1][ell][ell+m] * B_test[2] - coeffs_A[2][ell][ell+m] * B_test[1] ) * pow(ylm(t,p,ell,m),-1) * ylm_dtheta(t,p,ell,m)
			                              - ( coeffs_A[1][ell][ell+m] * B_test[1] + coeffs_A[2][ell][ell+m] * B_test[2] ) * pow(sin(t),-1) * std::complex<double>{0,(double)m};
			std::complex<double> Psi_term =   coeffs_A[2][ell][ell+m] * B_test[0] - A_test[0] * coeffs_B[2][ell][ell+m];
			std::complex<double> Phi_term = - coeffs_A[1][ell][ell+m] * B_test[0] + A_test[0] * coeffs_B[1][ell][ell+m];
			
			cross_product_from_VSH_coeffs_a = vector_sum( cross_product_from_VSH_coeffs_a, scalar_times_vector(   Y_term, VSH_Y  (t,p,ell,m) ) );
			cross_product_from_VSH_coeffs_a = vector_sum( cross_product_from_VSH_coeffs_a, scalar_times_vector( Psi_term, VSH_Psi(t,p,ell,m) ) );
			cross_product_from_VSH_coeffs_a = vector_sum( cross_product_from_VSH_coeffs_a, scalar_times_vector( Phi_term, VSH_Phi(t,p,ell,m) ) );
		}
	}
				

				
	
	//----- Calculate exact expression -----
	for( int L1=0; L1<=ell_max; L1++ ){
		for( int M1=-L1; M1<=L1; M1++ ){
			for( int L2=0; L2<=ell_max; L2++ ){
				for( int M2=-L2; M2<=L2; M2++ ){
			
					//--- Precalculate vectors and VSHs at the test point for these indices ---
					std::vector< std::complex<double> > Ar12 { coeffs_A[0][L1][L1+M1], coeffs_A[1][L1][L1+M1], coeffs_A[2][L1][L1+M1] };
					std::vector< std::complex<double> > Br12 { coeffs_B[0][L2][L2+M2], coeffs_B[1][L2][L2+M2], coeffs_B[2][L2][L2+M2] };
					std::vector< std::complex<double> > Y1   = VSH_Y  (t,p,L1,M1);
					std::vector< std::complex<double> > Y2   = VSH_Y  (t,p,L2,M2);
					std::vector< std::complex<double> > Psi1 = VSH_Psi(t,p,L1,M1);
					std::vector< std::complex<double> > Psi2 = VSH_Psi(t,p,L2,M2);
					std::vector< std::complex<double> > Phi1 = VSH_Phi(t,p,L1,M1);
					std::vector< std::complex<double> > Phi2 = VSH_Phi(t,p,L2,M2);
					std::complex<double> SH1 = ylm(t,p,L1,M1);
					std::complex<double> SH2 = ylm(t,p,L2,M2);
					
					//--- Add to each vector in the sum ---
					cross_product_from_VSH_coeffs_01 = vector_sum( cross_product_from_VSH_coeffs_01, scalar_times_vector( Ar12[0] * Br12[1] * SH1, Phi2 ) );
					cross_product_from_VSH_coeffs_02 = vector_sum( cross_product_from_VSH_coeffs_02, scalar_times_vector( Ar12[0] * Br12[2] *-SH1, Psi2 ) );
					cross_product_from_VSH_coeffs_10 = vector_sum( cross_product_from_VSH_coeffs_10, scalar_times_vector( Ar12[1] * Br12[0] *-SH2, Phi1 ) );
					cross_product_from_VSH_coeffs_11 = vector_sum( cross_product_from_VSH_coeffs_11, scalar_times_vector( Ar12[1] * Br12[1], cross_product( Psi1, Psi2 ) ) );
					cross_product_from_VSH_coeffs_12 = vector_sum( cross_product_from_VSH_coeffs_12, scalar_times_vector( Ar12[1] * Br12[2], cross_product( Psi1, Phi2 ) ) );
					cross_product_from_VSH_coeffs_20 = vector_sum( cross_product_from_VSH_coeffs_20, scalar_times_vector( Ar12[2] * Br12[0] * SH2, Psi1 ) );
					cross_product_from_VSH_coeffs_21 = vector_sum( cross_product_from_VSH_coeffs_21, scalar_times_vector( Ar12[2] * Br12[1], cross_product( Phi1, Psi2 ) ) );
					cross_product_from_VSH_coeffs_22 = vector_sum( cross_product_from_VSH_coeffs_22, scalar_times_vector( Ar12[2] * Br12[2], cross_product( Phi1, Phi2 ) ) );
				}
			}
		}
	}
	
	cross_product_from_VSH_coeffs = vector_sum( cross_product_from_VSH_coeffs, cross_product_from_VSH_coeffs_01 );
	cross_product_from_VSH_coeffs = vector_sum( cross_product_from_VSH_coeffs, cross_product_from_VSH_coeffs_02 );
	cross_product_from_VSH_coeffs = vector_sum( cross_product_from_VSH_coeffs, cross_product_from_VSH_coeffs_10 );
	cross_product_from_VSH_coeffs = vector_sum( cross_product_from_VSH_coeffs, cross_product_from_VSH_coeffs_11 );
	cross_product_from_VSH_coeffs = vector_sum( cross_product_from_VSH_coeffs, cross_product_from_VSH_coeffs_12 );
	cross_product_from_VSH_coeffs = vector_sum( cross_product_from_VSH_coeffs, cross_product_from_VSH_coeffs_20 );
	cross_product_from_VSH_coeffs = vector_sum( cross_product_from_VSH_coeffs, cross_product_from_VSH_coeffs_21 );
	cross_product_from_VSH_coeffs = vector_sum( cross_product_from_VSH_coeffs, cross_product_from_VSH_coeffs_22 );
	
	print_vector( cross_product_from_VSH_coeffs, "Cross product from VSH coeffs" );
	
	std::cout << "\n\n" << std::endl;
	
	int w = 20;
	
	print_vector_setw( cross_product_from_VSH_coeffs  , w, "Exact" );
	print_vector_setw( cross_product_from_VSH_coeffs_a, w, "New  " );
	
	
	
	//----- Output execution time and finish code -----
	auto time_stop = std::chrono::high_resolution_clock::now();
	double exec_time = std::chrono::duration_cast<std::chrono::nanoseconds>( time_stop - time_start ).count() * 1e-9;
	std::cout << "\nExecution time:\t" << floor(exec_time/60) << " m " << exec_time-60*floor(exec_time/60) << " s." << std::endl;
	
	std::cout << "\nFinished" << std::endl;
}