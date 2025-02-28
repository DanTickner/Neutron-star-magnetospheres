/*
OneDrive\PhD\Codes\20221108 Vector spherical harmonics\Headers_VSH\VSH.h

Generate the vector spherical harmonics Y_ell^m = Y_ell^m e_r, Psi_ell^m = r grad Y_ell^m, Phi_ell^m = e_r x ( r grad Y_ell^m ).
Combine codes VSH_Y_01_h, VSH_Psi_02.h, VSH_Phi_01.h to save importing three header files.
*/

std::vector< std::complex<double> > VSH_Y( double theta, double phi, int ell, int m ){
	return std::vector< std::complex<double> > { ylm(theta,phi,ell,m), 0, 0 };
}

std::vector< std::complex<double> > VSH_Psi( double theta, double phi, int ell, int m ){
	std::vector< std::complex<double> > Psi (3);
	Psi[1] = ylm_dtheta( theta, phi, ell, m );
	if( m != 0 ) {
		Psi[2] = std::complex<double>{0,(double)m} * ylm_over_sintheta( theta, phi, ell, m );
	}
	return Psi;
}

std::vector< std::complex<double> > VSH_Phi( double theta, double phi, int ell, int m ){
	std::vector< std::complex<double> > Phi (3);
	if( m != 0 ) {
		Phi[1] = -std::complex<double>{0,(double)m} * ylm_over_sintheta( theta, phi, ell, m );
	}
	Phi[2] = ylm_dtheta( theta, phi, ell, m );
	return Phi;
}