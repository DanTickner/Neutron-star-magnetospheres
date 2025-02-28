// VSH_Phi_01.h

// Generate the vector spherical harmonic Phi_ell^m = e_r x ( r grad Y_ell^m ).
// There will be an issue with m=1 when theta=0,pi.

std::vector< std::complex<double> > VSH_Phi( double theta, double phi, int ell, int m ){
	std::vector< std::complex<double> > Phi (3);
	if( m != 0 ) {
		Phi[1] = -std::complex<double>{0,(double)m} * ylm_over_sintheta( theta, phi, ell, m );
	}
	Phi[2] = ylm_dtheta( theta, phi, ell, m );
	return Phi;
}