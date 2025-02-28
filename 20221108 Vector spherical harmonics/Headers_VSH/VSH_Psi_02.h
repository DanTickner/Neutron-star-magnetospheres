// VSH_Psi_02.h

// Generate the vector spherical harmonic Psi_ell^m = r grad Y_ell^m.
// There will be an issue with m=1 when theta=0,pi.

std::vector< std::complex<double> > VSH_Psi( double theta, double phi, int ell, int m ){
	std::vector< std::complex<double> > Psi (3);
	Psi[1] = ylm_dtheta( theta, phi, ell, m );
	if( m != 0 ) {
		Psi[2] = std::complex<double>{0,(double)m} * ylm_over_sintheta( theta, phi, ell, m );
	}
	return Psi;
}