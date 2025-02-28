// VSH_Psi_01.h

// Generate the vector spherical harmonic Y_ell^m = Yscalar_ell^m(theta,phi) e_r.
// There will be an issue with theta=0 and theta=pi. Partially remedied in the next version.

std::vector< std::complex<double> > VSH_Psi( double theta, double phi, int ell, int m ){
	std::vector< std::complex<double> > Psi (3);
	Psi[1] = ylm_dtheta( theta, phi, ell, m );
	Psi[2] = std::complex<double>{0,(double)m}/sin(theta) * ylm( theta, phi, ell, m );
	return Psi;
}