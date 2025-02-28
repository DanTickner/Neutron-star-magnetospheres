// VSH_Y_01.h

// Generate the vector spherical harmonic Y_ell^m = Yscalar_ell^m(theta,phi) e_r.

std::vector< std::complex<double> > VSH_Y( double theta, double phi, int ell, int m ){
	return std::vector< std::complex<double> > { ylm(theta,phi,ell,m), 0, 0 };
}