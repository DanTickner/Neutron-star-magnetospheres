// VSH_series_evaluation_08.h
// Evaluate the VSH series of a vector A with no r-dependence at a point (theta,phi), given its series coefficients.
// Precalculate the SHs.

std::vector< std::complex<double> > evaluate_vsh_series( double t, double p, std::vector< std::vector< std::vector< std::complex<double> > > > coeffs_VSH_series ){
	// Generate the C^{r,ell}_m coefficients in the vector spherical harmonic series of a vector function A(theta,phi) with no r-dependence.
	
	int ell_max = coeffs_VSH_series[0].size() - 1;
	std::vector< std::complex<double> > ret { 0, 0, 0 };
	
	for( int ell=0; ell<=ell_max; ell++ ){
		for( int m=-ell; m<=ell; m++ ){
			std::complex<double> SH = ylm(t,p,ell,m);
			std::complex<double> SH_dtheta = ylm_dtheta(t,p,ell,m);
			std::vector< std::complex<double> > Y   { SH, 0, 0 };
			std::vector< std::complex<double> > Psi { 0, SH_dtheta, std::complex<double>{0,m/sin(t)}*SH };
			std::vector< std::complex<double> > Phi { 0, std::complex<double>{0,-m/sin(t)}*SH, SH_dtheta };
			ret[0] += coeffs_VSH_series[0][ell][ell+m] * Y[0];
			ret[1] += coeffs_VSH_series[1][ell][ell+m] * Psi[1] + coeffs_VSH_series[2][ell][ell+m] * Phi[1];
			ret[2] += coeffs_VSH_series[1][ell][ell+m] * Psi[2] + coeffs_VSH_series[2][ell][ell+m] * Phi[2];
		}
	}
	
	return ret;
}