// VSH_series_evaluation_10.h
// Evaluate the VSH series of a vector A with no r-dependence at a point (theta,phi), given its series coefficients.
// Separate the m=0 terms.

std::vector< std::complex<double> > evaluate_vsh_series( double t, double p, std::vector< std::vector< std::vector< std::complex<double> > > > coeffs_VSH_series ){
	// Generate the C^{r,ell}_m coefficients in the vector spherical harmonic series of a vector function A(theta,phi) with no r-dependence.
	
	int ell_max = coeffs_VSH_series[0].size() - 1;
	std::vector< std::complex<double> > ret { 0, 0, 0 };
	
	for( int ell=0; ell<=ell_max; ell++ ){
		
		//----- m=0 VSHs -----
		std::complex<double> SH = ylm(t,p,ell,0);
		std::complex<double> SH_dtheta = ylm_dtheta(t,p,ell,0);
		ret[0] += coeffs_VSH_series[0][ell][ell] * SH;
		ret[1] += coeffs_VSH_series[1][ell][ell] * SH_dtheta;
		ret[2] += coeffs_VSH_series[2][ell][ell] * SH_dtheta;
		
		//----- m!= VSHs -----
		for( int m=-ell; m<=ell; m++ ){
			if( m==0 ){ continue; }
			std::complex<double> SH = ylm(t,p,ell,m);
			std::complex<double> SH_dtheta = ylm_dtheta(t,p,ell,m);
			ret[0] += coeffs_VSH_series[0][ell][ell+m] * SH;
			ret[1] += coeffs_VSH_series[1][ell][ell+m] * SH_dtheta + coeffs_VSH_series[2][ell][ell+m] * std::complex<double>{0,-m/sin(t)}*SH;
			ret[2] += coeffs_VSH_series[1][ell][ell+m] * std::complex<double>{0,m/sin(t)}*SH + coeffs_VSH_series[2][ell][ell+m] * SH_dtheta;
		}
	}
	
	return ret;
}