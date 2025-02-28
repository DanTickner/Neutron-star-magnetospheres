// Vector_function_01.h
// Exact expressions for the three components of a spherical vector independent of r and its spatial derivatives, and the exact VSH series coefficients.
// A = 1 e_r + 2 e_theta + 1 e_phi.

double A_r_function( double t, double p ){
	return 1;
}




double A_t_function( double t, double p ){
	return 1;
}




double A_p_function( double t, double p ){
	return 1;
}




/*std::vector< std::vector< std::vector<double> > > generate_exact_coeffs( int n_max, int ell_max ){
	
	//--- Radial    component ---
	std::vector< std::vector<double> > coeffs_r ( n_max+1, std::vector<double> ( ell_max+1 ) );
	coeffs_r[2][0] = 1;
	
	//--- Polar     component ---
	std::vector< std::vector<double> > coeffs_t ( n_max+1, std::vector<double> ( ell_max+1 ) );
	coeffs_t[1][1] = 1;
	
	//--- Azimuthal component ---
	std::vector< std::vector<double> > coeffs_p ( n_max+1, std::vector<double> ( ell_max+1 ) );
	coeffs_p[1][2] = 1;
	
	return std::vector< std::vector< std::vector<double> > > { coeffs_r, coeffs_t, coeffs_p };
}*/




double div_A_function( double t, double p ){
	return 0;
}




std::vector<double> curl_A_function( double t, double p ){
	return std::vector<double> { 0, 0, 0 };
}