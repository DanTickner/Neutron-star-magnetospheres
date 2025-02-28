// Vector_function_02.h
// Exact expressions for the three components of a spherical vector independent of r and its spatial derivatives, and the exact VSH series coefficients.
// A = [ cos(theta) cos(phi) ] e_r + [ cos(theta) sin(phi) ] e_theta + [ P_2(cos(theta)) cos(phi) ] e_phi, where P_2 is the second Legendre polynomial.

double A_r_function( double t, double p ){
	return cos(t) * cos(p);
}




double A_t_function( double t, double p ){
	return cos(t) * sin(p);
}




double A_p_function( double t, double p ){
	return 0.5*(3*pow(cos(t),2)-1) * cos(p);
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
	//return 8*r - 2.0/r - sin(t) + pow(cos(t),2)/sin(t);
	return 1;
}




std::vector<double> curl_A_function( double t, double p ){
	std::vector<double> ret (3);
	ret[0] = -1.5*sin(2*t) + 0.5*cos(t)/sin(t)*(3*pow(cos(t),2)-1);
	ret[1] = -3*pow(cos(t),2) + 1;
	ret[2] = 2*cos(t);
	return ret;
}