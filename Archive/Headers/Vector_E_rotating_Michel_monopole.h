/* 
Vector_E_rotating_Michel_monopole.h

Parfrey PhD thesis section 4.2.1
Michel F. C., 1973, ApJL, 180, https://ui.adsabs.harvard.edu/abs/1973ApJ...180L.133M/abstract

B = f_0/r^2 e_r with f_0 an arbitrary scalar.
E = - ( Omega x r ) x B

f_0 was defined in the header file for B.
*/

double E_r_function( double r, double t ){
	return 0;
}

double E_t_function( double r, double t ){
	return - f_0 * Omega * sin(t) * pow( r, -1 );
}

double E_p_function( double r, double t ){
	return 0;
}




double E_r_L_function( double r, int L ){
	// E^{r,L}(r).
	return 0;
}

double E_1_L_function( double r, int L ){
	// E^{(1),L}(r).
	if( L == 1 ){
		return 2.0 * sqrt( pi / 3.0 ) * f_0 * Omega * pow( r, -1 );
	}
	return 0;
}

double E_2_L_function( double r, int L ){
	// E^{*(2),L}(r).
	return 0;
}



std::vector<double> E_check( double r, double t ){
	// Calculate E by crossing Omega, r and B, to check the given expression for E.
	std::vector<double> Omega_vector { Omega * cos(t), -Omega * sin(t), 0 };
	std::vector<double> r_vector { r, 0, 0 };
	std::vector<double> B_vector { B_r_function(r,t), B_t_function(r,t), B_p_function(r,t) };
	return scalar_times_vector( -1, cross_product( cross_product( Omega_vector, r_vector ), B_vector ) );
}
	