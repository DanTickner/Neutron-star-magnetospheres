/*
Functions discarded but may prove useful for future reference.
*/

//----- Discarded 20240405 -----
// The first derivative was improved by this, but the second derivative appears identical to the result from finite differencing.

std::vector<double> parabola_from_three_points( double x0, double y0, double x1, double y1, double x2, double y2 ){
	// For any three pairs of gridpoints, there exists a unique parabola ax^2+bx+c passing through them.
	// Return the coefficients a,b,c.
	
	double det_A = x0*x0 * ( x1 - x2 ) + x0 * ( x2*x2 - x1*x1 ) + x1 * x2 * ( x1 - x2 );

	double a = ( ( x1 - x2 ) * y0 + ( x2 - x0 ) * y1 + ( x0 - x1 ) * y2 ) / det_A;
	double b = ( ( x2*x2 - x1*x1 ) * y0 + ( x0*x0 - x2*x2 ) * y1 + ( x1*x1 - x0*x0 ) * y2 ) / det_A;
	double c = ( x1 * x2 * ( x1 - x2 ) * y0 + x0 * x2 * ( x2 - x0 ) * y1 + x0 * x1 * ( x0 - x1 ) * y2 ) / det_A;
	
	return std::vector<double> { a, b, c };
}

std::vector<double> radial_derivatives_1_FD_enhanced_endpoints( std::vector<double> &v ){
	// Calculate the first radial derivatives of a vector v representing a function evaluated at the radial gridpoints r[i] used in this code.
	// Not intended for derivatives with respect to an arbitrary coordinate.
	// Construct parabolas at inner and outer boundaries; symmetric derivatives for intermediate points.
	
	std::vector<double> v_dr ( n_r );
	
	std::vector<double> parabola_inner = parabola_from_three_points( r[0], v[0], r[1], v[1], r[2], v[2] );
	v_dr[0] = 2.0 * parabola_inner[0] * r[0] + parabola_inner[1];
	
	for( int i=1; i<n_r-1; i++ ){
		v_dr[i] = ( v[i+1] - v[i-1] ) / ( 2.0 * delta_r );
	}
	
	std::vector<double> parabola_outer = parabola_from_three_points( r[n_r-3], v[n_r-3], r[n_r-2], v[n_r-2], r.back(), v.back() );
	v_dr.back() = 2.0 * parabola_outer[0] * r.back() + parabola_outer[1];
	
	return v_dr;
}

std::vector<double> radial_derivatives_2_FD_enhanced_endpoints( std::vector<double> &v ){
	// Calculate the second radial derivatives of a vector v representing a function evaluated at the radial gridpoints r[i] used in this code.
	// Not intended for derivatives with respect to an arbitrary coordinate.
	// Construct parabolas at inner and outer boundaries; symmetric derivatives for intermediate points.
	
	std::vector<double> v_dr2 ( n_r );
	
	//std::cout << v_dr2[0] << std::endl;
	//v_dr2[0] = ( v[2] - 2.0*v[1] + v[0] ) / ( delta_r*delta_r );
	//std::cout << v_dr2[0] << std::endl;
	//v_dr2[0] = 0;
	std::vector<double> parabola_inner = parabola_from_three_points( r[0], v[0], r[1], v[1], r[2], v[2] );
	v_dr2[0] = 2.0 * parabola_inner[0];
	//std::cout << v_dr2[0] << std::endl;
	
	for( int i=1; i<n_r-1; i++ ){
		v_dr2[i] = ( v[i+1] - 2.0*v[i] + v[i-1] ) / ( delta_r*delta_r );
	}
	
	//std::cout << v_dr2.back() << std::endl;
	//v_dr2.back() = ( v.back() - 2.0*v[n_r-2] + v[n_r-3] ) / ( delta_r * delta_r );
	//std::cout << v_dr2.back() << std::endl;
	//v_dr2.back() = 0;
	std::vector<double> parabola_outer = parabola_from_three_points( r[n_r-3], v[n_r-3], r[n_r-2], v[n_r-2], r.back(), v.back() );
	v_dr2.back() = 2.0 * parabola_outer[0];
	//std::cout << v_dr2.back() << std::endl;
	
	return v_dr2;
}