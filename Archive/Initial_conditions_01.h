/*
Initial_conditions_01.h

Specify the magnetic field strength, rotation rate, extent of radial coordinates and maximum basis function indices.
Doing this guarantees that values are in agreement between all files, and prevents lengthy-in-lines-of-code requirement to ensure maximum indices are equal.
*/

const double pi = acos( -1 );


//----- Maximum indices -----
int n_max   = 20;	// Maximum order of Chebyshev polynomials.
int ell_max = 20;	// Maximum order of Legendre  polynomials.


//----- Calculated properties of the star -----
double Omega = 0.2;
double R_LC  = 1.0 / Omega;							// Radius of the light cylinder in code units.


//----- Extent of radial coordinates (in code units) -----
double r_min    = 1.0;
double r_max    = 2.0 * R_LC;	// 0.9 * R_LC;
double r_sponge = 0.99 * r_max;	// Radius beyond which we wish to add a sponge layer (Parfrey, 2012, PhD thesis, S3.10).
double rho      = ( r_max + r_min ) / ( r_max - r_min );
double r_factor = 2.0 / ( r_max - r_min );				// For rescaling Chebyshev series from [-1,1] to [r_min,r_max].


//----- Extent of polar coordinates -----
double t_min = 0;
double t_max = pi;


//----- Properties of the star -----
double R_star = 1.0e4;			// Radius of the star in metres.
double m_star = 1.0e32;			// Magnetic dipole moment of the star in A m^2 (10^3 emu).


//----- Constants -----
double c    = 299792458;			// Speed of light in m s^-1.
double mu_0 = 1.25663706212e-6;	// Permeability of free space in kg m s^-2 A^-2.


//----- Factors to convert between SI units and code units -----
double T_factor = R_star / c;
double B_factor = mu_0/(4.0*pi) * m_star * pow( R_star, -3 );
double E_factor = B_factor / c;




//----- Define functions common to all files -----
double Lambda( double r ){
	// Map r in [-1,1] to [r_min,r_max].
	return 0.5 * ( (r_max-r_min)*r + r_min+r_max );
}

double Lambda_inverse( double r ){
	// Map r in [r_min,r_max] to [-1,1].
	return ( 2.0 * r - r_min - r_max ) / ( r_max - r_min );
}