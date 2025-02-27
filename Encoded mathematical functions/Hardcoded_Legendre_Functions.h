#include <math.h>

//----- Legendre polynomials -----
double p0 ( double x ){
	return 1;
}
double p1 ( double x ){
	return x;
}
double p2 ( double x ){
	return ( 3*x*x - 1 ) * 0.5;
}
double p3 ( double x ){
	return ( 5*x*x*x - 3*x ) * 0.5;
}
double p4 ( double x ){
	return ( 35*pow(x,4) - 30*pow(x,2) + 3 ) * 0.125;
}
double p5 ( double x ){
	return ( 63*pow(x,5) - 70*pow(x,3) + 15*x ) * 0.125;
}
double p6 ( double x ){
	return ( 231*pow(x,6) - 315*pow(x,4) + 105*x*x - 5 ) * 0.0625;
}
double p7 ( double x ){
	return ( 429*pow(x,7) - 693*pow(x,5) + 315*pow(x,3) - 35*x ) * 0.0625;
}
double p8 ( double x ){
	return ( 6435*pow(x,8) - 12012*pow(x,6) + 6930*pow(x,4) - 1260*x*x + 35 ) / 128;
}
double p9 ( double x ){
	return ( 12155*pow(x,9) - 25740*pow(x,7) + 18018*pow(x,5) - 4620*pow(x,3) + 315*x ) / 128;
}
double p10( double x ){
	return ( 46189*pow(x,10) - 109395*pow(x,8) + 90090*pow(x,6) - 30030*pow(x,4) + 3465*x*x - 63 ) / 256;
}


double p0_costheta ( double theta ){
	double x = cos(theta); return 1;
}
double p1_costheta ( double theta ){
	double x = cos(theta); return x;
}
double p2_costheta ( double theta ){
	double x = cos(theta); return ( 3*x*x - 1 ) * 0.5;
}
double p3_costheta ( double theta ){
	double x = cos(theta); return ( 5*x*x*x - 3*x ) * 0.5;
}
double p4_costheta ( double theta ){
	double x = cos(theta); return ( 35*pow(x,4) - 30*pow(x,2) + 3 ) * 0.125;
}
double p5_costheta ( double theta ){
	double x = cos(theta); return ( 63*pow(x,5) - 70*pow(x,3) + 15*x ) * 0.125;
}
double p6_costheta ( double theta ){
	double x = cos(theta); return ( 231*pow(x,6) - 315*pow(x,4) + 105*x*x - 5 ) * 0.0625;
}
double p7_costheta ( double theta ){
	double x = cos(theta); return ( 429*pow(x,7) - 693*pow(x,5) + 315*pow(x,3) - 35*x ) * 0.0625;
}
double p8_costheta ( double theta ){
	double x = cos(theta); return ( 6435*pow(x,8) - 12012*pow(x,6) + 6930*pow(x,4) - 1260*x*x + 35 ) / 128;
}
double p9_costheta ( double theta ){
	double x = cos(theta); return ( 12155*pow(x,9) - 25740*pow(x,7) + 18018*pow(x,5) - 4620*pow(x,3) + 315*x ) / 128;
}
double p10_costheta( double theta ){
	double x = cos(theta); return ( 46189*pow(x,10) - 109395*pow(x,8) + 90090*pow(x,6) - 30030*pow(x,4) + 3465*x*x - 63 ) / 256;
}




//----- Legendre polynomials: First derivatives (unverified)-----
double p0_dx ( double x ){
	return 0;
}
double p1_dx ( double x ){
	return 1;
}
double p2_dx ( double x ){
	return ( 6*x - 0 ) * 0.5;
}
double p3_dx ( double x ){
	return ( 3*5*x*x - 3 ) * 0.5;
}
double p4_dx ( double x ){
	return ( 4*35*pow(x,3) - 2*30*x + 0 ) * 0.125;
}
double p5_dx ( double x ){
	return ( 5*63*pow(x,4) - 3*70*pow(x,2) + 15 ) * 0.125;
}
double p6_dx ( double x ){
	return ( 6*231*pow(x,5) - 4*315*pow(x,3) + 2*105*x  ) * 0.0625;
}
double p7_dx ( double x ){
	return ( 7*429*pow(x,6) - 5*693*pow(x,4) + 3*315*pow(x,2) - 35 ) * 0.0625;
}
double p8_dx ( double x ){
	return ( 8*6435*pow(x,7) - 6*12012*pow(x,5) + 4*6930*pow(x,3) - 2*1260*x  ) / 128;
}
double p9_dx ( double x ){
	return ( 9*12155*pow(x,8) - 7*25740*pow(x,6) + 5*18018*pow(x,4) - 3*4620*pow(x,2) + 315 ) / 128;
}
double p10_dx( double x ){
	return ( 10*46189*pow(x,9) - 8*109395*pow(x,7) + 6*90090*pow(x,5) - 4*30030*pow(x,3) + 2*3465*x 	) / 256;
}





//----- Legendre polynomials: Second derivatives (unverified)-----
double p0_dx2 ( double x ){
	return 0;
}
double p1_dx2 ( double x ){
	return 0;
}
double p2_dx2 ( double x ){
	return ( 6 ) * 0.5;
}
double p3_dx2 ( double x ){
	return ( 2*3*5*x ) * 0.5;
}
double p4_dx2 ( double x ){
	return ( 3*4*35*pow(x,2) - 2*30 ) * 0.125;
}
double p5_dx2 ( double x ){
	return ( 4*5*63*pow(x,3) - 2*3*70*x ) * 0.125;
}
double p6_dx2 ( double x ){
	return ( 5*6*231*pow(x,4) - 3*4*315*pow(x,2) + 2*105 ) * 0.0625;
}
double p7_dx2 ( double x ){
	return ( 6*7*429*pow(x,5) - 4*5*693*pow(x,3) + 2*3*315*x ) * 0.0625;
}
double p8_dx2 ( double x ){
	return ( 7*8*6435*pow(x,6) - 5*6*12012*pow(x,4) + 3*4*6930*pow(x,2) - 2*1260 ) / 128;
}
double p9_dx2 ( double x ){
	return ( 8*9*12155*pow(x,7) - 6*7*25740*pow(x,5) + 4*5*18018*pow(x,3) - 2*3*4620*x ) / 128;
}
double p10_dx2( double x ){
	return ( 9*10*46189*pow(x,8) - 7*8*109395*pow(x,6) + 5*6*90090*pow(x,4) - 3*4*30030*pow(x,2) + 2*3465 ) / 256;
}








//----- Associated Legendre polynomials -----
double alf00 (double x ){
	return 1;
}

double alf11 ( double x ){
	return -      pow( 1-x*x, 0.5 );
}
double alf10 ( double x ){ return p1( x );          }
double alf1m1( double x ){ return alf11(x) / -2;    }


double alf22 ( double x ){
	return 3    *    ( 1-x*x      );
}
double alf21 ( double x ){
	return -3   * pow( 1-x*x, 0.5 ) * x;

}
double alf20 ( double x ){ return p2( x );          }
double alf2m1( double x ){ return alf21(x) / -6;    }
double alf2m2( double x ){ return alf22(x) / 24;    }


double alf33 ( double x ){
	return -15  * pow( 1-x*x, 1.5 );
}
double alf32 ( double x ){
	return 15   *    ( 1-x*x      ) * x;
}
double alf31 ( double x ){
	return 1.5  * pow( 1-x*x, 0.5 ) * ( 1 - 5*x*x );
}
double alf30 ( double x ){ return p3( x );          }
double alf3m1( double x ){ return alf31(x) / -12;   }
double alf3m2( double x ){ return alf32(x) / 120;   }
double alf3m3( double x ){ return alf33(x) / -720;  }


double alf44 ( double x ){
	return 105  * pow( 1-x*x, 2 );
}
double alf43 ( double x ){
	return -105 * pow( 1-x*x, 1.5 ) * x;
}
double alf42 ( double x ){
	return 7.5  *    ( 1-x*x      ) * ( 7*x*x      - 1   );
}
double alf41 ( double x ){	
	return -2.5 * pow( 1-x*x, 0.5 ) * ( 7*pow(x,3) - 3*x );
}
double alf40 ( double x ){ return p4( x );          }
double alf4m1( double x ){ return alf41(x) / -20;   }
double alf4m2( double x ){ return alf42(x) / 360;   }
double alf4m3( double x ){ return alf43(x) / -5040; }
double alf4m4( double x ){ return alf44(x) / 40320; }




//----- Associated Legendre polynomials: First derivatives -----
double alf00_dx (double x ){
	return 0;
}

double alf11_dx ( double x ){
	return pow( 1-x*x, -0.5 ) * x;
}
double alf10_dx ( double x ){ return p1_dx( x );          }
double alf1m1_dx( double x ){ return alf11_dx(x) / -2;    }


double alf22_dx ( double x ){
	return -6 * x;
}
double alf21_dx ( double x ){
	return pow( 1-x*x, -0.5 ) * (2*x*x - 1 ) * 3;

}
double alf20_dx ( double x ){ return p2_dx( x );          }
double alf2m1_dx( double x ){ return alf21_dx(x) / -6;    }
double alf2m2_dx( double x ){ return alf22_dx(x) / 24;    }


double alf33_dx ( double x ){
	return pow( 1-x*x,  0.5 ) * 45 * x;
}
double alf32_dx ( double x ){
	return 15 * ( 1 - 3*x*x );
}
double alf31_dx ( double x ){
	return pow( 1-x*x, -0.5 ) * ( 15*x*x - 11 ) * 1.5*x;
}
double alf30_dx ( double x ){ return p3_dx( x );          }
double alf3m1_dx( double x ){ return alf31_dx(x) / -12;   }
double alf3m2_dx( double x ){ return alf32_dx(x) / 120;   }
double alf3m3_dx( double x ){ return alf33_dx(x) / -720;  }


double alf44_dx ( double x ){
	return 420 * x * ( x*x - 1 );
}
double alf43_dx ( double x ){
	return pow( 1-x*x, 0.5 ) * ( 4*x*x - 1 ) * 105;
}
double alf42_dx ( double x ){
	return 30 * x * ( -7*x*x + 4 );
}
double alf41_dx ( double x ){	
	return pow( 1-x*x, -0.5 ) * ( 28*pow(x,4) - 27*x*x + 3 ) * 2.5;
}
double alf40_dx ( double x ){ return p4_dx( x );          }
double alf4m1_dx( double x ){ return alf41_dx(x) / -20;   }
double alf4m2_dx( double x ){ return alf42_dx(x) / 360;   }
double alf4m3_dx( double x ){ return alf43_dx(x) / -5040; }
double alf4m4_dx( double x ){ return alf44_dx(x) / 40320; }




//----- Associated Legendre polynomials: Second derivatives -----
double alf00_dx2 (double x ){
	return 0;
}

double alf11_dx2 ( double x ){
	return pow( 1-x*x, -1.5 );
}
double alf10_dx2 ( double x ){ return p1_dx2( x );          }
double alf1m1_dx2( double x ){ return alf11_dx2(x) / -2;    }


double alf22_dx2 ( double x ){
	return -6;
}
double alf21_dx2 ( double x ){
	return pow( 1-x*x, -1.5 ) * ( -2*x*x + 3 ) * 3*x;

}
double alf20_dx2 ( double x ){ return p2_dx2( x );          }
double alf2m1_dx2( double x ){ return alf21_dx2(x) / -6;    }
double alf2m2_dx2( double x ){ return alf22_dx2(x) / 24;    }


double alf33_dx2 ( double x ){
	return pow( 1-x*x, -0.5 ) * ( 1 - 2*x*x ) * 45;
}
double alf32_dx2 ( double x ){
	return -90 * x;
}
double alf31_dx2 ( double x ){
	return pow( 1-x*x, -1.5 ) * ( -30*pow(x,4) + 45*x*x - 11 ) * 1.5;
}
double alf30_dx2 ( double x ){ return p3_dx2( x );          }
double alf3m1_dx2( double x ){ return alf31_dx2(x) / -12;   }
double alf3m2_dx2( double x ){ return alf32_dx2(x) / 120;   }
double alf3m3_dx2( double x ){ return alf33_dx2(x) / -720;  }


double alf44_dx2 ( double x ){
	return 420 * ( 3*x*x - 1 );
}
double alf43_dx2 ( double x ){
	return pow( 1-x*x, -0.5 ) * ( -4*x*x + 3 ) * 315*x;
}
double alf42_dx2 ( double x ){
	return 30 * ( -21*x*x + 4 );
}
double alf41_dx2 ( double x ){	
	return pow( 1-x*x, -1.5 ) * ( -84*pow(x,4) + 139*pow(x,2) - 51 ) * 2.5*x;
}
double alf40_dx2 ( double x ){ return p4_dx2( x );          }
double alf4m1_dx2( double x ){ return alf41_dx2(x) / -20;   }
double alf4m2_dx2( double x ){ return alf42_dx2(x) / 360;   }
double alf4m3_dx2( double x ){ return alf43_dx2(x) / -5040; }
double alf4m4_dx2( double x ){ return alf44_dx2(x) / 40320; }




//----- Single function giving all ALFs -----
double hardcoded_alf( double x, int ell, int m ){
	if( ell == 0 ){
		return alf00( x );
	}
	if( ell == 1 ){
		if( m == -1 ){ return alf1m1( x ); }
		if( m ==  0 ){ return alf10 ( x ); }
		if( m ==  1 ){ return alf11 ( x ); }
	}
	if( ell == 2 ){
		if( m == -2 ){ return alf2m2( x ); }
		if( m == -1 ){ return alf2m1( x ); }
		if( m ==  0 ){ return alf20 ( x ); }
		if( m ==  1 ){ return alf21 ( x ); }
		if( m ==  2 ){ return alf22 ( x ); }
	}
	if( ell == 3 ){
		if( m == -3 ){ return alf3m3( x ); }
		if( m == -2 ){ return alf3m2( x ); }
		if( m == -1 ){ return alf3m1( x ); }
		if( m ==  0 ){ return alf30 ( x ); }
		if( m ==  1 ){ return alf31 ( x ); }
		if( m ==  2 ){ return alf32 ( x ); }
		if( m ==  3 ){ return alf33 ( x ); }
	}
	if( ell == 4 ){
		if( m == -4 ){ return alf4m4( x ); }
		if( m == -3 ){ return alf4m3( x ); }
		if( m == -2 ){ return alf4m2( x ); }
		if( m == -1 ){ return alf4m1( x ); }
		if( m ==  0 ){ return alf40 ( x ); }
		if( m ==  1 ){ return alf41 ( x ); }
		if( m ==  2 ){ return alf42 ( x ); }
		if( m ==  3 ){ return alf43 ( x ); }
		if( m ==  4 ){ return alf44 ( x ); }
	}
	return 123;
}


double hardcoded_lp( double x, int ell ){
	if( ell ==  0 ){ return p0 ( x ); }
	if( ell ==  1 ){ return p1 ( x ); }
	if( ell ==  2 ){ return p2 ( x ); }
	if( ell ==  3 ){ return p3 ( x ); }
	if( ell ==  4 ){ return p4 ( x ); }
	if( ell ==  5 ){ return p5 ( x ); }
	if( ell ==  6 ){ return p6 ( x ); }
	if( ell ==  7 ){ return p7 ( x ); }
	if( ell ==  8 ){ return p8 ( x ); }
	if( ell ==  9 ){ return p9 ( x ); }
	if( ell == 10 ){ return p10( x ); }
	return 123;
}