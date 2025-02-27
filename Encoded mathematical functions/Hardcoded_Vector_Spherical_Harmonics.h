#include <vector>
//#include "Hardcoded_Spherical_Harmonics.h"

//const double pi = 3.14159265358979323846


//----- Y_ell^m -----
std::vector< std::complex<double> > VSH_Y_00 ( double theta, double phi ){
	std::complex<double> A = sqrt( 1.0 / ( 4 * pi ) );
	return std::vector< std::complex<double> > { A, 0, 0 };
}

std::vector< std::complex<double> > VSH_Y_11 ( double theta, double phi ){
	std::complex<double> A = - sqrt( 3.0 / ( 8 * pi ) ) * exp(std::complex<double>{0,phi});
	return std::vector< std::complex<double> > { A*sin(theta), 0, 0 };
}
std::vector< std::complex<double> > VSH_Y_10 ( double theta, double phi ){
	std::complex<double> A = sqrt( 3.0 / ( 4 * pi ) );
	return std::vector< std::complex<double> > { A*cos(theta), 0, 0 };
}
std::vector< std::complex<double> > VSH_Y_1m1( double theta, double phi ){
	std::complex<double> A = sqrt( 3.0 / ( 8 * pi ) ) * exp(std::complex<double>{0,-phi});
	return std::vector< std::complex<double> > { A*sin(theta), 0, 0 };
}

std::vector< std::complex<double> > VSH_Y_22 ( double theta, double phi ){
	std::complex<double> A = sqrt( 15.0 / ( 8 * pi ) ) * sin(theta) * exp(std::complex<double>{0,2*phi});
	return std::vector< std::complex<double> > { A*0.5*sin(theta), 0, 0 };
}
std::vector< std::complex<double> > VSH_Y_21 ( double theta, double phi ){
	std::complex<double> A = - sqrt( 15.0 / ( 8 * pi ) ) * exp(std::complex<double>{0,phi});
	return std::vector< std::complex<double> > { A*0.5*sin(2*theta), 0, 0 };
}
std::vector< std::complex<double> > VSH_Y_20 ( double theta, double phi ){
	std::complex<double> A = 0.25 * sqrt( 5.0 / ( pi ) );
	return std::vector< std::complex<double> > { A*(3*pow(cos(theta),2)-1), 0, 0 };
}
std::vector< std::complex<double> > VSH_Y_2m1( double theta, double phi ){
	std::complex<double> A = sqrt( 15.0 / ( 8 * pi ) ) * exp(std::complex<double>{0,-phi});
	return std::vector< std::complex<double> > { A*0.5*sin(2*theta), 0, 0 };
}
std::vector< std::complex<double> > VSH_Y_2m2( double theta, double phi ){
	std::complex<double> A = sqrt( 15.0 / ( 8 * pi ) ) * sin(theta) * exp(std::complex<double>{0,-2*phi});
	return std::vector< std::complex<double> > { A*0.5*sin(theta), 0, 0 };
}

std::vector< std::complex<double> > VSH_Y_33 ( double theta, double phi ){
	std::complex<double> A = -0.125 * sqrt( 35.0 / pi ) * pow(sin(theta),2) * exp(std::complex<double>{0,3*phi});
	return std::vector< std::complex<double> > { A*sin(theta), 0, 0 };
}
std::vector< std::complex<double> > VSH_Y_32 ( double theta, double phi ){
	std::complex<double> A = 0.25 * sqrt( 105.0 / ( 2.0 * pi ) ) * sin(theta) * exp(std::complex<double>{0,2*phi});
	return std::vector< std::complex<double> > { A*0.5*sin(2*theta), 0, 0 };
}
std::vector< std::complex<double> > VSH_Y_31 ( double theta, double phi ){
	std::complex<double> A = -0.125 * sqrt( 21.0 / pi ) * exp(std::complex<double>{0,phi});
	return std::vector< std::complex<double> > { A*sin(theta)*(5*pow(cos(theta),2)-1), 0, 0 };
}
std::vector< std::complex<double> > VSH_Y_30 ( double theta, double phi ){
	std::complex<double> A = 0.25 * sqrt( 7.0 / ( pi ) );
	return std::vector< std::complex<double> > { A*(5*pow(cos(theta),3)-3*cos(theta)), 0, 0 };
}
std::vector< std::complex<double> > VSH_Y_3m1( double theta, double phi ){
	std::complex<double> A = 0.125 * sqrt( 21.0 / pi ) * exp(std::complex<double>{0,-phi});
	return std::vector< std::complex<double> > { A*sin(theta)*(5*pow(cos(theta),2)-1), 0, 0 };
}
std::vector< std::complex<double> > VSH_Y_3m2( double theta, double phi ){
	std::complex<double> A = 0.25 * sqrt( 105.0 / ( 2.0 * pi ) ) * sin(theta) * exp(std::complex<double>{0,-2*phi});
	return std::vector< std::complex<double> > { A*0.5*sin(2*theta), 0, 0 };
}
std::vector< std::complex<double> > VSH_Y_3m3( double theta, double phi ){
	std::complex<double> A = 0.125 * sqrt( 35.0 / pi ) * pow(sin(theta),2) * exp(std::complex<double>{0,-3*phi});
	return std::vector< std::complex<double> > { A*sin(theta), 0, 0 };
}







//----- Psi_ell^m -----
std::vector< std::complex<double> > VSH_Psi_00 ( double theta, double phi ){
	std::complex<double> A = sqrt( 1.0 / ( 4 * pi ) );
	return std::vector< std::complex<double> > { 0, 0, 0 };
}

std::vector< std::complex<double> > VSH_Psi_11 ( double theta, double phi ){
	std::complex<double> A = - sqrt( 3.0 / ( 8 * pi ) ) * exp(std::complex<double>{0,phi});
	return std::vector< std::complex<double> > { 0, A*cos(theta), A*std::complex<double>{0,1} };
}
std::vector< std::complex<double> > VSH_Psi_10 ( double theta, double phi ){
	std::complex<double> A = sqrt( 3.0 / ( 4 * pi ) );
	return std::vector< std::complex<double> > { 0, A*-sin(theta), 0 };
}
std::vector< std::complex<double> > VSH_Psi_1m1( double theta, double phi ){
	std::complex<double> A = sqrt( 3.0 / ( 8 * pi ) ) * exp(std::complex<double>{0,-phi});
	return std::vector< std::complex<double> > { 0, A*cos(theta), A*std::complex<double>{0,-1} };
}

std::vector< std::complex<double> > VSH_Psi_22 ( double theta, double phi ){
	std::complex<double> A = sqrt( 15.0 / ( 8 * pi ) ) * sin(theta) * exp(std::complex<double>{0,2*phi});
	return std::vector< std::complex<double> > { 0, A*cos(theta), A*std::complex<double>{0,1} };
}
std::vector< std::complex<double> > VSH_Psi_21 ( double theta, double phi ){
	std::complex<double> A = - sqrt( 15.0 / ( 8 * pi ) ) * exp(std::complex<double>{0,phi});
	return std::vector< std::complex<double> > { 0, A*cos(2*theta), A*std::complex<double>{0,1}*cos(theta) };
}
std::vector< std::complex<double> > VSH_Psi_20 ( double theta, double phi ){
	std::complex<double> A = 0.25 * sqrt( 5.0 / ( pi ) );
	return std::vector< std::complex<double> > { 0, A*-3.0*sin(2*theta), 0 };
}
std::vector< std::complex<double> > VSH_Psi_2m1( double theta, double phi ){
	std::complex<double> A = sqrt( 15.0 / ( 8 * pi ) ) * exp(std::complex<double>{0,-phi});
	return std::vector< std::complex<double> > { 0, A*cos(2*theta), A*std::complex<double>{0,-1}*cos(theta) };
}
std::vector< std::complex<double> > VSH_Psi_2m2( double theta, double phi ){
	std::complex<double> A = sqrt( 15.0 / ( 8 * pi ) ) * sin(theta) * exp(std::complex<double>{0,-2*phi});
	return std::vector< std::complex<double> > { 0, A*cos(theta), A*std::complex<double>{0,-1} };
}

std::vector< std::complex<double> > VSH_Psi_33 ( double theta, double phi ){
	std::complex<double> A = -0.125 * sqrt( 35.0 / pi ) * pow(sin(theta),2) * exp(std::complex<double>{0,3*phi});
	return std::vector< std::complex<double> > { 0, A*3.0*cos(theta), A*3.0*std::complex<double>{0,1} };
}
std::vector< std::complex<double> > VSH_Psi_32 ( double theta, double phi ){
	std::complex<double> A = 0.25 * sqrt( 105.0 / ( 2.0 * pi ) ) * sin(theta) * exp(std::complex<double>{0,2*phi});
	return std::vector< std::complex<double> > { 0, A*(3.0*pow(cos(theta),2)-1), A*std::complex<double>{0,1}*2.0*cos(theta) };
}
std::vector< std::complex<double> > VSH_Psi_31 ( double theta, double phi ){
	std::complex<double> A = -0.125 * sqrt( 21.0 / pi ) * exp(std::complex<double>{0,phi});
	return std::vector< std::complex<double> > { 0, A*cos(theta)*(15*pow(cos(theta),2)-11), A*std::complex<double>{0,1}*(5*pow(cos(theta),2)-1) };
}
std::vector< std::complex<double> > VSH_Psi_30 ( double theta, double phi ){
	std::complex<double> A = 0.25 * sqrt( 7.0 / ( pi ) );
	return std::vector< std::complex<double> > { 0, A*3.0*sin(theta)*(1-5*pow(cos(theta),2)), 0 };
}
std::vector< std::complex<double> > VSH_Psi_3m1( double theta, double phi ){
	std::complex<double> A = 0.125 * sqrt( 21.0 / pi ) * exp(std::complex<double>{0,-phi});
	return std::vector< std::complex<double> > { 0, A*cos(theta)*(15*pow(cos(theta),2)-11), A*std::complex<double>{0,-1}*(5*pow(cos(theta),2)-1) };
}
std::vector< std::complex<double> > VSH_Psi_3m2( double theta, double phi ){
	std::complex<double> A = 0.25 * sqrt( 105.0 / ( 2.0 * pi ) ) * sin(theta) * exp(std::complex<double>{0,-2*phi});
	return std::vector< std::complex<double> > { 0, A*(3.0*pow(cos(theta),2)-1), A*std::complex<double>{0,-1}*2.0*cos(theta) };
}
std::vector< std::complex<double> > VSH_Psi_3m3( double theta, double phi ){
	std::complex<double> A = 0.125 * sqrt( 35.0 / pi ) * pow(sin(theta),2) * exp(std::complex<double>{0,-3*phi});
	return std::vector< std::complex<double> > { 0, A*3.0*cos(theta), A*3.0*std::complex<double>{0,-1} };
}




//----- Phi_ell^m -----

std::vector< std::complex<double> > VSH_Phi_00 ( double theta, double phi ){
	std::complex<double> A = sqrt( 1.0 / ( 4 * pi ) );
	return std::vector< std::complex<double> > { 0, 0, 0 };
}

std::vector< std::complex<double> > VSH_Phi_11 ( double theta, double phi ){
	std::complex<double> A = - sqrt( 3.0 / ( 8 * pi ) ) * exp(std::complex<double>{0,phi});
	return std::vector< std::complex<double> > { 0, A*-std::complex<double>{0,1}, A*cos(theta) };
}
std::vector< std::complex<double> > VSH_Phi_10 ( double theta, double phi ){
	std::complex<double> A = sqrt( 3.0 / ( 4 * pi ) );
	return std::vector< std::complex<double> > { 0, 0, A*-sin(theta) };
}
std::vector< std::complex<double> > VSH_Phi_1m1( double theta, double phi ){
	std::complex<double> A = sqrt( 3.0 / ( 8 * pi ) ) * exp(std::complex<double>{0,-phi});
	return std::vector< std::complex<double> > { 0, -A*std::complex<double>{0,-1}, A*cos(theta),  };
}

std::vector< std::complex<double> > VSH_Phi_22 ( double theta, double phi ){
	std::complex<double> A = sqrt( 15.0 / ( 8 * pi ) ) * sin(theta) * exp(std::complex<double>{0,2*phi});
	return std::vector< std::complex<double> > { 0, -A*std::complex<double>{0,1}, A*cos(theta) };
}
std::vector< std::complex<double> > VSH_Phi_21 ( double theta, double phi ){
	std::complex<double> A = - sqrt( 15.0 / ( 8 * pi ) ) * exp(std::complex<double>{0,phi});
	return std::vector< std::complex<double> > { 0, -A*std::complex<double>{0,1}*cos(theta), A*cos(2*theta) };
}
std::vector< std::complex<double> > VSH_Phi_20 ( double theta, double phi ){
	std::complex<double> A = 0.25 * sqrt( 5.0 / ( pi ) );
	return std::vector< std::complex<double> > { 0, 0, A*-3.0*sin(2*theta) };
}
std::vector< std::complex<double> > VSH_Phi_2m1( double theta, double phi ){
	std::complex<double> A = sqrt( 15.0 / ( 8 * pi ) ) * exp(std::complex<double>{0,-phi});
	return std::vector< std::complex<double> > { 0, -A*std::complex<double>{0,-1}*cos(theta), A*cos(2*theta) };
}
std::vector< std::complex<double> > VSH_Phi_2m2( double theta, double phi ){
	std::complex<double> A = sqrt( 15.0 / ( 8 * pi ) ) * sin(theta) * exp(std::complex<double>{0,-2*phi});
	return std::vector< std::complex<double> > { 0, -A*std::complex<double>{0,-1}, A*cos(theta) };
}

std::vector< std::complex<double> > VSH_Phi_33 ( double theta, double phi ){
	std::complex<double> A = -0.125 * sqrt( 35.0 / pi ) * pow(sin(theta),2) * exp(std::complex<double>{0,3*phi});
	return std::vector< std::complex<double> > { 0, -A*3.0*std::complex<double>{0,1}, A*3.0*cos(theta) };
}
std::vector< std::complex<double> > VSH_Phi_32 ( double theta, double phi ){
	std::complex<double> A = 0.25 * sqrt( 105.0 / ( 2.0 * pi ) ) * sin(theta) * exp(std::complex<double>{0,2*phi});
	return std::vector< std::complex<double> > { 0, A*std::complex<double>{0,1}*2.0*cos(theta), A*(3.0*pow(cos(theta),2)-1) };
}
std::vector< std::complex<double> > VSH_Phi_31 ( double theta, double phi ){
	std::complex<double> A = -0.125 * sqrt( 21.0 / pi ) * exp(std::complex<double>{0,phi});
	return std::vector< std::complex<double> > { 0, -A*std::complex<double>{0,1}*(5*pow(cos(theta),2)-1), A*cos(theta)*(15*pow(cos(theta),2)-11) };
}
std::vector< std::complex<double> > VSH_Phi_30 ( double theta, double phi ){
	std::complex<double> A = 0.25 * sqrt( 7.0 / ( pi ) );
	return std::vector< std::complex<double> > { 0, 0, A*3.0*sin(theta)*(1-5*pow(cos(theta),2)) };
}
std::vector< std::complex<double> > VSH_Phi_3m1( double theta, double phi ){
	std::complex<double> A = 0.125 * sqrt( 21.0 / pi ) * exp(std::complex<double>{0,-phi});
	return std::vector< std::complex<double> > { 0, -A*std::complex<double>{0,-1}*(5*pow(cos(theta),2)-1), A*cos(theta)*(15*pow(cos(theta),2)-11) };
}
std::vector< std::complex<double> > VSH_Phi_3m2( double theta, double phi ){
	std::complex<double> A = 0.25 * sqrt( 105.0 / ( 2.0 * pi ) ) * sin(theta) * exp(std::complex<double>{0,-2*phi});
	return std::vector< std::complex<double> > { 0, -A*std::complex<double>{0,-1}*2.0*cos(theta), A*(3.0*pow(cos(theta),2)-1) };
}
std::vector< std::complex<double> > VSH_Phi_3m3( double theta, double phi ){
	std::complex<double> A = 0.125 * sqrt( 35.0 / pi ) * pow(sin(theta),2) * exp(std::complex<double>{0,-3*phi});
	return std::vector< std::complex<double> > { 0, -A*3.0*std::complex<double>{0,-1}, A*3.0*cos(theta) };
}