// Hardcoded_Spherical_Harmonics.h

#include <math.h>
#include <complex>

//const double pi = 3.14159265358979323846;	// If using this header file in tandem with Fit_3D_Function.h, comment-out this line to prevent a redefinition error.

std::complex<double> y00 ( double theta, double phi ){
	return 0.5/sqrt(pi);
}


std::complex<double> y11 ( double theta, double phi ){
	return -0.5*sqrt(1.5/pi) * exp(std::complex<double>{0,phi}) * sin(theta);
}
std::complex<double> y10 ( double theta, double phi ){
	return 0.5 *sqrt(3.0/pi)                                    * cos(theta);
}
std::complex<double> y1m1( double theta, double phi ){
	return - std::conj( y11(theta,phi) );
}


std::complex<double> y22 ( double theta, double phi ){
	return 0.25 * sqrt(7.5/pi) * exp(std::complex<double>{0,2*phi}) * pow(sin(theta),2);
}
std::complex<double> y21 ( double theta, double phi ){
	return -0.5 * sqrt(7.5/pi) * exp(std::complex<double>{0,  phi}) *     sin(theta)    *         cos(theta);
}
std::complex<double> y20 ( double theta, double phi ){
	return 0.25 * sqrt(5.0/pi)                                                          * ( 3*pow(cos(theta),2) - 1 );
}
std::complex<double> y2m1( double theta, double phi ){
	return - std::conj( y21(theta,phi) );
}
std::complex<double> y2m2( double theta, double phi ){
	return   std::conj( y22(theta,phi) );
}


std::complex<double> y33 ( double theta, double phi ){
	return -0.125 * sqrt(35.0/pi) * exp(std::complex<double>{0,3*phi}) * pow(sin(theta),3);
}
std::complex<double> y32 ( double theta, double phi ){
	return 0.25   * sqrt(52.5/pi) * exp(std::complex<double>{0,2*phi}) * pow(sin(theta),2) *         cos(theta);
}
std::complex<double> y31 ( double theta, double phi ){
	return -0.125 * sqrt(21.0/pi) * exp(std::complex<double>{0,  phi}) *     sin(theta)    * ( 5*pow(cos(theta),2) - 1 );
}
std::complex<double> y30 ( double theta, double phi ){
	return 0.25   * sqrt(7.0 /pi)                                                          * ( 5*pow(cos(theta),3) - 3*cos(theta) );
}
std::complex<double> y3m1( double theta, double phi ){
	return - std::conj( y31(theta,phi) );
}
std::complex<double> y3m2( double theta, double phi ){
	return   std::conj( y32(theta,phi) );
}
std::complex<double> y3m3( double theta, double phi ){
	return - std::conj( y33(theta,phi) );
}


std::complex<double> y44 ( double theta, double phi ){
	return 0.1875 * sqrt(17.5/pi) * exp(std::complex<double>{0,4*phi}) * pow(sin(theta),4);
}
std::complex<double> y43 ( double theta, double phi ){
	return -0.375 * sqrt(35.0/pi) * exp(std::complex<double>{0,3*phi}) * pow(sin(theta),3) *          cos(theta);
}
std::complex<double> y42 ( double theta, double phi ){
	return 0.375  * sqrt(2.5 /pi) * exp(std::complex<double>{0,2*phi}) * pow(sin(theta),2) * ( 7 *pow(cos(theta),2) - 1 );
}
std::complex<double> y41 ( double theta, double phi ){
	return -0.375 * sqrt(5.0 /pi) * exp(std::complex<double>{0,  phi}) *     sin(theta)    * ( 7 *pow(cos(theta),3) - 3*cos(theta) );
}
std::complex<double> y40 ( double theta, double phi ){
	return 0.1875 / sqrt(     pi)                                                          * ( 35*pow(cos(theta),4) - 30*pow(cos(theta),2) + 3 );
}
std::complex<double> y4m1( double theta, double phi ){
	return - std::conj( y41(theta,phi) );
}
std::complex<double> y4m2( double theta, double phi ){
	return   std::conj( y42(theta,phi) );
}
std::complex<double> y4m3( double theta, double phi ){
	return - std::conj( y43(theta,phi) );
}
std::complex<double> y4m4( double theta, double phi ){
	return   std::conj( y44(theta,phi) );
}


std::complex<double> y55 ( double theta, double phi ){
	return -0.09375 * sqrt(77.0 /pi) * exp(std::complex<double>{0,5*phi}) * pow(sin(theta),5);
}
std::complex<double> y54 ( double theta, double phi ){
	return 0.1875   * sqrt(192.5/pi) * exp(std::complex<double>{0,4*phi}) * pow(sin(theta),4) *         cos(theta);
}
std::complex<double> y53 ( double theta, double phi ){
	return -0.03125 * sqrt(385.0/pi) * exp(std::complex<double>{0,3*phi}) * pow(sin(theta),3) * ( 9*pow(cos(theta),2) - 1 );
}
std::complex<double> y52 ( double theta, double phi ){
	return 0.125    * sqrt(577.5/pi) * exp(std::complex<double>{0,2*phi}) * pow(sin(theta),2) * ( 3*pow(cos(theta),3) - cos(theta) );
}
std::complex<double> y51 ( double theta, double phi ){
	return -0.0625  * sqrt(82.5 /pi) * exp(std::complex<double>{0,  phi}) *     sin(theta)    * ( 21*pow(cos(theta),4) - 14*pow(cos(theta),2) + 1 );
}
std::complex<double> y50 ( double theta, double phi ){
	return 0.0625   * sqrt(11.0 /pi)                                                          * ( 63*pow(cos(theta),5) - 70*pow(cos(theta),3) + 15*cos(theta) );
}
std::complex<double> y5m1( double theta, double phi ){
	return - std::conj( y51(theta,phi) );
}
std::complex<double> y5m2( double theta, double phi ){
	return   std::conj( y52(theta,phi) );
}
std::complex<double> y5m3( double theta, double phi ){
	return - std::conj( y53(theta,phi) );
}
std::complex<double> y5m4( double theta, double phi ){
	return   std::conj( y54(theta,phi) );
}
std::complex<double> y5m5( double theta, double phi ){
	return - std::conj( y55(theta,phi) );
}




//----- First theta-derivatives -----
std::complex<double> y00_dtheta ( double theta, double phi ){
	return 0;
}

std::complex<double> y11_dtheta ( double theta, double phi ){
	return -0.5*sqrt(1.5/pi) * exp(std::complex<double>{0,phi}) *  cos(theta);
}
std::complex<double> y10_dtheta ( double theta, double phi ){
	return 0.5 *sqrt(3.0/pi)                                    * -sin(theta);
}
std::complex<double> y1m1_dtheta( double theta, double phi ){
	return - std::conj( y11_dtheta(theta,phi) );
}

std::complex<double> y22_dtheta ( double theta, double phi ){
	return 0.25 * sqrt(7.5/pi) * exp(std::complex<double>{0,2*phi}) *      sin(2*theta);
}
std::complex<double> y21_dtheta ( double theta, double phi ){
	return -0.5 * sqrt(7.5/pi) * exp(std::complex<double>{0,  phi}) *      cos(2*theta);
}
std::complex<double> y20_dtheta ( double theta, double phi ){
	return 0.25 * sqrt(5.0/pi)                                      * -3.0*sin(2*theta);
}
std::complex<double> y2m1_dtheta( double theta, double phi ){
	return - std::conj( y21_dtheta(theta,phi) );
}
std::complex<double> y2m2_dtheta( double theta, double phi ){
	return   std::conj( y22_dtheta(theta,phi) );
}

std::complex<double> y33_dtheta ( double theta, double phi ){
	return -0.125 * sqrt(35.0/pi) * exp(std::complex<double>{0,3*phi})                     * ( -3*pow(cos(theta),3) + 3 *cos(theta) );
}
std::complex<double> y32_dtheta ( double theta, double phi ){
	return 0.25   * sqrt(52.5/pi) * exp(std::complex<double>{0,2*phi}) *     sin(theta)    * ( -3*pow(sin(theta),2) + 2             );
}
std::complex<double> y31_dtheta ( double theta, double phi ){
	return -0.125 * sqrt(21.0/pi) * exp(std::complex<double>{0,  phi})                     * ( 15*pow(cos(theta),3) - 11*cos(theta) );
}
std::complex<double> y30_dtheta ( double theta, double phi ){
	return 0.25*sqrt(7.0/pi) * ( 15*pow(sin(theta),3) - 12*sin(theta) );
}
std::complex<double> y3m1_dtheta( double theta, double phi ){
	return - std::conj( y31_dtheta(theta,phi) );
}
std::complex<double> y3m2_dtheta( double theta, double phi ){
	return   std::conj( y32_dtheta(theta,phi) );
}
std::complex<double> y3m3_dtheta( double theta, double phi ){
	return - std::conj( y33_dtheta(theta,phi) );
}

std::complex<double> y44_dtheta ( double theta, double phi ){
	return 0.1875 * sqrt(17.5/pi) * exp(std::complex<double>{0,4*phi}) * 4.0 *pow(sin(theta),3) *         cos(theta);
}
std::complex<double> y43_dtheta ( double theta, double phi ){
	return -0.375 * sqrt(35.0/pi) * exp(std::complex<double>{0,3*phi}) *      pow(sin(theta),2) * (-4*pow(sin(theta),2)+3);
}
std::complex<double> y42_dtheta ( double theta, double phi ){
	return 0.375  * sqrt(2.5 /pi) * exp(std::complex<double>{0,2*phi}) * 2.0 *    sin(theta*2 ) * (7 *pow(cos(theta),2)-4);
}
std::complex<double> y41_dtheta ( double theta, double phi ){
	return -0.375 * sqrt(5.0 /pi) * exp(std::complex<double>{0,  phi})                          * (28*pow(cos(theta),4)-27*pow(cos(theta),2)+3 );
}
std::complex<double> y40_dtheta ( double theta, double phi ){
	return 0.1875 / sqrt(     pi)                                      *  10.0*   sin(theta*2 ) * (-7*pow(cos(theta),2)+3);
}
std::complex<double> y4m1_dtheta( double theta, double phi ){
	return - std::conj( y41_dtheta(theta,phi) );
}
std::complex<double> y4m2_dtheta( double theta, double phi ){
	return   std::conj( y42_dtheta(theta,phi) );
}
std::complex<double> y4m3_dtheta( double theta, double phi ){
	return - std::conj( y43_dtheta(theta,phi) );
}
std::complex<double> y4m4_dtheta( double theta, double phi ){
	return   std::conj( y44_dtheta(theta,phi) );
}

std::complex<double> y55_dtheta ( double theta, double phi ){
	return -0.09375 * sqrt(77.0 /pi) * exp(std::complex<double>{0,5*phi}) * 5.0*pow(sin(theta),4) *          cos(theta);
}
std::complex<double> y54_dtheta ( double theta, double phi ){
	return 0.1875   * sqrt(192.5/pi) * exp(std::complex<double>{0,4*phi}) *     pow(sin(theta),3) * (-5 *pow(sin(theta),2)+4                       );
}
std::complex<double> y53_dtheta ( double theta, double phi ){
	return -0.03125 * sqrt(385.0/pi) * exp(std::complex<double>{0,3*phi}) *         cos(theta)    * (24 *pow(sin(theta),2)-45 *pow(sin(theta),4)   );
}
std::complex<double> y52_dtheta ( double theta, double phi ){
	return 0.125    * sqrt(577.5/pi) * exp(std::complex<double>{0,2*phi}) *         sin(theta)    * (15 *pow(cos(theta),4)-12 *pow(cos(theta),2)+1 );
}
std::complex<double> y51_dtheta ( double theta, double phi ){
	return -0.0625  * sqrt(82.5 /pi) * exp(std::complex<double>{0,  phi}) *         cos(theta)    * (105*pow(cos(theta),4)-126*pow(cos(theta),2)+29);
}
std::complex<double> y50_dtheta ( double theta, double phi ){
	return 0.0625*sqrt(11.0/pi) * 15*sin(theta)*(-21*pow(cos(theta),4)+14*pow(cos(theta),2)-1);
}
std::complex<double> y5m1_dtheta( double theta, double phi ){
	return - std::conj( y51_dtheta(theta,phi) );
}
std::complex<double> y5m2_dtheta( double theta, double phi ){
	return   std::conj( y52_dtheta(theta,phi) );
}
std::complex<double> y5m3_dtheta( double theta, double phi ){
	return - std::conj( y53_dtheta(theta,phi) );
}
std::complex<double> y5m4_dtheta( double theta, double phi ){
	return   std::conj( y54_dtheta(theta,phi) );
}
std::complex<double> y5m5_dtheta( double theta, double phi ){
	return - std::conj( y55_dtheta(theta,phi) );
}




//----- Second theta-derivatives -----
std::complex<double> y00_dtheta2( double theta, double phi ){
	return 0;
}

std::complex<double> y11_dtheta2 ( double theta, double phi ){
	return -0.5*sqrt(1.5/pi) * exp(std::complex<double>{0,phi}) * -sin(theta);
}
std::complex<double> y10_dtheta2 ( double theta, double phi ){
	return 0.5*sqrt(3.0/pi)                                     * -cos(theta);
}
std::complex<double> y1m1_dtheta2( double theta, double phi ){
	return - std::conj( y11_dtheta2(theta,phi) );
}

std::complex<double> y22_dtheta2 ( double theta, double phi ){
	return 0.25 * sqrt(7.5/pi) * exp(std::complex<double>{0,2*phi}) *  2.0*cos(2*theta);
}
std::complex<double> y21_dtheta2 ( double theta, double phi ){
	return -0.5 * sqrt(7.5/pi) * exp(std::complex<double>{0,  phi}) * -2.0*sin(2*theta);
}
std::complex<double> y20_dtheta2 ( double theta, double phi ){
	return 0.25 * sqrt(5.0/pi)                                      * -6.0*cos(2*theta);
}
std::complex<double> y2m1_dtheta2( double theta, double phi ){
	return - std::conj( y21_dtheta2(theta,phi) );
}
std::complex<double> y2m2_dtheta2( double theta, double phi ){
	return   std::conj( y22_dtheta2(theta,phi) );
}

std::complex<double> y33_dtheta2 ( double theta, double phi ){
	return -0.125 * sqrt(35.0/pi) * exp(std::complex<double>{0,3*phi}) *              ( -3 *pow(sin(theta),3) + 6*sin(theta)*pow(cos(theta),2) );
}
std::complex<double> y32_dtheta2 ( double theta, double phi ){
	return 0.25   * sqrt(52.5/pi) * exp(std::complex<double>{0,2*phi}) * cos(theta) * ( 9  *pow(cos(theta),2) - 7                              );
}
std::complex<double> y31_dtheta2 ( double theta, double phi ){
	return -0.125 * sqrt(21.0/pi) * exp(std::complex<double>{0,  phi}) *              ( 45 *pow(sin(theta),3) - 34*sin(theta)                  );
}
std::complex<double> y30_dtheta2 ( double theta, double phi ){
	return 0.25   * sqrt(7.0 /pi)                                                   * ( -45*pow(cos(theta),3) + 33*cos(theta)                  );
}
std::complex<double> y3m1_dtheta2( double theta, double phi ){
	return - std::conj( y31_dtheta2(theta,phi) );
}
std::complex<double> y3m2_dtheta2( double theta, double phi ){
	return   std::conj( y32_dtheta2(theta,phi) );
}
std::complex<double> y3m3_dtheta2( double theta, double phi ){
	return - std::conj( y33_dtheta2(theta,phi) );
}

std::complex<double> y44_dtheta2 ( double theta, double phi ){
	return 0.1875 * sqrt(17.5/pi) * exp(std::complex<double>{0,4*phi}) * 4.0 *pow(sin(theta),2) * (-4 *pow(sin(theta),2)+3                     );
}
std::complex<double> y43_dtheta2 ( double theta, double phi ){
	return -0.375 * sqrt(35.0/pi) * exp(std::complex<double>{0,3*phi}) *          sin(theta*2 ) * (-8 *pow(sin(theta),2)+3                     );
}
std::complex<double> y42_dtheta2 ( double theta, double phi ){
	return 0.375  * sqrt(2.5 /pi) * exp(std::complex<double>{0,2*phi}) * 4.0                    * (28 *pow(cos(theta),4)-29*pow(cos(theta),2)+4);
}
std::complex<double> y41_dtheta2 ( double theta, double phi ){
	return -0.375 * sqrt(5.0 /pi) * exp(std::complex<double>{0,  phi}) *          sin(theta*2 ) * (-56*pow(cos(theta),2)+27                    );
}
std::complex<double> y40_dtheta2 ( double theta, double phi ){
	return 0.1875 / sqrt(     pi)                                      * 20.0                   * (-28*pow(cos(theta),4)+27*pow(cos(theta),2)-3);
}
std::complex<double> y4m1_dtheta2( double theta, double phi ){
	return - std::conj( y41_dtheta2(theta,phi) );
}
std::complex<double> y4m2_dtheta2( double theta, double phi ){
	return   std::conj( y42_dtheta2(theta,phi) );
}
std::complex<double> y4m3_dtheta2( double theta, double phi ){
	return - std::conj( y43_dtheta2(theta,phi) );
}
std::complex<double> y4m4_dtheta2( double theta, double phi ){
	return   std::conj( y44_dtheta2(theta,phi) );
}

std::complex<double> y55_dtheta2 ( double theta, double phi ){
	return -0.09375 * sqrt(77.0 /pi) * exp(std::complex<double>{0,5*phi}) * 5.0 *pow(sin(theta),3)            * (-5  *pow(sin(theta),2)+4                       );
}
std::complex<double> y54_dtheta2 ( double theta, double phi ){
	return 0.1875   * sqrt(192.5/pi) * exp(std::complex<double>{0,4*phi}) *      pow(sin(theta),2)*cos(theta) * (-25 *pow(sin(theta),2)+12                      );
}
std::complex<double> y53_dtheta2 ( double theta, double phi ){
	return -0.03125 * sqrt(385.0/pi) * exp(std::complex<double>{0,3*phi}) * 3.0 *    sin(theta)               * (75  *pow(sin(theta),4)-84 *pow(sin(theta),2)+16);
}
std::complex<double> y52_dtheta2 ( double theta, double phi ){
	return 0.125    * sqrt(577.5/pi) * exp(std::complex<double>{0,2*phi}) *          cos(theta)               * (75  *pow(cos(theta),4)-96 *pow(cos(theta),2)+25);
}
std::complex<double> y51_dtheta2 ( double theta, double phi ){
	return -0.0625  * sqrt(82.5 /pi) * exp(std::complex<double>{0,  phi}) *          sin(theta)               * (-525*pow(cos(theta),4)+378*pow(cos(theta),2)-29);
}
std::complex<double> y50_dtheta2 ( double theta, double phi ){
	return 0.0625   * sqrt(11.0 /pi)                                      * 15.0*    cos(theta)               * (-105*pow(cos(theta),4)+126*pow(cos(theta),2)-29);
}
std::complex<double> y5m1_dtheta2( double theta, double phi ){
	return - std::conj( y51_dtheta2(theta,phi) );
}
std::complex<double> y5m2_dtheta2( double theta, double phi ){
	return   std::conj( y52_dtheta2(theta,phi) );
}
std::complex<double> y5m3_dtheta2( double theta, double phi ){
	return - std::conj( y53_dtheta2(theta,phi) );
}
std::complex<double> y5m4_dtheta2( double theta, double phi ){
	return   std::conj( y54_dtheta2(theta,phi) );
}
std::complex<double> y5m5_dtheta2( double theta, double phi ){
	return - std::conj( y55_dtheta2(theta,phi) );
}




//----- First phi-derivatives -----
double y00_dphi               ( double theta, double phi ){
	return 0;
}


std::complex<double> y11_dphi ( double theta, double phi ){
	return -0.5*sqrt(1.5/pi) * exp(std::complex<double>{0,phi+pi/2.0}) * sin(theta);
}
double               y10_dphi ( double theta, double phi ){
	return 0;
}
std::complex<double> y1m1_dphi( double theta, double phi ){
	return - std::conj( y11_dphi(theta,phi) );
}


std::complex<double> y22_dphi ( double theta, double phi ){
	return 0.25 * sqrt(7.5/pi) * 2 * exp(std::complex<double>{0,2*phi+pi/2.0}) * pow(sin(theta),2);
}
std::complex<double> y21_dphi ( double theta, double phi ){
	return -0.5 * sqrt(7.5/pi)     * exp(std::complex<double>{0,  phi+pi/2.0}) *     sin(theta)    *         cos(theta);
}
double               y20_dphi ( double theta, double phi ){
	return 0;
}
std::complex<double> y2m1_dphi( double theta, double phi ){
	return - std::conj( y21_dphi(theta,phi) );
}
std::complex<double> y2m2_dphi( double theta, double phi ){
	return   std::conj( y22_dphi(theta,phi) );
}


std::complex<double> y33_dphi ( double theta, double phi ){
	return -0.125 * sqrt(35.0/pi) * 3 * exp(std::complex<double>{0,3*phi+pi/2.0}) * pow(sin(theta),3);
}
std::complex<double> y32_dphi ( double theta, double phi ){
	return 0.25   * sqrt(52.5/pi) * 2 * exp(std::complex<double>{0,2*phi+pi/2.0}) * pow(sin(theta),2) *         cos(theta);
}
std::complex<double> y31_dphi ( double theta, double phi ){
	return -0.125 * sqrt(21.0/pi)     * exp(std::complex<double>{0,  phi+pi/2.0}) *     sin(theta)    * ( 5*pow(cos(theta),2) - 1 );
}
double               y30_dphi ( double theta, double phi ){
	return 0;
}
std::complex<double> y3m1_dphi( double theta, double phi ){
	return - std::conj( y31_dphi(theta,phi) );
}
std::complex<double> y3m2_dphi( double theta, double phi ){
	return   std::conj( y32_dphi(theta,phi) );
}
std::complex<double> y3m3_dphi( double theta, double phi ){
	return - std::conj( y33_dphi(theta,phi) );
}


std::complex<double> y44_dphi ( double theta, double phi ){
	return 0.1875 * sqrt(17.5/pi) * 4 * exp(std::complex<double>{0,4*phi+pi/2.0}) * pow(sin(theta),4);
}
std::complex<double> y43_dphi ( double theta, double phi ){
	return -0.375 * sqrt(35.0/pi) * 3 * exp(std::complex<double>{0,3*phi+pi/2.0}) * pow(sin(theta),3) *          cos(theta);
}
std::complex<double> y42_dphi ( double theta, double phi ){
	return 0.375  * sqrt(2.5 /pi) * 2 * exp(std::complex<double>{0,2*phi+pi/2.0}) * pow(sin(theta),2) * ( 7 *pow(cos(theta),2) - 1 );
}
std::complex<double> y41_dphi ( double theta, double phi ){
	return -0.375 * sqrt(5.0 /pi)     * exp(std::complex<double>{0,  phi+pi/2.0}) *     sin(theta)    * ( 7 *pow(cos(theta),3) - 3*cos(theta) );
}
double               y40_dphi ( double theta, double phi ){
	return 0;
}
std::complex<double> y4m1_dphi( double theta, double phi ){
	return - std::conj( y41_dphi(theta,phi) );
}
std::complex<double> y4m2_dphi( double theta, double phi ){
	return   std::conj( y42_dphi(theta,phi) );
}
std::complex<double> y4m3_dphi( double theta, double phi ){
	return - std::conj( y43_dphi(theta,phi) );
}
std::complex<double> y4m4_dphi( double theta, double phi ){
	return   std::conj( y44_dphi(theta,phi) );
}


std::complex<double> y55_dphi ( double theta, double phi ){
	return -0.09375 * sqrt(77.0 /pi) * 5 * exp(std::complex<double>{0,5*phi+pi/2.0}) * pow(sin(theta),5);
}
std::complex<double> y54_dphi ( double theta, double phi ){
	return 0.1875   * sqrt(192.5/pi) * 4 * exp(std::complex<double>{0,4*phi+pi/2.0}) * pow(sin(theta),4) *         cos(theta);
}
std::complex<double> y53_dphi ( double theta, double phi ){
	return -0.03125 * sqrt(385.0/pi) * 3 * exp(std::complex<double>{0,3*phi+pi/2.0}) * pow(sin(theta),3) * ( 9*pow(cos(theta),2) - 1 );
}
std::complex<double> y52_dphi ( double theta, double phi ){
	return 0.125    * sqrt(577.5/pi) * 2 * exp(std::complex<double>{0,2*phi+pi/2.0}) * pow(sin(theta),2) * ( 3*pow(cos(theta),3) - cos(theta) );
}
std::complex<double> y51_dphi ( double theta, double phi ){
	return -0.0625  * sqrt(82.5 /pi) *     exp(std::complex<double>{0,  phi+pi/2.0}) *     sin(theta)    * ( 21*pow(cos(theta),4) - 14*pow(cos(theta),2) + 1 );
}
double               y50_dphi ( double theta, double phi ){
	return 0;
}
std::complex<double> y5m1_dphi( double theta, double phi ){
	return - std::conj( y51_dphi(theta,phi) );
}
std::complex<double> y5m2_dphi( double theta, double phi ){
	return   std::conj( y52_dphi(theta,phi) );
}
std::complex<double> y5m3_dphi( double theta, double phi ){
	return - std::conj( y53_dphi(theta,phi) );
}
std::complex<double> y5m4_dphi( double theta, double phi ){
	return   std::conj( y54_dphi(theta,phi) );
}
std::complex<double> y5m5_dphi( double theta, double phi ){
	return - std::conj( y55_dphi(theta,phi) );
}




//----- Second phi-derivatives -----
double y00_dphi2               ( double theta, double phi ){
	return 0;
}


std::complex<double> y11_dphi2 ( double theta, double phi ){
	return -0.5*sqrt(1.5/pi) * -1 * exp(std::complex<double>{0,phi}) * sin(theta);
}
double               y10_dphi2 ( double theta, double phi ){
	return 0;
}
std::complex<double> y1m1_dphi2( double theta, double phi ){
	return - std::conj( y11_dphi2(theta,phi) );
}


std::complex<double> y22_dphi2 ( double theta, double phi ){
	return 0.25 * sqrt(7.5/pi) * -4 * exp(std::complex<double>{0,2*phi}) * pow(sin(theta),2);
}
std::complex<double> y21_dphi2 ( double theta, double phi ){
	return -0.5 * sqrt(7.5/pi) * -1 * exp(std::complex<double>{0,  phi}) *     sin(theta)    *         cos(theta);
}
double               y20_dphi2 ( double theta, double phi ){
	return 0;
}
std::complex<double> y2m1_dphi2( double theta, double phi ){
	return - std::conj( y21_dphi2(theta,phi) );
}
std::complex<double> y2m2_dphi2( double theta, double phi ){
	return   std::conj( y22_dphi2(theta,phi) );
}


std::complex<double> y33_dphi2 ( double theta, double phi ){
	return -0.125 * sqrt(35.0/pi) * -9 * exp(std::complex<double>{0,3*phi}) * pow(sin(theta),3);
}
std::complex<double> y32_dphi2 ( double theta, double phi ){
	return 0.25   * sqrt(52.5/pi) * -4 * exp(std::complex<double>{0,2*phi}) * pow(sin(theta),2) *         cos(theta);
}
std::complex<double> y31_dphi2 ( double theta, double phi ){
	return -0.125 * sqrt(21.0/pi) * -1 * exp(std::complex<double>{0,  phi}) *     sin(theta)    * ( 5*pow(cos(theta),2) - 1 );
}
double               y30_dphi2 ( double theta, double phi ){
	return 0;
}
std::complex<double> y3m1_dphi2( double theta, double phi ){
	return - std::conj( y31_dphi2(theta,phi) );
}
std::complex<double> y3m2_dphi2( double theta, double phi ){
	return   std::conj( y32_dphi2(theta,phi) );
}
std::complex<double> y3m3_dphi2( double theta, double phi ){
	return - std::conj( y33_dphi2(theta,phi) );
}


std::complex<double> y44_dphi2 ( double theta, double phi ){
	return 0.1875 * sqrt(17.5/pi) * -16 * exp(std::complex<double>{0,4*phi}) * pow(sin(theta),4);
}
std::complex<double> y43_dphi2 ( double theta, double phi ){
	return -0.375 * sqrt(35.0/pi) * -9  * exp(std::complex<double>{0,3*phi}) * pow(sin(theta),3) *          cos(theta);
}
std::complex<double> y42_dphi2 ( double theta, double phi ){
	return 0.375  * sqrt(2.5 /pi) * -4  * exp(std::complex<double>{0,2*phi}) * pow(sin(theta),2) * ( 7 *pow(cos(theta),2) - 1 );
}
std::complex<double> y41_dphi2 ( double theta, double phi ){
	return -0.375 * sqrt(5.0 /pi) * -1  * exp(std::complex<double>{0,  phi}) *     sin(theta)    * ( 7 *pow(cos(theta),3) - 3*cos(theta) );
}
double               y40_dphi2 ( double theta, double phi ){
	return 0;
}
std::complex<double> y4m1_dphi2( double theta, double phi ){
	return - std::conj( y41_dphi2(theta,phi) );
}
std::complex<double> y4m2_dphi2( double theta, double phi ){
	return   std::conj( y42_dphi2(theta,phi) );
}
std::complex<double> y4m3_dphi2( double theta, double phi ){
	return - std::conj( y43_dphi2(theta,phi) );
}
std::complex<double> y4m4_dphi2( double theta, double phi ){
	return   std::conj( y44_dphi2(theta,phi) );
}


std::complex<double> y55_dphi2 ( double theta, double phi ){
	return -0.09375 * sqrt(77.0 /pi) * -25 * exp(std::complex<double>{0,5*phi}) * pow(sin(theta),5);
}
std::complex<double> y54_dphi2 ( double theta, double phi ){
	return 0.1875   * sqrt(192.5/pi) * -16 * exp(std::complex<double>{0,4*phi}) * pow(sin(theta),4) *         cos(theta);
}
std::complex<double> y53_dphi2 ( double theta, double phi ){
	return -0.03125 * sqrt(385.0/pi) * -9  * exp(std::complex<double>{0,3*phi}) * pow(sin(theta),3) * ( 9*pow(cos(theta),2) - 1 );
}
std::complex<double> y52_dphi2 ( double theta, double phi ){
	return 0.125    * sqrt(577.5/pi) * -4  * exp(std::complex<double>{0,2*phi}) * pow(sin(theta),2) * ( 3*pow(cos(theta),3) - cos(theta) );
}
std::complex<double> y51_dphi2 ( double theta, double phi ){
	return -0.0625  * sqrt(82.5 /pi) * -1  * exp(std::complex<double>{0,  phi}) *     sin(theta)    * ( 21*pow(cos(theta),4) - 14*pow(cos(theta),2) + 1 );
}
double               y50_dphi2 ( double theta, double phi ){
	return 0;
}
std::complex<double> y5m1_dphi2( double theta, double phi ){
	return - std::conj( y51_dphi2(theta,phi) );
}
std::complex<double> y5m2_dphi2( double theta, double phi ){
	return   std::conj( y52_dphi2(theta,phi) );
}
std::complex<double> y5m3_dphi2( double theta, double phi ){
	return - std::conj( y53_dphi2(theta,phi) );
}
std::complex<double> y5m4_dphi2( double theta, double phi ){
	return   std::conj( y54_dphi2(theta,phi) );
}
std::complex<double> y5m5_dphi2( double theta, double phi ){
	return - std::conj( y55_dphi2(theta,phi) );
}




//----- Spherical harmonics divided by sin(theta) -----
std::complex<double> y00_over_sintheta ( double theta, double phi ){
	return 0.5/sqrt(pi) / sin(theta);
}


std::complex<double> y11_over_sintheta ( double theta, double phi ){
	return -0.5*sqrt(1.5/pi) * exp(std::complex<double>{0,phi});
}
std::complex<double> y10_over_sintheta ( double theta, double phi ){
	return 0.5 *sqrt(3.0/pi)                                    * tan( 0.5*pi - theta );	// cos(x)/sin(x) = cot(x) = tan(pi/2-x) and C++ does not feature a cotangent function
}
std::complex<double> y1m1_over_sintheta( double theta, double phi ){
	return - std::conj( y11_over_sintheta(theta,phi) );
}


std::complex<double> y22_over_sintheta ( double theta, double phi ){
	return 0.25 * sqrt(7.5/pi) * exp(std::complex<double>{0,2*phi}) *     sin(theta)   ;
}
std::complex<double> y21_over_sintheta ( double theta, double phi ){
	return -0.5 * sqrt(7.5/pi) * exp(std::complex<double>{0,  phi})                     *         cos(theta);
}
std::complex<double> y20_over_sintheta ( double theta, double phi ){
	return 0.25 * sqrt(5.0/pi)                                                          * ( 3*pow(cos(theta),2) - 1 ) / sin(theta);
}
std::complex<double> y2m1_over_sintheta( double theta, double phi ){
	return - std::conj( y21_over_sintheta(theta,phi) );
}
std::complex<double> y2m2_over_sintheta( double theta, double phi ){
	return   std::conj( y22_over_sintheta(theta,phi) );
}


std::complex<double> y33_over_sintheta ( double theta, double phi ){
	return -0.125 * sqrt(35.0/pi) * exp(std::complex<double>{0,3*phi}) * pow(sin(theta),2);
}
std::complex<double> y32_over_sintheta ( double theta, double phi ){
	return 0.25   * sqrt(52.5/pi) * exp(std::complex<double>{0,2*phi}) *     sin(theta)    *         cos(theta);
}
std::complex<double> y31_over_sintheta ( double theta, double phi ){
	return -0.125 * sqrt(21.0/pi) * exp(std::complex<double>{0,  phi})                     * ( 5*pow(cos(theta),2) - 1 );
}
std::complex<double> y30_over_sintheta ( double theta, double phi ){
	return 0.25   * sqrt(7.0 /pi)                                                          * ( 5*pow(cos(theta),3) - 3*cos(theta) ) / sin(theta);
}
std::complex<double> y3m1_over_sintheta( double theta, double phi ){
	return - std::conj( y31_over_sintheta(theta,phi) );
}
std::complex<double> y3m2_over_sintheta( double theta, double phi ){
	return   std::conj( y32_over_sintheta(theta,phi) );
}
std::complex<double> y3m3_over_sintheta( double theta, double phi ){
	return - std::conj( y33_over_sintheta(theta,phi) );
}


std::complex<double> y44_over_sintheta ( double theta, double phi ){
	return 0.1875 * sqrt(17.5/pi) * exp(std::complex<double>{0,4*phi}) * pow(sin(theta),3);
}
std::complex<double> y43_over_sintheta ( double theta, double phi ){
	return -0.375 * sqrt(35.0/pi) * exp(std::complex<double>{0,3*phi}) * pow(sin(theta),2) *          cos(theta);
}
std::complex<double> y42_over_sintheta ( double theta, double phi ){
	return 0.375  * sqrt(2.5 /pi) * exp(std::complex<double>{0,2*phi}) *     sin(theta)    * ( 7 *pow(cos(theta),2) - 1 );
}
std::complex<double> y41_over_sintheta ( double theta, double phi ){
	return -0.375 * sqrt(5.0 /pi) * exp(std::complex<double>{0,  phi})                     * ( 7 *pow(cos(theta),3) - 3*cos(theta) );
}
std::complex<double> y40_over_sintheta ( double theta, double phi ){
	return 0.1875 / sqrt(     pi)                                                          * ( 35*pow(cos(theta),4) - 30*pow(cos(theta),2) + 3 ) / sin(theta);
}
std::complex<double> y4m1_over_sintheta( double theta, double phi ){
	return - std::conj( y41_over_sintheta(theta,phi) );
}
std::complex<double> y4m2_over_sintheta( double theta, double phi ){
	return   std::conj( y42_over_sintheta(theta,phi) );
}
std::complex<double> y4m3_over_sintheta( double theta, double phi ){
	return - std::conj( y43_over_sintheta(theta,phi) );
}
std::complex<double> y4m4_over_sintheta( double theta, double phi ){
	return   std::conj( y44_over_sintheta(theta,phi) );
}


std::complex<double> y55_over_sintheta ( double theta, double phi ){
	return -0.09375 * sqrt(77.0 /pi) * exp(std::complex<double>{0,5*phi}) * pow(sin(theta),4);
}
std::complex<double> y54_over_sintheta ( double theta, double phi ){
	return 0.1875   * sqrt(192.5/pi) * exp(std::complex<double>{0,4*phi}) * pow(sin(theta),3) *         cos(theta);
}
std::complex<double> y53_over_sintheta ( double theta, double phi ){
	return -0.03125 * sqrt(385.0/pi) * exp(std::complex<double>{0,3*phi}) * pow(sin(theta),2) * ( 9*pow(cos(theta),2) - 1 );
}
std::complex<double> y52_over_sintheta ( double theta, double phi ){
	return 0.125    * sqrt(577.5/pi) * exp(std::complex<double>{0,2*phi}) *     sin(theta)    * ( 3*pow(cos(theta),3) - cos(theta) );
}
std::complex<double> y51_over_sintheta ( double theta, double phi ){
	return -0.0625  * sqrt(82.5 /pi) * exp(std::complex<double>{0,  phi})                     * ( 21*pow(cos(theta),4) - 14*pow(cos(theta),2) + 1 );
}
std::complex<double> y50_over_sintheta ( double theta, double phi ){
	return 0.0625   * sqrt(11.0 /pi)                                                          * ( 63*pow(cos(theta),5) - 70*pow(cos(theta),3) + 15*cos(theta) ) / sin(theta);
}
std::complex<double> y5m1_over_sintheta( double theta, double phi ){
	return - std::conj( y51_over_sintheta(theta,phi) );
}
std::complex<double> y5m2_over_sintheta( double theta, double phi ){
	return   std::conj( y52_over_sintheta(theta,phi) );
}
std::complex<double> y5m3_over_sintheta( double theta, double phi ){
	return - std::conj( y53_over_sintheta(theta,phi) );
}
std::complex<double> y5m4_over_sintheta( double theta, double phi ){
	return   std::conj( y54_over_sintheta(theta,phi) );
}
std::complex<double> y5m5_over_sintheta( double theta, double phi ){
	return - std::conj( y55_over_sintheta(theta,phi) );
}




//----- Real spherical harmonics -----
double y00_real( double theta, double phi ){ return y00(theta,phi).real(); }
double y10_real( double theta, double phi ){ return y10(theta,phi).real(); }
double y20_real( double theta, double phi ){ return y20(theta,phi).real(); }
double y30_real( double theta, double phi ){ return y30(theta,phi).real(); }
double y40_real( double theta, double phi ){ return y40(theta,phi).real(); }
double y50_real( double theta, double phi ){ return y50(theta,phi).real(); }

double y11_real ( double theta, double phi ){ return (                             sqrt(0.5) * ( y1m1(theta,phi) - y11(theta,phi) ) ).real(); }
double y1m1_real( double theta, double phi ){ return ( std::complex<double>{0,1} * sqrt(0.5) * ( y1m1(theta,phi) + y11(theta,phi) ) ).real(); }

double y22_real ( double theta, double phi ){ return (                             sqrt(0.5) * ( y2m2(theta,phi) + y22(theta,phi) ) ).real(); }
double y21_real ( double theta, double phi ){ return (                             sqrt(0.5) * ( y2m1(theta,phi) - y21(theta,phi) ) ).real(); }
double y2m1_real( double theta, double phi ){ return ( std::complex<double>{0,1} * sqrt(0.5) * ( y2m1(theta,phi) + y21(theta,phi) ) ).real(); }
double y2m2_real( double theta, double phi ){ return ( std::complex<double>{0,1} * sqrt(0.5) * ( y2m2(theta,phi) - y22(theta,phi) ) ).real(); }

double y33_real ( double theta, double phi ){ return (                             sqrt(0.5) * ( y3m3(theta,phi) - y33(theta,phi) ) ).real(); }
double y32_real ( double theta, double phi ){ return (                             sqrt(0.5) * ( y3m2(theta,phi) + y32(theta,phi) ) ).real(); }
double y31_real ( double theta, double phi ){ return (                             sqrt(0.5) * ( y3m1(theta,phi) - y31(theta,phi) ) ).real(); }
double y3m1_real( double theta, double phi ){ return ( std::complex<double>{0,1} * sqrt(0.5) * ( y3m1(theta,phi) + y31(theta,phi) ) ).real(); }
double y3m2_real( double theta, double phi ){ return ( std::complex<double>{0,1} * sqrt(0.5) * ( y3m2(theta,phi) - y32(theta,phi) ) ).real(); }
double y3m3_real( double theta, double phi ){ return ( std::complex<double>{0,1} * sqrt(0.5) * ( y3m3(theta,phi) + y33(theta,phi) ) ).real(); }

double y44_real ( double theta, double phi ){ return (                             sqrt(0.5) * ( y4m4(theta,phi) + y44(theta,phi) ) ).real(); }
double y43_real ( double theta, double phi ){ return (                             sqrt(0.5) * ( y4m3(theta,phi) - y43(theta,phi) ) ).real(); }
double y42_real ( double theta, double phi ){ return (                             sqrt(0.5) * ( y4m2(theta,phi) + y42(theta,phi) ) ).real(); }
double y41_real ( double theta, double phi ){ return (                             sqrt(0.5) * ( y4m1(theta,phi) - y41(theta,phi) ) ).real(); }
double y4m1_real( double theta, double phi ){ return ( std::complex<double>{0,1} * sqrt(0.5) * ( y4m1(theta,phi) + y41(theta,phi) ) ).real(); }
double y4m2_real( double theta, double phi ){ return ( std::complex<double>{0,1} * sqrt(0.5) * ( y4m2(theta,phi) - y42(theta,phi) ) ).real(); }
double y4m3_real( double theta, double phi ){ return ( std::complex<double>{0,1} * sqrt(0.5) * ( y4m3(theta,phi) + y43(theta,phi) ) ).real(); }
double y4m4_real( double theta, double phi ){ return ( std::complex<double>{0,1} * sqrt(0.5) * ( y4m4(theta,phi) - y44(theta,phi) ) ).real(); }