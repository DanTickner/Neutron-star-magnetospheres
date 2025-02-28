// cd OneDrive\PhD\Codes\Test mathematical expressions
// g++ 20220908_atan2_arguments_multiplied_by_constant.cpp -o 20220908_atan2_arguments_multiplied_by_constant

// atan2( ay, ax ) = atan2(y,x) if a>0, atan2(y,x)-pi if a<0 and ygeq0, atan2(y,x)+pi if a<0 and y<0, undefined if a=0

#include <iostream>
#include <math.h>

int main(){
	
	//----- Define variables -----
	const double pi = 3.14159265358979323846;
	double x_pos =  0.5;
	double x_nil =  0  ;
	double x_neg = -0.3;
	double y_pos =  0.5;
	double y_nil =  0  ;
	double y_neg = -0.3;
	
	double a =  2.34923;
	double b = -2.34923;
	
	std::cout.precision( 3 );
	
	
	//----- atan2( y, x ) -----
	std::cout <<         atan2(  y_pos,  x_pos ) <<"\t"<< atan2(  y_pos,  x_nil ) <<"\t"<< atan2(  y_pos,  x_neg ) <<"\t"<< 
	                     atan2(  y_nil,  x_pos ) <<"\t"<< atan2(  y_nil,  x_nil ) <<"\t"<< atan2(  y_nil,  x_neg ) <<"\t"<<
				         atan2(  y_neg,  x_pos ) <<"\t"<< atan2(  y_neg,  x_nil ) <<"\t"<< atan2(  y_neg,  x_neg ) <<       std::endl;
	
	//----- atan2( y, ax ) -----
	std::cout << "\n" << atan2( y_pos, a*x_pos ) <<"\t"<< atan2( y_pos, a*x_nil ) <<"\t"<< atan2( y_pos, a*x_neg ) <<"\t"<< 
	                     atan2( y_nil, a*x_pos ) <<"\t"<< atan2( y_nil, a*x_nil ) <<"\t"<< atan2( y_nil, a*x_neg ) <<"\t"<<
				         atan2( y_neg, a*x_pos ) <<"\t"<< atan2( y_neg, a*x_nil ) <<"\t"<< atan2( y_neg, a*x_neg ) <<       std::endl;
						 
	/*std::cout <<         atan2(  y_pos,  x_pos ) - pi <<"\t"<< atan2(  y_pos,  x_nil ) - pi <<"\t"<< atan2(  y_pos,  x_neg ) - pi <<"\t"<< 
	                     atan2(  y_nil,  x_pos ) - pi <<"\t"<< atan2(  y_nil,  x_nil ) - pi <<"\t"<< atan2(  y_nil,  x_neg ) - pi <<"\t"<<
				         atan2(  y_neg,  x_pos ) + pi <<"\t"<< atan2(  y_neg,  x_nil ) + pi <<"\t"<< atan2(  y_neg,  x_neg ) + pi <<       std::endl;*/
						 
	
	//----- atan2( by, bx ) -----
	std::cout << "\n" << atan2( b*y_pos, b*x_pos ) <<"\t"<< atan2( b*y_pos, b*x_nil ) <<"\t"<< atan2( b*y_pos, b*x_neg ) <<"\t"<< 
	                     atan2( b*y_nil, b*x_pos ) <<"\t"<< atan2( b*y_nil, b*x_nil ) <<"\t"<< atan2( b*y_nil, b*x_neg ) <<"\t"<<
				         atan2( b*y_neg, b*x_pos ) <<"\t"<< atan2( b*y_neg, b*x_nil ) <<"\t"<< atan2( b*y_neg, b*x_neg ) <<       std::endl;
						 
	std::cout         << atan2( y_pos, x_pos ) - pi <<"\t"<< atan2( y_pos, x_nil ) - pi <<"\t"<< atan2( y_pos, x_neg ) - pi <<"\t"<< 
	                     atan2( y_nil, x_pos ) - pi <<"\t"<< atan2( y_nil, x_nil ) - pi <<"\t"<< atan2( y_nil, x_neg ) - pi <<"\t"<<
				         atan2( y_neg, x_pos ) + pi <<"\t"<< atan2( y_neg, x_nil ) + pi <<"\t"<< atan2( y_neg, x_neg ) + pi <<       std::endl;
	
	return 0;
}