// cd OneDrive\PhD\Codes\Test mathematical expressions
// g++ 20220908_atan2_negative_arguments.cpp -o 20220908_atan2_negative_arguments

/*
atan2( -y, -x ) = atan2( y, x ) - pi if y geq 0 else atan( y, x ) + pi
atan2(  y, -x ) = ?
atan2( -y,  x ) = ??
*/

#include <iostream>
#include <math.h>

//----- Function definitions -----
const double pi = 3.14159265358979323846;

double atan2_neg_neg( double y, double x );
double atan2_pos_neg( double y, double x );
double atan2_neg_pos( double y, double x );

int main(){
	
	//----- Define variables -----
	double x_pos =  0.5;
	double x_nil =  0  ;
	double x_neg = -0.3;
	double y_pos =  0.5;
	double y_nil =  0  ;
	double y_neg = -0.3;
	
	std::cout.precision( 3 );
	
	
	//----- atan2( y, x ) -----
	std::cout <<         atan2        (  y_pos,  x_pos )      <<"\t"<< atan2        (  y_nil,  x_pos )      <<"\t"<< atan2        (  y_neg,  x_pos )      <<"\t"<< 
	                     atan2        (  y_pos,  x_nil )      <<"\t"<< atan2        (  y_nil,  x_nil )      <<"\t"<< atan2        (  y_neg,  x_nil )      <<"\t"<<
				         atan2        (  y_pos,  x_neg )      <<"\t"<< atan2        (  y_nil,  x_neg )      <<"\t"<< atan2        (  y_neg,  x_neg )      <<       std::endl;
	
	//----- atan2( -y, -x ) -----
	std::cout << "\n" << atan2        ( -y_pos, -x_pos )      <<"\t"<< atan2        ( -y_nil, -x_pos )      <<"\t"<< atan2        ( -y_neg, -x_pos )      <<"\t"<< 
	                     atan2        ( -y_pos, -x_nil )      <<"\t"<< atan2        ( -y_nil, -x_nil )      <<"\t"<< atan2        ( -y_neg, -x_nil )      <<"\t"<<
				         atan2        ( -y_pos, -x_neg )      <<"\t"<< atan2        ( -y_nil, -x_neg )      <<"\t"<< atan2        ( -y_neg, -x_neg )      <<       std::endl;
						 
	std::cout <<         atan2        (  y_pos,  x_pos ) - pi <<"\t"<< atan2        (  y_nil,  x_pos ) - pi <<"\t"<< atan2        (  y_neg,  x_pos ) + pi <<"\t"<< 
	                     atan2        (  y_pos,  x_nil ) - pi <<"\t"<< atan2        (  y_nil,  x_nil ) - pi <<"\t"<< atan2        (  y_neg,  x_nil ) + pi <<"\t"<<
				         atan2        (  y_pos,  x_neg ) - pi <<"\t"<< atan2        (  y_nil,  x_neg ) - pi <<"\t"<< atan2        (  y_neg,  x_neg ) + pi <<       std::endl;
	
	std::cout <<         atan2_neg_neg(  y_pos,  x_pos )      <<"\t"<< atan2_neg_neg(  y_nil,  x_pos )      <<"\t"<< atan2_neg_neg(  y_neg,  x_pos )      <<"\t"<< 
						 atan2_neg_neg(  y_pos,  x_nil )      <<"\t"<< atan2_neg_neg(  y_nil,  x_nil )      <<"\t"<< atan2_neg_neg(  y_neg,  x_nil )      <<"\t"<< 
						 atan2_neg_neg(  y_pos,  x_neg )      <<"\t"<< atan2_neg_neg(  y_nil,  x_neg )      <<"\t"<< atan2_neg_neg(  y_neg,  x_neg )      <<       std::endl;
	
	/*//----- atan2( y, -x ) -----
	std::cout << "\n" << atan2        (  y_pos, -x_pos )      <<"\t"<< atan2        (  y_nil, -x_pos )      <<"\t"<< atan2        (  y_neg, -x_pos )      <<"\t"<< 
	                     atan2        (  y_pos, -x_nil )      <<"\t"<< atan2        (  y_nil, -x_nil )      <<"\t"<< atan2        (  y_neg, -x_nil )      <<"\t"<<
				         atan2        (  y_pos, -x_neg )      <<"\t"<< atan2        (  y_nil, -x_neg )      <<"\t"<< atan2        (  y_neg, -x_neg )      <<       std::endl;
	
	std::cout <<         -atan2       (  y_pos,  x_pos ) + pi <<"\t"<< -atan2       (  y_nil,  x_pos ) + pi <<"\t"<< -atan2       (  y_neg,  x_pos ) - pi <<"\t"<< 
	                     -atan2       (  y_pos,  x_nil ) + pi <<"\t"<< -atan2       (  y_nil,  x_nil ) + pi <<"\t"<< -atan2       (  y_neg,  x_nil ) - pi <<"\t"<<
				         -atan2       (  y_pos,  x_neg ) + pi <<"\t"<< -atan2       (  y_nil,  x_neg ) + pi <<"\t"<< -atan2       (  y_neg,  x_neg ) - pi <<       std::endl;
	
	std::cout <<         atan2_pos_neg(  y_pos,  x_pos )      <<"\t"<< atan2_pos_neg(  y_nil,  x_pos )      <<"\t"<< atan2_pos_neg(  y_neg,  x_pos )      <<"\t"<< 
						 atan2_pos_neg(  y_pos,  x_nil )      <<"\t"<< atan2_pos_neg(  y_nil,  x_nil )      <<"\t"<< atan2_pos_neg(  y_neg,  x_nil )      <<"\t"<< 
						 atan2_pos_neg(  y_pos,  x_neg )      <<"\t"<< atan2_pos_neg(  y_nil,  x_neg )      <<"\t"<< atan2_pos_neg(  y_neg,  x_neg )      <<       std::endl;*/
	
	/*//----- atan2( -y, x ) -----
	std::cout << "\n" << atan2        (  -y_pos, x_pos )      <<"\t"<< atan2        ( -y_nil, x_pos )      <<"\t"<< atan2        ( -y_neg, x_pos )      <<"\t"<< 
	                     atan2        (  -y_pos, x_nil )      <<"\t"<< atan2        ( -y_nil, x_nil )      <<"\t"<< atan2        ( -y_neg, x_nil )      <<"\t"<<
				         atan2        (  -y_pos, x_neg )      <<"\t"<< atan2        ( -y_nil, x_neg )      <<"\t"<< atan2        ( -y_neg, x_neg )      <<       std::endl;
	
	std::cout <<         atan2_neg_pos(  y_pos,  x_pos )      <<"\t"<< atan2_neg_pos(  y_nil,  x_pos )      <<"\t"<< atan2_neg_pos(  y_neg,  x_pos )      <<"\t"<< 
						 atan2_neg_pos(  y_pos,  x_nil )      <<"\t"<< atan2_neg_pos(  y_nil,  x_nil )      <<"\t"<< atan2_neg_pos(  y_neg,  x_nil )      <<"\t"<< 
						 atan2_neg_pos(  y_pos,  x_neg )      <<"\t"<< atan2_neg_pos(  y_nil,  x_neg )      <<"\t"<< atan2_neg_pos(  y_neg,  x_neg )      <<       std::endl;*/
	
	
	
	
	
	return 0;
}


double atan2_neg_neg( double y, double x ){
	// Guessed atan2(-y,-x). If y=0 and x>0, this returns pi whereas the inbuilt atan2 returns -pi. The correct result is pi.
	
	/*//----- Option 1 (both options should work) -----
	if     ( ( y>0 ) xor ( y==0 and x<=0 ) ){
		return atan2( y, x ) - pi;
	}
	else if( ( y<0 ) xor ( y==0 and x>=0 ) ){
		return atan2( y, x ) + pi;
	}*/
	
	//----- Option 2 (both options should work) -----
	if     ( ( x<=0 and y>=0 ) xor ( x>0 and y>0 ) ){
		return atan2( y, x ) - pi;
	}
	else if( ( x>=0 and y<=0 ) xor ( x<0 and y<0 ) ){
		return atan2( y, x ) + pi;
	}
	
	else {
		std::cout << "y=" << y << ", x=" << x <<" do not match the given cases." << std::endl;
		return 999;
	}
}

double atan2_pos_neg( double y, double x ){
	// Guessed atan2(y,-x).
	if     ( y>=0 ){
		return -atan2( y, x ) + pi;
	}
	else if( y<0  ){
		return -atan2( y, x ) - pi;
	}
	else {
		//std::cout << "y=" << y << ", x=" << x <<" do not match the given cases." << std::endl;
		return 999;
	}
}

double atan2_neg_pos( double y, double x ){
	// Guessed atan2(-y,x).
	if     ( y==0 and x<0 ){
		return -atan2( y, x ) - 2.0*pi;
	}
	return -atan2( y, x );
}