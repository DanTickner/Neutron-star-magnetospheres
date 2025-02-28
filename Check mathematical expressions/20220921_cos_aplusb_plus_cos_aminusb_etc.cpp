// cd OneDrive\PhD\Codes\Test mathematical expressions
// g++ 20220921_cos_aplusb_plus_cos_aminusb_etc.cpp -o 20220921_cos_aplusb_plus_cos_aminusb_etc

/*
cos(a+b) + cos(a+b) = 2 cos(a) cos(b)
cos(a-b) - cos(a+b) = 2 sin(a) sin(b)
sin(a+b) + sin(a-b) = 2 sin(a) cos(b)
sin(a+b) - sin(a-b) = 2 cos(a) sin(b)
sin(a+b) + cos(a+b) = sin(a) [ cos(b) - sin(b) ] + cos(a) [ cos(b) + sin(b) ]
sin(a+b) - cos(a+b) = sin(a) [ cos(b) + sin(b) ] + cos(a) [ sin(b) - cos(b) ]
etc
*/

#include <iostream>
#include <math.h>

int main(){
	
	//----- Define variables -----
	double a = 0.3;
	double b = 0.5;
	
	
	//----- Output results for only the half-integers -----
	std::cout << cos(a+b) + cos(a-b) << "\t" << 2 * cos(a) * cos(b) << std::endl;
	std::cout << cos(a-b) - cos(a+b) << "\t" << 2 * sin(a) * sin(b) << std::endl;
	std::cout << sin(a+b) + sin(a-b) << "\t" << 2 * sin(a) * cos(b) << std::endl;
	std::cout << sin(a+b) - sin(a-b) << "\t" << 2 * cos(a) * sin(b) << std::endl;
	std::cout << sin(a+b) + cos(a+b) << "\t" << sin(a) * ( cos(b) - sin(b) ) + cos(a) * ( cos(b) + sin(b) ) << std::endl;
	std::cout << sin(a+b) - cos(a+b) << "\t" << sin(a) * ( cos(b) + sin(b) ) + cos(a) * ( sin(b) - cos(b) ) << std::endl;
	
	
	
	return 0;
}