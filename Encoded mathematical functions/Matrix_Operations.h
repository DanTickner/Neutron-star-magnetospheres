/*
Matrix_Operations.h

Matrix calculator


*/


#include <vector>
#include <iostream>
#include <iomanip>


//----- Determinants -----
double det_3x3( std::vector< std::vector<double> > A ){
	return A[0][0] * ( A[1][1]*A[2][2] - A[1][2]*A[2][1] ) + A[0][1] * ( A[1][2]*A[2][0] - A[1][0]*A[2][2] ) + A[0][2] * ( A[1][0]*A[2][1] - A[1][1]*A[2][0] );
}



//----- Inverse matrices -----

std::vector< std::vector<double> > inverse_3x3( std::vector< std::vector<double> > A ){
	double det = det_3x3( A );
	
	std::vector< std::vector<double> > inverse ( 3, std::vector<double> ( 3 ) );
	
	if( det == 0 ){
		std::cout << "Determinant is zero so inverse does not exist. Returning the zero matrix." << std::endl;
	}
	
	else {
		inverse[0] = { ( A[1][1]*A[2][2]-A[1][2]*A[2][1] )/det, ( A[0][2]*A[2][1]-A[0][1]*A[2][2] )/det, ( A[0][1]*A[1][2]-A[0][2]*A[1][1] )/det };
		inverse[1] = { ( A[1][2]*A[2][0]-A[1][0]*A[2][2] )/det, ( A[0][0]*A[2][2]-A[0][2]*A[2][0] )/det, ( A[0][2]*A[1][0]-A[0][0]*A[1][2] )/det };
		inverse[2] = { ( A[1][0]*A[2][1]-A[1][1]*A[2][0] )/det, ( A[0][1]*A[2][0]-A[0][0]*A[2][1] )/det, ( A[0][0]*A[1][1]-A[0][1]*A[1][0] )/det };
	}
	
	return inverse;
}




//----- Cout matrices -----
void cout_matrix_3x3( std::vector< std::vector<double> > A, int w=15 ){
	std::cout << std::endl;
	std::cout <<std::left<<std::setw(w)<< A[0][0] <<std::left<<std::setw(w)<< A[0][1] <<std::left<<std::setw(w)<< A[0][2] << std::endl;
	std::cout <<std::left<<std::setw(w)<< A[1][0] <<std::left<<std::setw(w)<< A[1][1] <<std::left<<std::setw(w)<< A[1][2] << std::endl;
	std::cout <<std::left<<std::setw(w)<< A[2][0] <<std::left<<std::setw(w)<< A[2][1] <<std::left<<std::setw(w)<< A[2][2] << std::endl;
	std::cout << std::endl;
}