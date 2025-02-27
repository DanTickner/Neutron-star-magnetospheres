/*
cd OneDrive\PhD\Codes\Check vector expressions
g++ Check_cross_product_of_cross_product.cpp -o Check_cross_product_of_cross_product --std=c++11
Check_cross_product_of_cross_product

Check the expression in Cartesian coordinates
( A x B ) x B = ( Az*Bx*Bz - Ax*Bz^2 - Ax*By^2 - Ay*Bx*By ) i + ( Ax*By*Bx - Ay*Bx^2 - Ay*Bz^2 + Az*By*Bz ) j + ( Ay*Bz*By - Az*By^2 - Az*Bx^2 - Ax*Bz*Bx ) k
*/

#include <iostream>
#include <iomanip>
#include <math.h>
#include <vector>
#include <chrono>

#include "../Generate_Spherical_Harmonics.h"
#include "../Generate_Associated_Legendre_Functions.h"


int main(){
	
	//----- Define variables -----
	
	double Ax = 1.234;
	double Ay = 2.345;
	double Az = 4.321;
	double Bx = 5.321;
	double By = 6.324;
	double Bz = 0.432;
	
	std::vector<double> A { Ax, Ay, Az };
	std::vector<double> B { Bx, By, Bz };
	
	
	//----- Calculate -----
	std::vector<double> AxB = cross_product( A, B );
	std::vector<double> AxBxB_calculated = cross_product( AxB, B );
	
	
	cout_vector( AxBxB_calculated, "Calculated" );
	
	
	//----- Guessed -----
	
	//--- V1 ---
	//double AxBxB_guess_x = ( Az*Bx - Ax*Bz ) * Bz - ( Ax*By - Ay*Bx ) * By;
	//double AxBxB_guess_y = ( Ax*By - Ay*Bx ) * Bx - ( Ay*Bz - Az*By ) * Bz;	
	//double AxBxB_guess_z = ( Ay*Bz - Az*By ) * By - ( Az*Bx - Ax*Bz ) * Bx;
	
	//--- V2 ---
	double AxBxB_guess_x = ( Az*Bx*Bz - Ax*pow(Bz,2) - Ax*pow(By,2) + Ay*Bx*By );
	double AxBxB_guess_y = ( Ax*By*Bx - Ay*pow(Bx,2) - Ay*pow(Bz,2) + Az*By*Bz );
	double AxBxB_guess_z = ( Ay*Bz*By - Az*pow(By,2) - Az*pow(Bx,2) + Ax*Bz*Bx );
	
	std::vector<double> AxBxB_guess { AxBxB_guess_x, AxBxB_guess_y, AxBxB_guess_z };
	
	cout_vector( AxBxB_guess, "Guessed" );
	
	std::cout << "\ndone" << std::endl;
	
	
	return 0;
	
}