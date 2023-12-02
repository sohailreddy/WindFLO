/*

	Copyright 2019, Sohail R. Reddy
	email: sredd001@fiu.edu
	www.sohailreddy.com

*/


#include <cmath>
#include <vector>

#include "userDefined.hpp"



/*
		USER DEFINED WAKE MODEL
*/


std::vector<double> UserDefWakeModel::GetWakeDiameter(std::vector<double> position){


//	Example with Jensen's Model
	std::vector<double> wakeDiameter(2,0.0);

	double x = std::abs( T_position[0] - position[0] );
	wakeDiameter[0] = T_diameter * (1.0 + 2.0 * wakeExpansionCoeff[0] *
	                  (x / T_diameter ) ); 

	wakeDiameter[1] = wakeDiameter[0]; 
		
	return wakeDiameter;
}


std::vector<double> UserDefWakeModel::GetWakeVelocity(std::vector<double> position){

//	Example with Jensen's model

	std::vector<double> wakeVelocity(3,0.0);

	double x = std::abs( T_position[0] - position[0] );
		
	double numerator = 1.0 - std::sqrt(1.0 - T_Ct);
	double denominator = std::pow((1.0 + ( (2.0 * wakeExpansionCoeff[0] * x) /T_diameter)),2.0);

	double velocityDeficit = numerator / denominator;

	for(int i = 0; i < 3; i++){
		wakeVelocity[i] = (1.0 - velocityDeficit) * T_velocity[i];
	}
	
	return wakeVelocity;
}





void UserDefWakeMergeModel(std::vector<double> &velocitySum,
											 std::vector<double> ambientWind, 
											 std::vector<double> wakeVelocity,
											 double areaRatio,bool reverse){
											 
//	Example with quadratic superposition scheme											 
	if(reverse){
		for(int i = 0; i < 3; i++){
			velocitySum[i] = ambientWind[i] - std::sqrt(velocitySum[i]);
		}
		return;
	}								
				
	for(int i = 0; i < 3; i++){
		velocitySum[i] += areaRatio*(std::pow(ambientWind[i] - wakeVelocity[i],2.0));
	}
	
}


// Cost Model

double computeCost(double coe, double D, double H){

//	univariate power model
//	double ratedPower = 0.4368 * std::pow(D, 1.8803);

//	univariate polynomial model
//	double ratedPower = 0.0975 * std::pow(D, 2.0) + 16.133*D - 245.66;

//	bivariate power model
	double ratedPower = 1.510 * std::pow(D, 1.521) * std::pow(H, 0.076);
	
//	bivariate polynomial model	
//	double ratedPower = 0.122 * std::pow(D, 2.0) + 0.139 * std::pow(H, 2.0) 
//					  - 0.14*D*H + 22.037*D - 8.715*H - 169.204;

	return coe * ratedPower;
}

