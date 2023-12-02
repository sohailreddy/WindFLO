/*

	Copyright 2019, Sohail R. Reddy
	email: sredd001@fiu.edu
	www.sohailreddy.com

*/


#include "costModels.hpp"

double CostModel::ComputeCost(double D, double H){

//	univariate power model
//	double ratedPower = 0.4368 * std::pow(D, 1.8803);

//	univariate polynomial model
//	double ratedPower = 0.0975 * std::pow(D, 2.0) + 16.133*D - 245.66;

//	bivariate power model
//	double ratedPower = 1.510 * std::pow(D, 1.521) * std::pow(H, 0.076);
	
//	bivariate polynomial model	
//	double ratedPower = 0.122 * std::pow(D, 2.0) + 0.139 * std::pow(H, 2.0) 
//					  - 0.14*D*H + 22.037*D - 8.715*H - 169.204;

	double cost = computeCost(coe, D,H);

	return cost;
}
