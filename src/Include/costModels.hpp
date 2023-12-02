/*

	Copyright 2019, Sohail R. Reddy
	email: sredd001@fiu.edu
	www.sohailreddy.com

*/

#ifndef COSTMODELS_HEADER
#define COSTMODELS_HEADER

#include <vector>
#include <string>
#include <cmath>
#include <iostream>


class CostModel{
	public:
		CostModel(double inCOE){ coe = inCOE; }
		
		double coe;	
		virtual double ComputeCost(double D, double H);
};



double computeCost(double coe, double D, double H);

#endif