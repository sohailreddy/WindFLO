/*

	Copyright 2019, Sohail R. Reddy
	email: sredd001@fiu.edu
	www.sohailreddy.com

*/


#ifndef AMBIENTWINDMODELS_HEADER
#define AMBIENTWINDMODELS_HEADER

#include <vector>
#include <string>
#include <algorithm>
#include <iostream>
#include <cmath>


class AmbientWindModel
{

	public:
	
	AmbientWindModel(){};
	AmbientWindModel(std::string inModel, std::vector<double> *inVelocity, 
					 double inSurfaceRoughness, double inReferenceHeight);

//	Variables
	std::string model;
	std::vector<double> *modelVelocity;
	double surfaceRoughness;
	double referenceHeight;
	int gaussOrder;
	double vonKarmanConstant = 0.40;
	

//	Functions
	std::vector<double> ComputeVelocity(double height);
	std::vector<double> ConstantProfile(double height = 0.0);	
	std::vector<double> LogLawProfile(double height);
	std::vector<double> PowerLawProfile(double height);	
	std::vector<double> DeavesHarrisProfile(double height);	



	std::vector<double> AverageVelocity(std::vector<double> center, double radius);

};









#endif