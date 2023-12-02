/*

	Copyright 2019, Sohail R. Reddy
	email: sredd001@fiu.edu
	www.sohailreddy.com

*/


#ifndef TURBINE_HEADER
#define TURBINE_HEADER

#include <vector>
#include <string>
#include <algorithm>
#include <iostream>
#include <cmath>

#include "CubicSpline.hpp"
#include "ambientWindModels.hpp"
#include "wakeModels.hpp"
#include "costModels.hpp"


//	Constants

//#define PI 3.14159265359

class TurbineProps
{
	public : 
	TurbineProps(char inName[], int inTurbineNum, double inPosition[], double inHeight, 
				  double inRadius, double inDiameter, bool inFictitious);
				  
	TurbineProps(char inName[], double turbineParams[]);

				  
	 TurbineProps(){};
	 
	 ~TurbineProps(){
	 };

public:

	std::string name;				// name of turbine
	int turbineNum;					// turbine number
	std::vector<double> position ;	// position of turbine
	std::vector<double> orientation ;	// orientation of turbine	
	std::vector<double> velocity ;	// average velocity seen by turbine
		
	double height;					// height of turbine
	double radius;					// radius of turbine
	double diameter;				// diameter of turbine
	double area;					// rotor swept area
	double ratedPower;				// rated power of turbine

	spline CpSpline;			// Cubic spline interpolation for Cp Curve
	spline CtSpline;			// Cubic spline interpolation for Cp Curve
		
	WakeModel	*wakeModel;			// Wake model
	WakeMergeModel	wakeMergeModel;	// Wake merge model	
	AmbientWindModel ambientWindModel;	// Ambient wind model
	
	bool fictitious;				// is it fictitious
	bool yaw;						// can it yaw
	int rank;						// rank of the turbine	
	
	std::vector<TurbineProps*> influencedBy;	// influencing turbine	

	

	double power;					// power generated


//	double rpm;						// rpm	
//	double reynoldsNumber;			// Reynolds number	
//	double inducFactor;				// induction factor
//	double TSR;						// tip-speed-ratio
	


//	Functions
	void setCpCurve(std::vector<double> &V, std::vector<double> &Cp);
	void setCtCurve(std::vector<double> &V, std::vector<double> &Ct);

	void SetAmbientWindModel(std::string modelName, std::vector<double> *modelVelocity, 
							double surfaceRoughness, double referenceHeight);
														
	void SetWakeModel(std::string modelName, std::vector<double> wakeExpansionCoeff, double turbulenceIntensity);
	void SetWakeMergeModel(std::string modelName);	
	std::vector<double> GetHubCoordinates();
	std::vector<double> ComputeAverageAmbientVelocity();
	double computeCt(double uInf = 0.0);
	double computeCp(double uInf = 0.0);
	double computeRatedPower(double rho);
	double computePower(double rho,double uInf = 0.0);
	void RunCheck();	
	
	
	
	
	double dummy;
};






#endif