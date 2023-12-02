/*

	Copyright 2019, Sohail R. Reddy
	email: sredd001@fiu.edu
	www.sohailreddy.com

*/


#ifndef WINFLO_HEADER
#define WINFLO_HEADER

#include <vector>
#include <string>
#include <cstring>
#include <cmath>
#include <iostream>


#include "ambientWindModels.hpp"
#include "wakeModels.hpp"
#include "turbine.hpp"
#include "costModels.hpp"

//#include "Input.hpp"

//	Constants

//#define PI 3.14159265359



//	std::string offBodyPoints;					// off body points at which to compute velocities if needed






class WindFLO
{

public:
	
	WindFLO(){};
	WindFLO(std::string inFileName);
	~WindFLO();
	
	std::string fileName;
	int nTurbines;
	bool batch;
	int nWindRose;

//	Results of WindFLO
	double landUsed;
	double farmEfficiency;
	double farmPower;
	double farmCost;
	double AEP;


//	Atmospheric parameters
	double rho;
	double turbulenceIntensity;


//	Ambient wind parameters	
	std::string windModel;		// wind model
	std::vector<double> modelVelocity;	// can be for constant, friction velocity (for log) or ref velocity (for power)
	double surfaceRoughness;	// for log profile
	double referenceHeight;		// for power profile

	
//	Wake modeling parameters
	std::string wakeModel;					// wake model
	std::string wakeMergeModel;				// wake merge model
	std::vector<double> wakeExpansionCoeff;	// wake expansion coefficient ... 2 for XA model
	

//	Numerical parameters
	int gaussOrder;
	int monteCarloPts;


//	Cost parameters
	double coe;

	
//	Turbine parameters
	std::vector<TurbineProps> turbines;
	

//	ConvexHull points
	std::vector<std::vector<double> > convexHull;
	

//	Rotation Matrix
	std::vector<std::vector<double> > rotationMat;	
	std::vector<std::vector<double> > rotationMat_inv;	
	
	
//	Functions
	void setup();
	void run();
	void Clean();
	void CleanAll();	
	void GetInput();
	void InitializeOutputs();
	void GetRunParameters();
	void GetAtmParameters();
	void GetAmbientWindParameters();
	void GetWakeParameters();
	void GetNumericalParameters();
	void GetCostParameters();
	void GetTurbine(int ith);
	void ComputeRotationMatrix();
	void PerformCoordinateTransformation(bool flag = true);
	void ComputeRank();
	void UnSortTurbines();
	void ComputeVelocities();
	void ComputeTurbineElevation();
	void CheckOrientationOfTurbines();
	void ComputeTurbinePower(double scale = 1.0);
	void ComputeLandUsage();
	void ComputeFarmCost();
	void ComputeFarmEfficiency();
	void ComputeFarmPower();
	void ComputeAEP(double probability = 1.0, double dVdTheta = 1.0);
	void GetWindRoseVelocity(int &ith, std::vector<double> &velocity, double &probability, double &dVdTheta);

	
	void WriteResults(std::string filename = "WindFLO.res");
	
};

extern "C" {

void ReadInput_(char fileName[]);
void GetRunParameters_(double runParams[]);
void GetAtmParameters_(double atmParams[]);
void GetAmbientWindParameters_(char modelName[], double ambientParams[]);
void GetWakeParameters_(char wakeModelName[], char wakeMergeModelName[], double wakeParams[]);
void GetNumericalParameters_(int numericalParams[]);
void GetCostParameters_(double costParams[]);
void GetTurbine_(int *ith, char name[], double turbineParams[], double Ctx[], double Cty[], double Cpx[], double Cpy[]);
void GetElevation_(double *x, double *y, double *z);
void GetWindRoseParameters_(int *ith, double windRoseParams[]);
void GetRotationMatrix_(double a[], double b[], double rotMat[]);
void MatInv3_(double a[]);
void ComputeArcLength_(double a[], double b[], double *length);
void CleanAll_();
void Clean_();
}

#endif