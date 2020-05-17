/*

	Copyright 2019, Sohail R. Reddy
	email: sredd001@fiu.edu
	www.sohailreddy.com

*/


#ifndef API_HEADER
#define API_HEADER

#include <vector>
#include <string>
#include <cstring>
#include <cmath>
#include <iostream>




class WindFLO;
class spline;

class TurbineAPI
{
	public : 
	 TurbineAPI(){};
	 ~TurbineAPI(){};

public:

	std::string *name;				// name of turbine
	std::vector<double> *position ;	// position of turbine
	std::vector<double> *orientation ;	// orientation of turbine	
	std::vector<double> *velocity ;	// average velocity seen by turbine
		
	double *height;					// height of turbine
	double *radius;					// radius of turbine
	double *area;					// rotor swept area
	double *ratedPower;				// rated power of turbine

	bool *fictitious;				// is it fictitious
	bool *yaw;						// can it yaw
	double *power;					// power generated
	
	spline *CpCurve;
	spline *CtCurve;

	void SetCpCurve(std::vector<double> &V, std::vector<double> &Cp);	
	void SetCtCurve(std::vector<double> &V, std::vector<double> &Ct);
	
};




class WindFLOAPI
{

public:
	
	WindFLOAPI(std::string inFileName);
	~WindFLOAPI();
	
	int *nTurbines;

//	Results of WindFLO
	double *landUsed;
	double *farmEfficiency;
	double *farmPower;
	double *farmCost;
	double *AEP;

//	Atmospheric parameters
	double *rho;

//	Ambient wind parameters	
	std::vector<double> *modelVelocity;	// can be for constant, friction velocity (for log) or ref velocity (for power)
	
//	Cost parameters
	double *coe;

//	Turbine parameters
	std::vector<TurbineAPI> turbines;
	
//	ConvexHull points
	std::vector<std::vector<double> > *convexHull;
	
//	Functions
	void run();
	void write(std::string filename = "WindFLO.res");

	WindFLO *windFLO;		
};




void FormWindFLOAPIConnections(WindFLOAPI *inWindFLO);
extern "C" {
void Python_WindFLO_API(char* infile, char* outfile, int n, double velocities[],
						double power[], double ratedPower[], double outputs[]);
}




#endif