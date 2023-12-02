/*

	Copyright 2019, Sohail R. Reddy
	email: sredd001@fiu.edu
	www.sohailreddy.com

*/
#include "API.hpp"
#include "WindFLO.hpp"


WindFLOAPI::WindFLOAPI(std::string inFileName){

	this->windFLO = nullptr;	
	windFLO = new WindFLO(inFileName);
	this->windFLO->setup();			
	FormWindFLOAPIConnections(this);
}



WindFLOAPI::~WindFLOAPI(){
	Clean_();
	if(this->windFLO != nullptr) delete this->windFLO;	
}

void WindFLOAPI::run(){
	windFLO->run();	
}

void WindFLOAPI::write(std::string filename){
	windFLO->WriteResults(filename);
}


void FormWindFLOAPIConnections(WindFLOAPI *inWindFLO){

	inWindFLO->nTurbines = &inWindFLO->windFLO->nTurbines;

//	Results of WindFLO
	inWindFLO->landUsed = &inWindFLO->windFLO->landUsed;
	inWindFLO->farmEfficiency = &inWindFLO->windFLO->farmEfficiency;
	inWindFLO->farmPower = &inWindFLO->windFLO->farmPower;
	inWindFLO->farmCost = &inWindFLO->windFLO->farmCost;
	inWindFLO->AEP = &inWindFLO->windFLO->AEP;

//	Atmospheric parameters
	inWindFLO->rho = &inWindFLO->windFLO->rho;

//	Ambient wind parameters	
	inWindFLO->modelVelocity = &inWindFLO->windFLO->modelVelocity;
		
//	Cost parameters
	inWindFLO->coe = &inWindFLO->windFLO->coe;

//	ConvexHull points
	inWindFLO->convexHull = &inWindFLO->windFLO->convexHull;

//	Turbine parameters
	inWindFLO->turbines.resize(*inWindFLO->nTurbines);

	for(int i = 0; i < *inWindFLO->nTurbines; i++){
	
		inWindFLO->turbines[i].name = &inWindFLO->windFLO->turbines[i].name;	
		inWindFLO->turbines[i].position = &inWindFLO->windFLO->turbines[i].position;	
		inWindFLO->turbines[i].orientation = &inWindFLO->windFLO->turbines[i].orientation;	
		inWindFLO->turbines[i].velocity = &inWindFLO->windFLO->turbines[i].velocity;	
		
		inWindFLO->turbines[i].height = &inWindFLO->windFLO->turbines[i].height;	
		inWindFLO->turbines[i].radius = &inWindFLO->windFLO->turbines[i].radius;	
		inWindFLO->turbines[i].area = &inWindFLO->windFLO->turbines[i].area;	
		inWindFLO->turbines[i].ratedPower = &inWindFLO->windFLO->turbines[i].ratedPower;	

		inWindFLO->turbines[i].fictitious = &inWindFLO->windFLO->turbines[i].fictitious;	
		inWindFLO->turbines[i].yaw = &inWindFLO->windFLO->turbines[i].yaw;	
		inWindFLO->turbines[i].power = &inWindFLO->windFLO->turbines[i].power;	

		inWindFLO->turbines[i].CpCurve = &inWindFLO->windFLO->turbines[i].CpSpline;
		inWindFLO->turbines[i].CtCurve = &inWindFLO->windFLO->turbines[i].CtSpline;
	
	}	

}


void TurbineAPI::SetCpCurve(std::vector<double> &V, std::vector<double> &Cp){
	CpCurve->set_points(V,Cp);		
}	
void TurbineAPI::SetCtCurve(std::vector<double> &V, std::vector<double> &Ct){
	CtCurve->set_points(V,Ct);
}




void Python_WindFLO_API(char infile[], char outfile[], int n, double velocities[], 
						double power[], double ratedPower[], double outputs[]){

//	std::string inFile = "windFLO.inp";		
	std::string inFile = infile;
	std::string outFile = outfile;
				
	WindFLOAPI PyWindFLO(inFile);	
	PyWindFLO.run();	
	if(outFile != "") {
		PyWindFLO.write(outFile);
	}
	
	int k = 0;
	for (int i = 0; i < n; i++){
		for (int j = 0; j < 3; j++, k++){	
			velocities[k] = (*PyWindFLO.turbines[i].velocity)[j]; 
		}
		power[i] = *PyWindFLO.turbines[i].power;
		ratedPower[i] = *PyWindFLO.turbines[i].ratedPower;
	}
	
	//	Results of WindFLO
	outputs[0] = *PyWindFLO.AEP;
	outputs[1] = *PyWindFLO.farmPower;
	outputs[2] = *PyWindFLO.farmEfficiency;
	outputs[3] = *PyWindFLO.farmCost;	
	outputs[4] = *PyWindFLO.landUsed;
	
	return;
}



