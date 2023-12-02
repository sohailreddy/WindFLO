/*

	Copyright 2019, Sohail R. Reddy
	email: sredd001@fiu.edu
	www.sohailreddy.com

*/

#include "IO.hpp"
#include "turbine.hpp"
#include "input.hpp"


TurbineProps::TurbineProps(char inName[], int inTurbineNum, double inPosition[], double inHeight, 
				  		   double inRadius, double inDiameter, bool inFictitious){
				  	  
		name = inName;
		turbineNum = inTurbineNum;
		height = inHeight;
		radius = inRadius;
		diameter = inDiameter;
		fictitious = inFictitious;
						  
		position.resize(3);
		velocity.resize(3);		
		
		position[0] = inPosition[0]; position[1] = inPosition[1]; position[2] = inPosition[2];
		velocity[0] = 0.0; velocity[1] = 0.0; velocity[2] = 0.0;
		
		area = PI * radius * radius;
		power = 0.0;
}


TurbineProps::TurbineProps(char inName[], double turbineParams[]){

	name = inName;

	turbineNum = turbineParams[0];
	
	position.resize(3);

	position[0] = turbineParams[1];
	position[1] = turbineParams[2];
	position[2] = turbineParams[3];	
	
	height = turbineParams[4];
	
	radius = turbineParams[5];
	diameter = turbineParams[6];
	fictitious = (int)turbineParams[7];
	ratedPower = turbineParams[8];
	power = 0.0;

	orientation.resize(3);
	orientation[0] = turbineParams[9];
	orientation[1] = turbineParams[10];
	orientation[2] = turbineParams[11];
	
	yaw = (int)turbineParams[12];

	area = PI * radius * radius;	

	velocity.resize(3);		
	velocity[0] = 0.0; velocity[1] = 0.0; velocity[2] = 0.0;	
}

/*	
TurbineProps::TurbineProps(const TurbineProps &inTurbine){
	
	name = inTurbine.name;
	turbineNum = inTurbine.turbineNum;
	height = inTurbine.height;
	radius = inTurbine.radius;
	diameter = inTurbine.diameter;
	area = inTurbine.area;
	rank = inTurbine.rank;	
	
	rpm = inTurbine.rpm;
	fictitious = inTurbine.fictitious;
	position = inTurbine.position;
	velocity = inTurbine.velocity;

	CpSpline = inTurbine.CpSpline;
	CtSpline = inTurbine.CtSpline;
	
	wakeModel = inTurbine.wakeModel;
		
	wakeMergeModel = inTurbine.wakeMergeModel;
	ambientWindModel = inTurbine.ambientWindModel;
	
	influencedBy = inTurbine.influencedBy;		
}
*/
	
void TurbineProps::setCpCurve(std::vector<double> &V, std::vector<double> &Cp){
	CpSpline.set_points(V,Cp);		
}
void TurbineProps::setCtCurve(std::vector<double> &V, std::vector<double> &Ct){
	CtSpline.set_points(V,Ct);				
}
	
	
void TurbineProps::SetWakeModel(std::string modelName, 
								std::vector<double> wakeExpansionCoeff,
								double turbulenceIntensity){

//	if(wakeModel != nullptr) delete wakeModel;
								
	wakeModel = wakeModel->getWakeModelPointer(modelName);
//	wakeModel = getWakeModelPointer(modelName);
	wakeModel->SetParent(this);
	wakeModel->SetParameter(wakeExpansionCoeff,turbulenceIntensity);
								
}

void TurbineProps::SetWakeMergeModel(std::string modelName){								
	wakeMergeModel = WakeMergeModel(this,modelName);
}


void TurbineProps::SetAmbientWindModel(std::string modelName, std::vector<double> *modelVelocity, 
										double surfaceRoughness, double referenceHeight) {
										
	ambientWindModel = AmbientWindModel(modelName, modelVelocity, 
				  				   surfaceRoughness, referenceHeight);				  				   
}

	
void TurbineProps::RunCheck(){

	if(this != wakeModel->parentTurbine){
		wakeModel->SetParent(this);
	}

	if(this != wakeMergeModel.parentTurbine){
		wakeMergeModel.SetParent(this);
	}

}	
	
std::vector<double> TurbineProps::GetHubCoordinates(){

	std::vector<double> hubCoordinates(3,0.0);

	hubCoordinates = this->position;
	hubCoordinates[2] += this->height; 
	
	return hubCoordinates;
}


std::vector<double> TurbineProps::ComputeAverageAmbientVelocity(){

	std::vector<double> hubCenter = this->GetHubCoordinates();
	hubCenter[2] -= this->position[2];
	
	std::vector<double> ambientWind = ambientWindModel.AverageVelocity(hubCenter, this->radius);

	return ambientWind;
}

double TurbineProps::computeCt(double uInf){

	double Unorm = 0.0;
	if(uInf == 0.0){
		Unorm = std::sqrt(velocity[0]*velocity[0] + velocity[1]*velocity[1] + velocity[2]*velocity[2]);
	}
	else{
		Unorm = uInf;
	}

	return this->CtSpline(Unorm);
}


double TurbineProps::computeCp(double uInf){

	double Unorm = 0.0;
	if(uInf == 0.0){
		Unorm = std::sqrt(velocity[0]*velocity[0] + velocity[1]*velocity[1] + velocity[2]*velocity[2]);
	}
	else{
		Unorm = uInf;
	}
	
	return this->CpSpline(Unorm);
}

double TurbineProps::computeRatedPower(double rho){


	double tmpRatedPower = 0.0;
	double delta = (CpSpline.right - CpSpline.left)/40.0;
	
	for(int i = 0; i <= 40; i++){
		double tmpPower = this->computePower(rho,  CpSpline.left + (delta * i) );
		tmpRatedPower = std::max(tmpRatedPower, tmpPower);
	}
	return tmpRatedPower;
}

double TurbineProps::computePower(double rho, double uInf){
	double Cp = computeCp(uInf);

	double Unorm = 0.0;
	if(uInf == 0.0){
		Unorm = std::sqrt(velocity[0]*velocity[0] + velocity[1]*velocity[1] + velocity[2]*velocity[2]);
	}
	else{
		Unorm = uInf;
	}

	double tmpPower = 0.5 * Cp * rho * std::pow(Unorm,3.0) * PI * std::pow(this->radius,2.0);
	
	return tmpPower;
}

