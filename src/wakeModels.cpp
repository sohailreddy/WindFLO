/*

	Copyright 2019, Sohail R. Reddy
	email: sredd001@fiu.edu
	www.sohailreddy.com

*/


#include <algorithm>
#include <math.h>

#include "IO.hpp"
#include "turbine.hpp"
#include "wakeModels.hpp"
#include "sobol.hpp"
#include "input.hpp"
#include "userDefined.hpp"

std::vector<double> GetArcLengthVector(std::vector<double> A, std::vector<double> B);


bool ichar_equals(char a, char b)
{
    return std::tolower(static_cast<unsigned char>(a)) ==
           std::tolower(static_cast<unsigned char>(b));
}


bool iequals(const std::string& a, const std::string& b)
{
    return a.size() == b.size() &&
           std::equal(a.begin(), a.end(), b.begin(), ichar_equals);
}


void CheckNanVelocity(std::vector<double> &inVelocity){

	for(int i = 0; i < inVelocity.size(); i++){
		if(inVelocity[i] != inVelocity[i]) inVelocity[i] = 0.0;
	}
	return;
}


double ComputeOverLappingArea(std::vector<double> centerOne, std::vector<double> centerTwo,
							   double radiusOne, double radiusTwo) {

	double distance = std::sqrt( std::pow(centerTwo[1] - centerOne[1],2.0) + std::pow(centerTwo[2] - centerOne[2],2.0) );
	 	 
	if(distance > (radiusOne + radiusTwo)){
		//two circles dont intersect
		return 0.0;
	}
	else if (distance <= std::abs( radiusOne - radiusTwo ) ){
		// complete overlap ... return area of smaller circle
		return std::min(PI*radiusOne*radiusOne, PI*radiusTwo*radiusTwo);
	}

	double a_distanceCenterFirst = ((radiusOne*radiusOne) - (radiusTwo*radiusTwo) + (distance*distance)) / (2.0 * distance);
	double b_distanceCenterSecond = distance - a_distanceCenterFirst;
	double h_height = std::sqrt((radiusOne*radiusOne) - (a_distanceCenterFirst*a_distanceCenterFirst));

	double alpha = fmod((std::atan2(h_height, a_distanceCenterFirst) * 2.0 + 2.0 * PI) , (2.0 * PI));
	double beta = fmod((std::atan2(h_height, b_distanceCenterSecond) * 2.0 + 2.0 * PI) , (2.0 * PI));
	
	double areaOne = (radiusOne*radiusOne) / 2.0 * (alpha - std::sin(alpha));
	double areaTwo = (radiusTwo*radiusTwo) / 2.0 * (beta - std::sin(beta));
	
	return areaOne + areaTwo;
}

WakeModel *getWakeModelPointer(std::string modelName){

	if (iequals(modelName, "Jensen")) {
		return new JensenModel();
	}
	else if (iequals(modelName, "Frandsen")){
		return new FrandsenModel();
	}
	else if (iequals(modelName, "Larsen")){
		return new LarsenModel();
	}
	else if (iequals(modelName, "Ishihara")){
		return new IshiharaModel();
	}
	else if (iequals(modelName, "BP")){
		return new BPModel();
	}
	else if (iequals(modelName, "XA")){
		return new XAModel();
	}
	else if (iequals(modelName, "UserDef")){
		return new UserDefModel();
	}	
	else{
		std::cout << "Incorrect wake model" << std::endl;
		exit(10);
	}
	
}

WakeModel *WakeModel::getWakeModelPointer(std::string modelName){

	if (iequals(modelName, "Jensen")) {
		return new JensenModel();
	}
	else if (iequals(modelName, "Frandsen")){
		return new FrandsenModel();
	}
	else if (iequals(modelName, "Larsen")){
		return new LarsenModel();
	}
	else if (iequals(modelName, "Ishihara")){
		return new IshiharaModel();
	}
	else if (iequals(modelName, "BP")){
		return new BPModel();
	}
	else if (iequals(modelName, "XA")){
		return new XAModel();
	}
	else if (iequals(modelName, "UserDef")){
		return new UserDefModel();
	}		
	else{
		std::cout << "Incorrect wake model" << std::endl;
		exit(10);
	}
	
}


double ComputeBeta(double Ct){
	double beta = 0.0;
	beta = (1.0 + std::sqrt(1.0 - Ct)) / (2.0 * std::sqrt(1.0 - Ct));
	return beta;
}

WakeModel::WakeModel(){


	modelName = "WakeModel";
	scaleFactor = 1.0;
	turbulenceIntensity = 0.0;
	wakeExpansionCoefficient.resize(3);
	
	wakeExpansionCoefficient[0] = 0.0;
	wakeExpansionCoefficient[1] = 0.0;
	wakeExpansionCoefficient[2] = 0.0;		

}

/*
WakeModel::WakeModel(const WakeModel *inWakeModel){
	turbulenceIntensity = inWakeModel->turbulenceIntensity;
	wakeExpansionCoefficient = inWakeModel->wakeExpansionCoefficient;	
	parentTurbine = inWakeModel->parentTurbine;	
}
*/


void WakeModel::SetParameter(std::vector<double> in_WakeExpansionCoefficient, double in_TurbulenceIntensity){
	wakeExpansionCoefficient = in_WakeExpansionCoefficient;
	turbulenceIntensity = in_TurbulenceIntensity;
	scaleFactor = 1.0;
}

void WakeModel::SetParent(TurbineProps *inTurbine){
	parentTurbine = inTurbine;
}


std::vector<double> WakeModel::GetWakeDiameter(std::vector<double> position){

}

std::vector<double> WakeModel::GetWakeVelocity(std::vector<double> position){

}

bool WakeModel::CheckIfInsideWakeRegion(std::vector<double> position){


	std::vector<double> centerOfWake = parentTurbine->GetHubCoordinates();
	centerOfWake[0] = position[0];
	
	std::vector<double> wakeDiameter = this->GetWakeDiameter(position); 
	wakeDiameter[0] *= 0.5 ;
	wakeDiameter[1] *= 0.5 ;
	
	double sum = std::pow( ( position[1] - centerOfWake[1] ) / wakeDiameter[0] , 2.0 ) +
			  	 std::pow( ( position[2] - centerOfWake[2] ) / wakeDiameter[1] , 2.0 ) ;
			  	 
	return sum <= 1.0 ? 1 : 0;		  	 
}

double WakeModel::ComputeWakeOverLapArea(std::vector<double> position, double radius){

	std::vector<double> hubLocation = parentTurbine->GetHubCoordinates();
	hubLocation[0] = position[0];
	
	std::vector<double> wakeDiameter = GetWakeDiameter(position);					
	double overlapArea = ComputeOverLappingArea(hubLocation, position,
									            wakeDiameter[0]*0.5, radius);									            
	return overlapArea;
}


std::vector<double> WakeModel::AverageWakeVelocity(std::vector<double> center, double radius) {

	std::vector<double> integratedVelocity(3,0.0);

	if(radius <= 0.0) {
		integratedVelocity = this->GetWakeVelocity(center);
		return integratedVelocity;
	}


	int dim_num = 2;
	long long int seed = 2;
	double rsobol[2];

	std::vector<double> pointOnDisk(3,0.0);
	std::vector<double> wakeVelocity(3,0.0);	
	std::vector<double> ambientWind = parentTurbine->ComputeAverageAmbientVelocity();
	double r, theta;
	
	pointOnDisk[0] = center[0];
		
	int nPointsInside = 0;		
	for(int i = 0; i < monteCarloPts; i++){

		i8_sobol( dim_num, &seed, rsobol );	
	
		r = rsobol[1] * radius;
		theta = 2.0 * PI * rsobol[0];
	
		pointOnDisk[1] = r * std::cos(theta) + center[1];
		pointOnDisk[2] = r * std::sin(theta) + center[2]; 
		
		
		if(this->CheckIfInsideWakeRegion(pointOnDisk)) {
	
			wakeVelocity = this->GetWakeVelocity(pointOnDisk);

			integratedVelocity[0] += wakeVelocity[0];
			integratedVelocity[1] += wakeVelocity[1];
			integratedVelocity[2] += wakeVelocity[2];
			
			nPointsInside += 1;
//			nPointsInside = monteCarloPts;
			
		}
		
	}	

	if(nPointsInside == 0) nPointsInside = 1;

	integratedVelocity[0] *= (1.0 / nPointsInside);
	integratedVelocity[1] *= (1.0 / nPointsInside);
	integratedVelocity[2] *= (1.0 / nPointsInside);

	for(int i = 0; i < 3; i++){
		if(integratedVelocity[i] == 0.0){
			integratedVelocity[i] = ambientWind[i];
		}
	}


	return integratedVelocity;
}

/*
		JENSEN WAKE MODEL
*/


std::vector<double> JensenModel::GetWakeDiameter(std::vector<double> position){

	std::vector<double> wakeDiameter(2,0.0);

	double x = std::abs( parentTurbine->position[0] - position[0] );
	wakeDiameter[0] = parentTurbine->diameter * (1.0 + 2.0 * wakeExpansionCoefficient[0] *
	                  (x / parentTurbine->diameter ) ); 

	wakeDiameter[1] = wakeDiameter[0]; 
	
	return wakeDiameter;
}



std::vector<double> JensenModel::GetWakeVelocity(std::vector<double> position){

	std::vector<double> wakeVelocity(3,0.0);

	double x = std::abs( parentTurbine->position[0] - position[0] );
	double Ct = parentTurbine->computeCt();
	
	
	double numerator = 1.0 - std::sqrt(1.0 - Ct);
	double denominator = std::pow((1.0 + ( (2.0 * wakeExpansionCoefficient[0] * x) /parentTurbine->diameter)),2.0);

	double velocityDeficit = numerator / denominator;

	for(int i = 0; i < 3; i++){
		wakeVelocity[i] = (1.0 - velocityDeficit) * parentTurbine->velocity[i];
	}
	
	return wakeVelocity;
}



/*
		FRANDSEN WAKE MODEL
*/


std::vector<double> FrandsenModel::GetWakeDiameter(std::vector<double> position){

	std::vector<double> wakeDiameter(2,0.0);

	
	double x = std::abs( parentTurbine->position[0] - position[0] );
	double Ct = parentTurbine->computeCt();
	double beta = ComputeBeta(Ct);

	double k = 2.0;

	wakeExpansionCoefficient[0] = std::pow(beta,0.5*k) * (  std::pow( 1.0 + 2.0 * 0.05 * (x/parentTurbine->diameter) ,k) - 1.0)
								  *( parentTurbine->diameter/x);

	wakeExpansionCoefficient[1] = wakeExpansionCoefficient[0];

	wakeDiameter[0] = parentTurbine->diameter * std::pow( std::pow(beta,0.5*k) + wakeExpansionCoefficient[0] 
					* ( x /parentTurbine->diameter),1.0/k ); 
					
	wakeDiameter[1] = wakeDiameter[0]; 
	
	return wakeDiameter;
}


std::vector<double> FrandsenModel::GetWakeVelocity(std::vector<double> position){

	std::vector<double> wakeVelocity(3,0.0);

	double x = std::abs( parentTurbine->position[0] - position[0] );
	double Ct = parentTurbine->computeCt();
	
	std::vector<double> wakeDiameter = this->GetWakeDiameter(position); 
	
	double wakeArea = PI * wakeDiameter[0] * wakeDiameter[0] / 4.0;

	double velocityDeficit = 0.5 * (1.0 - std::sqrt(1.0 - 2.0 * ( parentTurbine->area / wakeArea) * Ct));
		
	for(int i = 0; i < 3; i++){
		wakeVelocity[i] = (1.0 - velocityDeficit) * parentTurbine->velocity[i];
	}
	
	return wakeVelocity;
}





/*
		LARSEN WAKE MODEL
*/


std::vector<double> LarsenModel::GetWakeDiameter(std::vector<double> position){

	std::vector<double> wakeDiameter(2,0.0);
	
	double x = std::abs( parentTurbine->position[0] - position[0] );
	double Ct = parentTurbine->computeCt();
	double beta = ComputeBeta(Ct);

	double a1 = 0.435449861;
	double a2 = 0.797853685;
	double a3 = -0.124807893;
	double a4 = 0.136821858;
	double b1 = 15.6298;


	double R96D = a1 * std::exp( a2 * Ct*Ct + a3*Ct + a4 ) * (b1 * turbulenceIntensity + 1.0) *  parentTurbine->diameter;  

	double x0 = (9.6*parentTurbine->diameter) / ( std::pow( (2.0*R96D) / (std::sqrt(beta)*parentTurbine->diameter),3.0) - 1.0 );

	double c1 = std::pow(105.0 / (2*PI), -0.5) * std::pow(0.5 * std::sqrt(beta) * parentTurbine->diameter, 2.5)
				* std::pow( Ct * parentTurbine->area * x0, -5.0/6.0 );


	wakeDiameter[0] = 2.0 * std::pow( (105.0 * c1*c1)/(2*PI),0.2) * std::pow( Ct * parentTurbine->area * (x + x0), 1.0/3.0 ); 				
	wakeDiameter[1] = wakeDiameter[0]; 
		
	return wakeDiameter;
}



std::vector<double> LarsenModel::GetWakeVelocity(std::vector<double> position){

	std::vector<double> wakeVelocity(3,0.0);

	double Ct = parentTurbine->computeCt();
	double beta = ComputeBeta(Ct);

	double a1 = 0.435449861;
	double a2 = 0.797853685;
	double a3 = -0.124807893;
	double a4 = 0.136821858;
	double b1 = 15.6298;

	std::vector<double> hubCoordinates = parentTurbine->GetHubCoordinates();

	double x = std::abs( hubCoordinates[0] - position[0] );
	double y = std::abs( hubCoordinates[1] - position[1] );
	double z = std::abs( hubCoordinates[2] - position[2] );	
	double r = std::sqrt( y*y + z*z );

	double R96D = a1 * std::exp( a2 * Ct*Ct + a3*Ct + a4 ) * (b1 * turbulenceIntensity + 1.0) *  parentTurbine->diameter;  

	double x0 = (9.6*parentTurbine->diameter) / ( std::pow( (2.0*R96D) / (std::sqrt(beta)*parentTurbine->diameter),3.0) - 1.0 );

	double c1 = std::pow(105.0 / (2*PI), -0.5) * std::pow(0.5 * std::sqrt(beta) * parentTurbine->diameter, 2.5)
				* std::pow( Ct * parentTurbine->area * x0, -5.0/6.0 );


	double velocityDeficit = (-1.0/9.0) * std::pow( Ct * parentTurbine->area * std::pow(x + x0,-2.0)  ,1.0/3.0)
							* std::pow( std::pow(r,1.5) * std::pow(3.0 * c1*c1 * Ct * parentTurbine->area * (x+x0),-0.5)
							- std::pow( 35.0 / (2.0*PI), 3.0/10.0) * std::pow(3.0 * c1 * c1,-0.2 ) ,2.0);
		
	for(int i = 0; i < 3; i++){
		wakeVelocity[i] = (1.0 - velocityDeficit) * parentTurbine->velocity[i];
	}

	return wakeVelocity;
}


	

/*
		ISHIHARA WAKE MODEL
*/



std::vector<double> IshiharaModel::GetWakeDiameter(std::vector<double> position){

	std::vector<double> wakeDiameter(2,0.0);
	
	double x = std::abs( parentTurbine->position[0] - position[0] );
	double Ct = parentTurbine->computeCt();

	double k = 0.11 * std::pow(Ct, 1.07) * std::pow(turbulenceIntensity, 0.2);
	double epsilon = 0.23 * std::pow(Ct, -0.25) * std::pow(turbulenceIntensity, 0.17);

	scaleFactor = 2.0;
	wakeDiameter[0] =  scaleFactor * parentTurbine->diameter * (k  * (x/parentTurbine->diameter) + epsilon );
	wakeDiameter[1] = wakeDiameter[0]; 
		
	return wakeDiameter;
}





std::vector<double> IshiharaModel::GetWakeVelocity(std::vector<double> position){

	std::vector<double> wakeVelocity(3,0.0);

	std::vector<double> hubCoordinates = parentTurbine->GetHubCoordinates();

	double x = std::abs( hubCoordinates[0] - position[0] );
	double y = std::abs( hubCoordinates[1] - position[1] );
	double z = std::abs( hubCoordinates[2] - position[2] );

	double r = std::sqrt( y*y + z*z );
	double Ct = parentTurbine->computeCt();

	double a = 0.93 * std::pow(Ct, -0.75) * std::pow(turbulenceIntensity, 0.17);
	double b = 0.42 * std::pow(Ct, 0.6) * std::pow(turbulenceIntensity, 0.2);
	double c = 0.15 * std::pow(Ct, -0.25) * std::pow(turbulenceIntensity, -0.7);

	double wakeDiameter = this->GetWakeDiameter(position)[0] / scaleFactor;
	
	double velocityDeficit = (1.0 / std::pow(a + b*(x/parentTurbine->diameter) 
						+ c*std::pow(1.0 + x/parentTurbine->diameter,-2.0),2.0 ) ) 		
	 				    * std::exp( -0.5* (r*r) / (wakeDiameter * wakeDiameter) );
		
	for(int i = 0; i < 3; i++){	
		wakeVelocity[i] = (1.0 - velocityDeficit) * parentTurbine->velocity[i];		
	}

	return wakeVelocity;
}


/*

std::vector<double> IshiharaModel::GetWakeDiameter(std::vector<double> position){

	std::vector<double> wakeDiameter(2,0.0);
	
	double x = std::abs( parentTurbine->position[0] - position[0] );
	double Ct = parentTurbine->computeCt();

	double k1 = 0.27;
	double k2 = 6.0;
	double k3 = 0.004;

	double Iw = (k3 * Ct) / std::max(turbulenceIntensity,0.03) 
				* (1.0 - std::exp ( -4.0 * std::pow( x / (10.0 * parentTurbine->diameter) ,2.0) ) );

	double p = k2 * (turbulenceIntensity + Iw);

	wakeDiameter[0] = (k1/0.833) * std::pow(Ct, 0.25) * std::pow(parentTurbine->diameter,1.0-0.5*p)
					  * std::pow(x, 0.5 * p);
	wakeDiameter[1] = wakeDiameter[0]; 
		
	return wakeDiameter;
}





std::vector<double> IshiharaModel::GetWakeVelocity(std::vector<double> position){

	std::vector<double> wakeVelocity(3,0.0);

	std::vector<double> hubCoordinates = parentTurbine->GetHubCoordinates();

	double x = std::abs( hubCoordinates[0] - position[0] );
	double y = std::abs( hubCoordinates[1] - position[1] );
	double z = std::abs( hubCoordinates[2] - position[2] );

	double r = std::sqrt( y*y + z*z );
	double Ct = parentTurbine->computeCt();

	double k1 = 0.27;
	double k2 = 6.0;
	double k3 = 0.004;

	double Iw = (k3 * Ct) / std::max(turbulenceIntensity,0.03) 
				* (1.0 - std::exp ( -4.0 * std::pow( x / (10.0 * parentTurbine->diameter) ,2.0) ) );

	double p = k2 * (turbulenceIntensity + Iw);

	std::vector<double> wakeDiameter = this->GetWakeDiameter(position);
	std::vector<double> ambientWind = parentTurbine->ComputeAverageAmbientVelocity();
	
	for(int i = 0; i < 3; i++){	
	
		wakeVelocity[i] = (ambientWind[i] * std::sqrt(Ct)/32.0) * std::pow(1.666/k1,2.0)
						  * std::pow(x/parentTurbine->diameter,-1.0*p)
						  * std::exp( - (r*r) / (wakeDiameter[0] * wakeDiameter[0]) );
		
	}

	return wakeVelocity;
}

*/

/*
		BASTANKHAH WAKE MODEL
*/


std::vector<double> BPModel::GetWakeDiameter(std::vector<double> position){

	std::vector<double> wakeDiameter(2,0.0);
	
	double x = std::abs( parentTurbine->position[0] - position[0] );
	double Ct = parentTurbine->computeCt();

	double beta = ComputeBeta(Ct);
	double epsilon = 0.2 * std::sqrt(beta);

	scaleFactor = 2.0;
	wakeDiameter[0] =  scaleFactor * parentTurbine->diameter * (wakeExpansionCoefficient[0] 
					  * (x/parentTurbine->diameter) + epsilon );
	wakeDiameter[1] = wakeDiameter[0]; 
	
	return wakeDiameter;
}




std::vector<double> BPModel::GetWakeVelocity(std::vector<double> position){

	std::vector<double> wakeVelocity(3,0.0);

	std::vector<double> hubCoordinates = parentTurbine->GetHubCoordinates();

	double x = std::abs( hubCoordinates[0] - position[0] );
	double y = std::abs( hubCoordinates[1] - position[1] );
	double z = std::abs( hubCoordinates[2] - position[2] );

	double r = std::sqrt( y*y + z*z );
	double Ct = parentTurbine->computeCt();

	double beta = ComputeBeta(Ct);
	double epsilon = 0.2 * std::sqrt(beta);
	double sigma = this->GetWakeDiameter(position)[0] / scaleFactor;

//	double C_x = 1.0 - std::sqrt(1.0 - ( Ct / (8.0 * std::pow(sigma / parentTurbine->diameter ,2.0) ) ));
	double C_x = 1.0 - std::sqrt(1.0 - ( ( Ct * std::pow(parentTurbine->diameter ,2.0 )) / (8.0 * sigma * sigma ) ) ) ;
	double velocityDeficit = C_x * std::exp( - (r*r) / (2.0 * sigma*sigma) );
	
	for(int i = 0; i < 3; i++){	
		wakeVelocity[i] = (1.0 - velocityDeficit) * parentTurbine->velocity[i];		
	}

	return wakeVelocity;
}



/*
		XIE AND ARCHER WAKE MODEL
*/


std::vector<double> XAModel::GetWakeDiameter(std::vector<double> position){

	std::vector<double> wakeDiameter(2,0.0);
	
	double x = std::abs( parentTurbine->position[0] - position[0] );
	double Ct = parentTurbine->computeCt();

	double beta = ComputeBeta(Ct);
	double epsilon = 0.2 * std::sqrt(beta);

	scaleFactor = 2.0;
	wakeDiameter[0] =  scaleFactor*parentTurbine->diameter * (wakeExpansionCoefficient[0] 
					  * (x/parentTurbine->diameter) + epsilon );
					  
	wakeDiameter[1] = scaleFactor*parentTurbine->diameter * (wakeExpansionCoefficient[1] 
					  * (x/parentTurbine->diameter) + epsilon );
		
	return wakeDiameter;
}


std::vector<double> XAModel::GetWakeVelocity(std::vector<double> position){

	std::vector<double> wakeVelocity(3,0.0);

	std::vector<double> hubCoordinates = parentTurbine->GetHubCoordinates();
	
	double x = std::abs( hubCoordinates[0] - position[0] );
	double y = std::abs( hubCoordinates[1] - position[1] );
	double z = std::abs( hubCoordinates[2] - position[2] );

	double Ct = parentTurbine->computeCt();

	double beta = ComputeBeta(Ct);
	double epsilon = 0.2 * std::sqrt(beta);
	
	std::vector<double> sigma = this->GetWakeDiameter(position);
	sigma[0] /= scaleFactor;	sigma[1] /= scaleFactor;


	double deltaHub = 1.0 - std::sqrt(1.0 - ( Ct * std::pow( parentTurbine->diameter,2.0 ) ) / ( 8.0 * sigma[0] * sigma[1]) );
	double velocityDeficit = deltaHub * std::exp( -0.5 * (std::pow( y / sigma[0],2.0) + std::pow( z / sigma[1],2.0)) );
	
	for(int i = 0; i < 3; i++){	
		wakeVelocity[i] = (1.0 - velocityDeficit) * parentTurbine->velocity[i];		
	}
	
	return wakeVelocity;
}



/*
		USER DEFINED WAKE MODEL
*/


std::vector<double> UserDefModel::GetWakeDiameter(std::vector<double> position){


	std::vector<double> wakeDiameter(2,0.0);
	
	UserDefWakeModel userDef(wakeExpansionCoefficient,turbulenceIntensity);
	userDef.T_Ct = parentTurbine->computeCt();
	userDef.T_hubCoordinates = parentTurbine->GetHubCoordinates();
	userDef.T_position = parentTurbine->position;
	userDef.T_velocity = parentTurbine->velocity;				

	userDef.T_height = parentTurbine->height;
	userDef.T_radius = parentTurbine->radius;
	userDef.T_diameter = parentTurbine->diameter;
	userDef.T_area = parentTurbine->area;
	userDef.T_ratedPower = parentTurbine->ratedPower;

	wakeDiameter = userDef.GetWakeDiameter(position);
	
	return wakeDiameter;
}


std::vector<double> UserDefModel::GetWakeVelocity(std::vector<double> position){


	std::vector<double> wakeVelocity(3,0.0);

	UserDefWakeModel userDef(wakeExpansionCoefficient,turbulenceIntensity);
	userDef.T_Ct = parentTurbine->computeCt();
	userDef.T_hubCoordinates = parentTurbine->GetHubCoordinates();
	userDef.T_position = parentTurbine->position;
	userDef.T_velocity = parentTurbine->velocity;				

	userDef.T_height = parentTurbine->height;
	userDef.T_radius = parentTurbine->radius;
	userDef.T_diameter = parentTurbine->diameter;
	userDef.T_area = parentTurbine->area;
	userDef.T_ratedPower = parentTurbine->ratedPower;
	
	
	wakeVelocity = userDef.GetWakeVelocity(position);
	
	
	return wakeVelocity;
}





/*


	WAKE MERGING FUNCTIONS

*/

WakeMergeModel::WakeMergeModel(std::string inMethod){
	method = inMethod;
	
	if (iequals(inMethod, "Linear")) {
		MergeScheme = &WakeMergeModel::LinearMerge;
	}
	else if (iequals(inMethod, "Quadratic")){
		MergeScheme = &WakeMergeModel::QuadraticMerge;
	}
	else if (iequals(inMethod, "Energy")){
		MergeScheme = &WakeMergeModel::EnergyMerge;
	}
	else if (iequals(inMethod, "DWM")){
		MergeScheme = &WakeMergeModel::DWMMerge;
	}
	else if (iequals(inMethod, "UserDef")){
		MergeScheme = &WakeMergeModel::UserDefMerge;				
	}
	else{
		std::cout << "Incorrect wake model" << std::endl;
		exit(10);
	}
}

WakeMergeModel::WakeMergeModel(TurbineProps *inParentTurbine, std::string inMethod){
	
	method = inMethod;
	this->SetParent(inParentTurbine);
	
	if (iequals(inMethod, "Linear")) {
		MergeScheme = &WakeMergeModel::LinearMerge;
	}
	else if (iequals(inMethod, "Quadratic")){
		MergeScheme = &WakeMergeModel::QuadraticMerge;
	}
	else if (iequals(inMethod, "Energy")){
		MergeScheme = &WakeMergeModel::EnergyMerge;
	}
	else if (iequals(inMethod, "DWM")){
		MergeScheme = &WakeMergeModel::DWMMerge;
	}
	else if (iequals(inMethod, "UserDef")){
		MergeScheme = &WakeMergeModel::UserDefMerge;				
	}
	else{
		std::cout << "Incorrect wake model" << std::endl;
		exit(10);
	}
}



void WakeMergeModel::SetParent(TurbineProps *inTurbine){
	parentTurbine = inTurbine;
}

std::vector<double> WakeMergeModel::Merge(){


	std::vector<double> hubCenter = parentTurbine->GetHubCoordinates();
	std::vector<double> ambientWind = parentTurbine->ComputeAverageAmbientVelocity();
	
	std::vector<double> modHubCenter(3,0.0);
	std::vector<double> sumOfDeficits(3,0.0);
	std::vector<double> wakeVelocity(3,0.0);

	double overlapRatio = 0.0;

	if(parentTurbine->influencedBy.size() == 0){
		return ambientWind;
	}

	for(int i = 0; i < parentTurbine->influencedBy.size(); i++){
		
		if(parentTurbine->influencedBy[i]->fictitious){
			continue;
		}
		
		modHubCenter = GetArcLengthVector(parentTurbine->influencedBy[i]->position,hubCenter);
	
		overlapRatio = parentTurbine->influencedBy[i]->wakeModel->ComputeWakeOverLapArea(modHubCenter, parentTurbine->radius) / parentTurbine->area;	
		wakeVelocity = parentTurbine->influencedBy[i]->wakeModel->AverageWakeVelocity(modHubCenter, parentTurbine->radius);
		
		CheckNanVelocity(wakeVelocity);
						
		(this->*MergeScheme)(sumOfDeficits, ambientWind, wakeVelocity, overlapRatio,false);
	}

	(this->*MergeScheme)(sumOfDeficits, ambientWind, wakeVelocity,overlapRatio, true);
	
	return sumOfDeficits;
}



void WakeMergeModel::LinearMerge(std::vector<double> &velocitySum,
												std::vector<double> ambientWind, 
												std::vector<double> wakeVelocity, 
												double areaRatio,bool reverse){
	if(reverse){
		for(int i = 0; i < 3; i++){
			velocitySum[i] = ambientWind[i] - velocitySum[i];
		}
		return;
	}								
				
	for(int i = 0; i < 3; i++){
		velocitySum[i] += areaRatio*(ambientWind[i] - wakeVelocity[i]);
	}
}

void WakeMergeModel::QuadraticMerge(std::vector<double> &velocitySum,
												std::vector<double> ambientWind, 
												std::vector<double> wakeVelocity,
												double areaRatio,bool reverse){
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

void WakeMergeModel::EnergyMerge(std::vector<double> &velocitySum,
												std::vector<double> ambientWind, 
												std::vector<double> wakeVelocity,
												double areaRatio,bool reverse){
	if(reverse){
		for(int i = 0; i < 3; i++){
			velocitySum[i] = std::sqrt( std::abs( std::pow(ambientWind[i],2.0) - velocitySum[i] ));			
		}
		return;
	}								
				
	for(int i = 0; i < 3; i++){
		velocitySum[i] += areaRatio*(std::pow(ambientWind[i],2.0) - std::pow(wakeVelocity[i],2.0));
	}
	

	
}

void WakeMergeModel::DWMMerge(std::vector<double> &velocitySum,
											 std::vector<double> ambientWind, 
											 std::vector<double> wakeVelocity,
											 double areaRatio,bool reverse){
	if(reverse){
		for(int i = 0; i < 3; i++){
			velocitySum[i] = ambientWind[i] - velocitySum[i];
		}
		return;
	}								
				
	for(int i = 0; i < 3; i++){
		velocitySum[i] = areaRatio*(std::max(velocitySum[i], ambientWind[i] - wakeVelocity[i]));
	}
}




/*
		USER DEFINED WAKE MERGE MODEL
*/


void WakeMergeModel::UserDefMerge(std::vector<double> &velocitySum,
											 std::vector<double> ambientWind, 
											 std::vector<double> wakeVelocity,
											 double areaRatio,bool reverse){

	UserDefWakeMergeModel(velocitySum,ambientWind, wakeVelocity, areaRatio, reverse);
	return;
}



