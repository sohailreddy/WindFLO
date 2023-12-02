/*

	Copyright 2019, Sohail R. Reddy
	email: sredd001@fiu.edu
	www.sohailreddy.com

*/



#include "ambientWindModels.hpp"
#include "input.hpp"

extern double PI;

bool iequals(const std::string& a, const std::string& b);

void getGaussPoints(int order, std::vector<double> &rTheta, std::vector<double> &weights);

AmbientWindModel::AmbientWindModel(std::string inModel, std::vector<double> *inVelocity, 
								   double inSurfaceRoughness, double inReferenceHeight){

	modelVelocity = inVelocity;
	
	if (iequals(inModel, "log")) {
		model = "log";
		surfaceRoughness = inSurfaceRoughness;	
	}
	else if (iequals(inModel, "power")){
		model = "power";	
		referenceHeight = inReferenceHeight;
	}
	else if (iequals(inModel, "constant")){
		model = "constant";	
	}	
	else if (iequals(inModel, "deaves-harris")){
		model = "deaves-harris";	
		surfaceRoughness = inSurfaceRoughness;	
		referenceHeight = inReferenceHeight;
	}		
	else{
		std::cout << "Incorrect Ambient Wind Model \n" ;
		exit(10);
	}
}



std::vector<double> AmbientWindModel::ComputeVelocity(double height){

	if (iequals(model, "log")) {
		return this->LogLawProfile(height);
	}
	else if (iequals(model, "power")){
		return this->PowerLawProfile(height);
	}
	else if (iequals(model, "constant")){
		return this->ConstantProfile();
	}
	else if (iequals(model, "deaves-harris")){
		return this->DeavesHarrisProfile(height);	
	}		

	

}

std::vector<double> AmbientWindModel::ConstantProfile(double height){
	return *modelVelocity;
}


std::vector<double> AmbientWindModel::LogLawProfile(double height){

	std::vector<double> velocity(3,0.0);

	for(int i = 0; i < 3; i++){
		velocity[i] = ((*modelVelocity)[i] / vonKarmanConstant) * std::log(height / surfaceRoughness);
	}
	
	return velocity;
}


std::vector<double> AmbientWindModel::PowerLawProfile(double height){

	std::vector<double> velocity(3,0.0);

	for(int i = 0; i < 3; i++){
		velocity[i] = (*modelVelocity)[i] * std::pow( height / referenceHeight, 0.143 ) ;
	}
	
	return velocity;
}

std::vector<double> AmbientWindModel::DeavesHarrisProfile(double height){

	std::vector<double> velocity(3,0.0);

	for(int i = 0; i < 3; i++){
		velocity[i] = ( (*modelVelocity)[i] / vonKarmanConstant )
						* ( std::log( height / surfaceRoughness )
						+ 5.75 * ( height / referenceHeight )
						- 1.88 * std::pow( height / referenceHeight, 2.0 )
						- 1.33 * std::pow( height / referenceHeight, 3.0 )
						+ 0.25 * std::pow( height / referenceHeight, 4.0 )	) ;
	}
	
	return velocity;
}




std::vector<double> AmbientWindModel::AverageVelocity(std::vector<double> center, double radius){

	std::vector<double> averageVelocity(3,0.0);

	std::vector<double> rTheta;
	std::vector<double> weights;

	if(radius <= 0.0){
		averageVelocity = this->ComputeVelocity(center[2]);
		return averageVelocity;
	}

	getGaussPoints(gaussOrder, rTheta, weights);
	int nPoints = rTheta.size();

	double r, theta, height;
	std::vector<double> windProfile(3,0.0);

	for(int i = 0; i < nPoints ; i++){
		for(int j = 0; j < nPoints ; j++){
	
			r = 0.5 * (rTheta[i] + 1.0) * radius;
			theta = PI * (rTheta[j] + 1.0);
	
			height = r * std::sin(theta) + center[2];
			windProfile = this->ComputeVelocity(height);
					
			for(int k = 0; k < 3; k++){
				averageVelocity[k] += weights[i]*weights[j]*windProfile[k]*(rTheta[i] + 1.0);
			}
	
		}
	}

	for(int k = 0; k < 3; k++){
		averageVelocity[k] *= 0.25;
	}


	return averageVelocity;
}


void getGaussPoints(int order, std::vector<double> &rTheta, std::vector<double> &weights){

	if(order == 1) {
		rTheta.resize(1); weights.resize(1);
		rTheta[0] = 0.0; weights[0] = 2.0;
	}
	else if (order == 2){
		rTheta.resize(2); 		weights.resize(2);
		rTheta[0] = -0.57735; 	weights[0] = 1.0;
		rTheta[1] = 0.57735	; 	weights[1] = 1.0;			
	}
	else if (order == 3){
		rTheta.resize(3); weights.resize(3);
		rTheta[0] = -0.774597; 	weights[0] = 0.555556;	
		rTheta[1] = 0.0		; 	weights[1] = 0.888889;
		rTheta[2] = 0.774597 ; 	weights[2] = 0.555556;
	}
	else if (order == 4){
		rTheta.resize(4); weights.resize(4);
		rTheta[0] = -0.339981043584856	; weights[0] = 0.652145154862546;	
		rTheta[1] = -0.861136311594053	; weights[1] = 0.347854845137454;
		rTheta[2] = 0.339981043584856	; weights[2] = 0.652145154862546;
		rTheta[3] = 0.861136311594053	; weights[3] = 0.347854845137454;
	}
	else if (order == 5){
		rTheta.resize(5); weights.resize(5);
		rTheta[0] = 0.0;		 	weights[0] = 0.568889;	
		rTheta[1] = 0.538469; 		weights[1] = 0.478629;
		rTheta[2] = -0.538469; 		weights[2] = 0.478629;
		rTheta[3] = 0.90618; 		weights[3] = 0.236927;
		rTheta[4] = -0.90618; 		weights[4] = 0.236927;
	}
	else{
		std::cout << "Incorrect order of Gauss Integration" << std::endl;
		exit(10);
	}
	
	return;
}

		
