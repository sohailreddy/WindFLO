/*

	Copyright 2019, Sohail R. Reddy
	email: sredd001@fiu.edu
	www.sohailreddy.com

*/


#include "WindFLO.hpp"
#include "input.hpp"
#include "convexHull.hpp"
#include "IO.hpp"
#include "math.hpp"

#include <iostream>
#include <fstream>


//std::vector<TurbineProps> turbines;

/*
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
*/

double norm2(std::vector<double> a);
double dot_product(std::vector<double> a, std::vector<double> b);
std::vector<double> matmul(std::vector<std::vector<double> > a, std::vector<double> x);



extern double undefined = 9876543210.0;
extern int undefined_i = 987654321;
extern double PI = 3.14159265359;



WindFLO::WindFLO(std::string inFileName){

	if(inFileName != ""){
		fileName = inFileName;
	}
	else {
		fileName = "WindFLO.dat";
	}
	GetInput();		
}

WindFLO::~WindFLO(){
	//CleanAll();
}


void WindFLO::setup(){
	turbines.clear();
	turbines.resize(nTurbines);
	for (int i = 0; i < nTurbines; i++)
		GetTurbine(i);
}


void WindFLO::Clean(){
	Clean_();
}
void WindFLO::CleanAll(){
	CleanAll_();
}


void WindFLO::GetInput(){

	char cFileName[100] = "";
	strcpy(cFileName, fileName.c_str());
	ReadInput_(cFileName);
	
	GetRunParameters();
	GetAtmParameters();
	GetAmbientWindParameters();
	GetWakeParameters();
	GetNumericalParameters();
	GetCostParameters();	
	
}



void WindFLO::GetRunParameters(){
	double runParams[50];
	GetRunParameters_(runParams);
	
	nTurbines = (int)runParams[0];
	batch = (int)runParams[1];
	nWindRose = (int)runParams[2];	
}

void WindFLO::GetAtmParameters(){
	double atmParams[50];
	GetAtmParameters_(atmParams);
	rho = atmParams[0];
	turbulenceIntensity = atmParams[1];
}

void WindFLO::GetAmbientWindParameters(){
	char modelName[100] = "";
	double ambientParams[50];

	GetAmbientWindParameters_(modelName, ambientParams);
	windModel = modelName;
	
	modelVelocity.resize(3);
	modelVelocity[0] = ambientParams[0];
	modelVelocity[1] = ambientParams[1];
	modelVelocity[2] = ambientParams[2];
	
	surfaceRoughness = ambientParams[3];
	referenceHeight = ambientParams[4];
			
	return;
}

void WindFLO::GetWakeParameters(){

	char wakeModelName[100] = "";	
	char wakeMergeModelName[100] = "";	
	double wakeParams[50];

	GetWakeParameters_(wakeModelName, wakeMergeModelName , wakeParams);
	
	wakeModel = wakeModelName;
	wakeMergeModel = wakeMergeModelName;
	
	wakeExpansionCoeff.resize(2);
	wakeExpansionCoeff[0] = wakeParams[0];
	wakeExpansionCoeff[1] = wakeParams[1];	
}

void WindFLO::GetNumericalParameters(){

	int numericalParams[50];
	GetNumericalParameters_(numericalParams);
	
	gaussOrder = numericalParams[0];
	monteCarloPts = numericalParams[1];
}

void WindFLO::GetCostParameters(){

	double costParams[50];
	GetCostParameters_(costParams);
	coe = costParams[0];
}
	
	
	
	
void WindFLO::GetTurbine(int ith){

	int ithTurbine = ith + 1;	

	char name[100] = ""; 
	double turbineParams[50];
	double Ctx[5000], Cty[5000];
	double Cpx[5000], Cpy[5000];
	
	
	GetTurbine_(&ithTurbine, name, turbineParams, Ctx, Cty, Cpx, Cpy);
	
	turbines[ith] = TurbineProps(name, turbineParams);
	std::vector<double> x, y;
	for(int i = 0; i < 5000; i++ ){
		if(Ctx[i] == undefined || Cty[i] == undefined) break;
		x.push_back(Ctx[i]);
		y.push_back(Cty[i]);
	}
	turbines[ith].CtSpline.set_points(x,y);	

	x.clear();
	y.clear();

	for(int i = 0; i < 5000; i++ ){
		if(Cpx[i] == undefined || Cpy[i] == undefined) break;
		x.push_back(Cpx[i]);
		y.push_back(Cpy[i]);
	}

	turbines[ith].CpSpline.set_points(x,y);	

	GetElevation_(&turbines[ith].position[0], 
				  &turbines[ith].position[1], 
				  &turbines[ith].position[2]);

	turbines[ith].SetWakeModel(wakeModel,wakeExpansionCoeff, turbulenceIntensity );
	turbines[ith].wakeModel->monteCarloPts = monteCarloPts;
	
	
	turbines[ith].SetWakeMergeModel(wakeMergeModel);
	
	turbines[ith].SetAmbientWindModel(windModel, &modelVelocity, 
				  				   surfaceRoughness, referenceHeight);	
				  				   		
	turbines[ith].ambientWindModel.gaussOrder = gaussOrder;
	
	if(turbines[ith].ratedPower <= 0.0) {
		turbines[ith].ratedPower = turbines[ith].computeRatedPower(rho);
	}
	
	
	// used to unsort turbines after rank sorting
	turbines[ith].dummy = ith;

}


void WindFLO::run(){


	double probability = 1.0;
	int i = 0;	
	AEP = 0.0;
	double dVdTheta = 1.0;
	InitializeOutputs();
	do {
		GetWindRoseVelocity(i, modelVelocity, probability,dVdTheta);

		i += 1;
		
		ComputeTurbineElevation();
		PerformCoordinateTransformation(true);
		CheckOrientationOfTurbines();
		ComputeRank();
		ComputeVelocities();
		PerformCoordinateTransformation(false);
		
		if(i == nWindRose){
			ComputeTurbinePower(1.0/nWindRose);
		} else{
			ComputeTurbinePower();
		}
		
		ComputeAEP(probability);
		
	} while (i < nWindRose);

	ComputeFarmPower();
	ComputeFarmEfficiency();
	ComputeLandUsage();
	ComputeFarmCost();
	UnSortTurbines();
	
}


void WindFLO::InitializeOutputs(){
	for(int i = 0; i < turbines.size() ; i++){
		turbines[i].power = 0.0;
	}
	landUsed = 0.0;
	farmEfficiency = 0.0;
	farmPower = 0.0 ;
	farmCost = 0.0;
	AEP = 0.0;
}


void WindFLO::GetWindRoseVelocity(int &ith, std::vector<double> &velocity, double &probability, double &dVdTheta){

	double windRoseParams[50];
	
	if(nWindRose == 0) {
		velocity = modelVelocity;
		probability = 1.0;
		dVdTheta = 1.0;	
		return;
	}
	
	int ithRose = ith + 1;
	GetWindRoseParameters_(&ithRose,windRoseParams);
	velocity[0] = windRoseParams[0];
	velocity[1] = windRoseParams[1];
	velocity[2] = windRoseParams[2];
	probability = windRoseParams[3];
	dVdTheta = windRoseParams[4];
		
	return;
}

void WindFLO::CheckOrientationOfTurbines(){

	double vnorm = norm2(this->modelVelocity);

	for(int i = 0; i < nTurbines ; i++){

//	if flow is coming from behind the turbine and it cannot yaw
//	it is set as ficticous	
		if(!turbines[i].yaw){
			double ans = dot_product( this->modelVelocity, turbines[i].orientation );
			if(ans < 0.0){
				turbines[i].fictitious = true;
			}
		}

//	if it can yaw then set its orientation (normal) to be the 
//	opposite direction of the wind direction
		else {
			for(int j = 0; j < this->modelVelocity.size() ; j++){
				turbines[i].orientation[j] = -this->modelVelocity[j] /vnorm;
			}
		}		
	}
	
	return;
}



void WindFLO::ComputeRotationMatrix(){

//	Compute rotation matrix	and rotate velocity
	int k = 0;
	double rotMat[9];
	double a[3] = {modelVelocity[0],modelVelocity[1],modelVelocity[2]};	
	double x[3] = {1.0,0.0,0.0};
		
	GetRotationMatrix_(a,x,rotMat);
	rotationMat.resize(3);
	for(int i = 0; i < 3 ; i++){
		rotationMat[i].resize(3);
		for(int j = 0; j < 3 ; j++, k++){
			if(abs(rotMat[k]) < 1.0e-10) rotMat[k] = 0.0;
			rotationMat[i][j] = rotMat[k];
		}
	}
	
//	Compute rotation matrix for the inverse transformation
	MatInv3_(rotMat);
	rotationMat_inv.resize(3);	
	k = 0;
	for(int i = 0; i < 3 ; i++){
		rotationMat_inv[i].resize(3);
		for(int j = 0; j < 3 ; j++, k++){
			if(abs(rotMat[k]) < 1.0e-10) rotMat[k] = 0.0;
			rotationMat_inv[i][j] = rotMat[k];
		}
	}
	
	return;
}

void WindFLO::PerformCoordinateTransformation(bool flag){
	
	if(flag){	
		ComputeRotationMatrix();
		for(int i = 0; i < nTurbines ; i++){
			turbines[i].position = matmul(rotationMat, turbines[i].position);
			turbines[i].orientation = matmul(rotationMat, turbines[i].orientation);
			turbines[i].velocity = matmul(rotationMat, turbines[i].velocity);			
		}			
		this->modelVelocity = matmul(rotationMat, this->modelVelocity);		

	}
	else{
		for(int i = 0; i < nTurbines ; i++){
			turbines[i].position = matmul(rotationMat_inv, turbines[i].position);
			turbines[i].orientation = matmul(rotationMat_inv, turbines[i].orientation);
			turbines[i].velocity = matmul(rotationMat_inv, turbines[i].velocity);			
		}
		this->modelVelocity = matmul(rotationMat_inv, this->modelVelocity);
	}
	
	return;
}


void WindFLO::ComputeTurbineElevation(){

	for(int i = 0; i < nTurbines ; i++){
		GetElevation_(&turbines[i].position[0], 
					  &turbines[i].position[1], 
					  &turbines[i].position[2]);
	}
	return;
}

void WindFLO::ComputeRank(){

	auto comp = [](TurbineProps &a, TurbineProps &b) {
		return (a.position[0] < b.position[0]);
	};
	std::sort(turbines.begin(), turbines.end(), comp);

	int rank = 0;
	turbines[0].rank = rank;
	turbines[0].RunCheck();		
	for(int i = 1; i < nTurbines ; i++){
		if(turbines[i].position[0] > turbines[i-1].position[0]){
			rank += 1;
		}		
	
		turbines[i].rank = rank;
		turbines[i].RunCheck();			
	}
	
	return;
}

void WindFLO::UnSortTurbines(){

	auto comp = [](TurbineProps &a, TurbineProps &b) {
		return (a.dummy < b.dummy);
	};
	std::sort(turbines.begin(), turbines.end(), comp);

	for(int i = 0; i < nTurbines ; i++){
		turbines[i].dummy = 0.0;
		turbines[i].RunCheck();			
	}
	return;
}


void WindFLO::ComputeVelocities(){

	for(int i = 0; i < nTurbines ; i++){
	
		if(turbines[i].rank == 0){
			turbines[i].velocity = turbines[i].ComputeAverageAmbientVelocity();
			
		}
		else{
		
			std::vector<double> hubLocation = turbines[i].GetHubCoordinates();
			double radius = turbines[i].radius;

			for(int j = 0; j < i ; j++){
				if( turbines[i].rank > turbines[j].rank && (!turbines[j].fictitious) ){
					
					if( turbines[j].wakeModel->ComputeWakeOverLapArea(hubLocation,radius) > 0.0 ) {
						turbines[i].influencedBy.push_back( &turbines[j] );
					}
				}
			}

			turbines[i].velocity = turbines[i].wakeMergeModel.Merge();
			turbines[i].influencedBy.clear();
		}	
				
	}
	
	return;
}


void WindFLO::ComputeLandUsage(){

	std::vector<std::vector<double> > xy;
	
	for(int i = 0; i < nTurbines ; i++){
		if(!turbines[i].fictitious){
			xy.push_back(turbines[i].position);	
		}
	}

	convexHull = computeConvexHull(xy);
	convexHull.push_back(convexHull[0]);
	landUsed = ComputeConvexHullArea(convexHull);
}


void WindFLO::ComputeFarmCost(){

	CostModel costModel(coe);
	farmCost = 0.0;
	for(int i = 0; i < nTurbines ; i++){
		if(!turbines[i].fictitious) {	
			farmCost += costModel.ComputeCost(turbines[i].diameter,turbines[i].height);
		}
	}	
}



void WindFLO::ComputeTurbinePower(double scale){

	for(int i = 0; i < nTurbines ; i++){
		if(!turbines[i].fictitious) {
			turbines[i].power += turbines[i].computePower(rho);
			turbines[i].power *= scale;
		}
	}
}


void WindFLO::ComputeFarmPower(){

	farmPower = 0.0;
	for(int i = 0; i < nTurbines ; i++){
		if(!turbines[i].fictitious) {
			farmPower += turbines[i].power;
		}
	}	
}

void WindFLO::ComputeAEP(double probability, double dVdTheta ){


	double tmpFarmPower = 0.0;
	for(int i = 0; i < nTurbines ; i++){
		if(!turbines[i].fictitious) {
			tmpFarmPower += turbines[i].computePower(rho);
		}
	}	

	AEP += (tmpFarmPower*probability)*(365.0*24.0) * dVdTheta;
}


void WindFLO::ComputeFarmEfficiency(){

	double totalRatedPower = 0.0;
	farmEfficiency = 0.0;
	double totalFarmPower = 0.0;

	for(int i = 0; i < nTurbines ; i++){
		if(!turbines[i].fictitious) {
			totalFarmPower += turbines[i].power;
			totalRatedPower += turbines[i].ratedPower;
		}
	}	
	farmEfficiency = totalFarmPower / totalRatedPower;
}



void WindFLO::WriteResults(std::string filename){

	if(filename == "") filename = "WindFLO.res";

	std::ofstream resFile;
	resFile.open (filename);
	resFile << "$Farm \n";
	resFile <<"nTurbines = " << nTurbines << "\n";	
	resFile <<"AEP = " << AEP << "\n";	
	resFile <<"Power = " << farmPower << "\n";	
	resFile <<"Efficiency = " << farmEfficiency << "\n";	
	resFile <<"Cost = " << farmCost << "\n";	
	resFile <<"Land Used = " << landUsed << "\n";	
	

	resFile << "$Turbines \n";	
	resFile << "i, x, y, z, theta, phi, chi, ux, uy, uz, H, R, A, P_r, P, fict \n";
	for(int i = 0; i < nTurbines; i++){
	
		resFile << i+1 << ", ";	
		resFile << turbines[i].position[0] << ", ";
		resFile << turbines[i].position[1] << ", ";
		resFile << turbines[i].position[2] << ", ";

		resFile << turbines[i].orientation[0] << ", ";
		resFile << turbines[i].orientation[1] << ", ";
		resFile << turbines[i].orientation[2] << ", ";


		resFile << turbines[i].velocity[0] << ", ";
		resFile << turbines[i].velocity[1] << ", ";
		resFile << turbines[i].velocity[2] << ", ";

		resFile << turbines[i].height << ", ";
		resFile << turbines[i].radius << ", ";
		resFile << turbines[i].area << ", ";
		resFile << turbines[i].ratedPower << ", ";
		
		resFile << turbines[i].power << ", ";
		resFile << turbines[i].fictitious << " \n";
		
	}
	
	resFile << "$ConvexHull \n";	
	resFile << "i, x, y, z \n";
	for(int i = 0; i < convexHull.size() ; i++){
		resFile << i+1 << ", ";		
		resFile << 	convexHull[i][0] << ", ";
		resFile << 	convexHull[i][1] << " \n";
//		resFile << 	convexHull[i][2] << " \n";
	}
	resFile << "$End \n";	
	
	resFile.close();

}


/*

	Math routines

*/


std::vector<double> GetArcLengthVector(std::vector<double> A, std::vector<double> B){

	double a[3] = {A[0],A[1],A[2]};
	double b[3] = {B[0],B[1],B[2]};
	double length = 0.0;
	
	ComputeArcLength_(a,b,&length);
	
	std::vector<double> c(3,0.0);
	
	c[0] = A[0] + length;
	c[1] = B[1];
	c[2] = B[2];
		
	return c;
}

double norm2(std::vector<double> a){
	double ans = 0.0;
	for(int i = 0; i < a.size(); i++){
		ans += (a[i] * a[i]);
	}
	return sqrt(ans);
}

std::vector<double> matmul(std::vector<std::vector<double> > a, std::vector<double> x){

	int n = x.size();
	std::vector<double> b(n,0.0);
		
	for (int i=0; i<n; i++){
		b[i] = 0.0;
		for (int j=0; j<n; j++){
			b[i]+=( a[j][i]*x[j]);
		}
		if(abs(b[i]) < 1.0e-15) b[i] = 0.0;
	}
	return b;
}

double dot_product(std::vector<double> a, std::vector<double> b){
	int n = a.size();
	double ans = 0.0;
	for(int i = 0; i < n; i++){
		ans += (a[i] * b[i]);
	}
	return ans;
}










