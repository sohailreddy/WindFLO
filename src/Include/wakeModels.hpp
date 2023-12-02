/*

	Copyright 2019, Sohail R. Reddy
	email: sredd001@fiu.edu
	www.sohailreddy.com

*/


#ifndef WAKEMODELS_HEADER
#define WAKEMODELS_HEADER

#include <vector>
#include <string>
#include <algorithm>
#include <iostream>



class TurbineProps;


							   
class WakeModel
{
	public :
		WakeModel();
//		WakeModel(const WakeModel *inWakeModel);		
		
//	Variables
	TurbineProps *parentTurbine;
	double turbulenceIntensity;
	double scaleFactor;
	std::vector<double> wakeExpansionCoefficient;
	std::string modelName;
	int monteCarloPts;
	
//	Functions
	WakeModel *getWakeModelPointer(std::string modelName);	
	virtual void SetParameter(std::vector<double> in_WakeExpansionCoefficient, double in_TurbulenceIntensity);
	virtual void SetParent(TurbineProps *inTurbine);
	virtual std::vector<double> GetWakeDiameter(std::vector<double> position);
	virtual std::vector<double> GetWakeVelocity(std::vector<double> position);
	virtual bool CheckIfInsideWakeRegion(std::vector<double> position);
	virtual double ComputeWakeOverLapArea(std::vector<double> position, double radius);
	virtual std::vector<double> AverageWakeVelocity(std::vector<double> center, double radius);
};



class JensenModel : public WakeModel
{
	public:
	JensenModel() : WakeModel() {modelName = "Jensen";};
	
	virtual std::vector<double> GetWakeDiameter(std::vector<double> position) override;
	virtual std::vector<double> GetWakeVelocity(std::vector<double> position) override;
};


class FrandsenModel : public WakeModel
{
	public:
	explicit FrandsenModel() : WakeModel() {modelName = "Frandsen";};
	
	virtual std::vector<double> GetWakeDiameter(std::vector<double> position) override;
	virtual std::vector<double> GetWakeVelocity(std::vector<double> position) override;
};

class LarsenModel : public WakeModel
{
	public:
	LarsenModel() : WakeModel() {modelName = "Larsen";};
	
	virtual std::vector<double> GetWakeDiameter(std::vector<double> position) override;
	virtual std::vector<double> GetWakeVelocity(std::vector<double> position) override;
};


class IshiharaModel : public WakeModel
{
	public:
	IshiharaModel() : WakeModel() {modelName = "Ishihara";};
	
	virtual std::vector<double> GetWakeDiameter(std::vector<double> position) override;
	virtual std::vector<double> GetWakeVelocity(std::vector<double> position) override;
};


class BPModel : public WakeModel
{
	public:
	BPModel() : WakeModel() {modelName = "BP";};
	
	virtual std::vector<double> GetWakeDiameter(std::vector<double> position) override;
	virtual std::vector<double> GetWakeVelocity(std::vector<double> position) override;
};


class XAModel : public WakeModel
{
	public:
	XAModel() : WakeModel() {modelName = "XA";};
	
	virtual std::vector<double> GetWakeDiameter(std::vector<double> position) override;
	virtual std::vector<double> GetWakeVelocity(std::vector<double> position) override;
};

class UserDefModel : public WakeModel
{
	public:
	UserDefModel() : WakeModel() {modelName = "UserDef";};
	
	virtual std::vector<double> GetWakeDiameter(std::vector<double> position) override;
	virtual std::vector<double> GetWakeVelocity(std::vector<double> position) override;
};


class WakeMergeModel
{
	public:
	WakeMergeModel(){};	
	WakeMergeModel(TurbineProps *inParentTurbine, std::string inMethod);
	WakeMergeModel(std::string inMethod);	

//	Variables
	TurbineProps *parentTurbine;
	std::string method;
		
//	Functions	
	virtual void SetParent(TurbineProps *inTurbine);
	virtual std::vector<double> Merge();
	
	void (WakeMergeModel::*MergeScheme)(std::vector<double> &,
						std::vector<double> , 
						std::vector<double> , double, bool );
											
	virtual void LinearMerge(std::vector<double> &velocitySum,
							 std::vector<double> ambientWind, 
							 std::vector<double> wakeVelocity,
							 double areaRatio = 1.0,bool reverse = false);
											
	virtual void QuadraticMerge(std::vector<double> &velocitySum,
								std::vector<double> ambientWind, 
								std::vector<double> wakeVelocity,
								double areaRatio = 1.0,bool reverse = false);
											   
	virtual void EnergyMerge(std::vector<double> &velocitySum,
							 std::vector<double> ambientWind, 
							 std::vector<double> wakeVelocity,
							 double areaRatio = 1.0,bool reverse = false);
											
	virtual void DWMMerge(std::vector<double> &velocitySum,
						  std::vector<double> ambientWind, 
						  std::vector<double> wakeVelocity,
						  double areaRatio = 1.0,bool reverse = false);
						  
						  
	virtual void UserDefMerge(std::vector<double> &velocitySum,
						  std::vector<double> ambientWind, 
						  std::vector<double> wakeVelocity,
						  double areaRatio = 1.0,bool reverse = false);
						  
};




WakeModel *getWakeModelPointer(std::string modelName);
double ComputeBeta(double Ct);
double ComputeOverLappingArea(std::vector<double> centerOne, std::vector<double> centerTwo,
							   double radiusOne, double radiusTwo);





#endif