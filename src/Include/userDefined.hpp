/*

	Copyright 2019, Sohail R. Reddy
	email: sredd001@fiu.edu
	www.sohailreddy.com

*/

#include <vector>
#include <cmath>
#include <string>
#include <iostream>



class UserDefWakeModel {


	public:
		UserDefWakeModel(std::vector<double> inWakeExpCoeff, double inTurbIntensity){
			wakeExpansionCoeff = inWakeExpCoeff;
			turbIntensity = inTurbIntensity;
		}
		
//		Wake model properties		
		std::vector<double> wakeExpansionCoeff;
		double turbIntensity;
		
//		Turbine properties				
		std::vector<double> T_hubCoordinates; // location of hub
		std::vector<double> T_position;		// position of turbine
		std::vector<double> T_velocity;		// velocity of turbine
		double T_Ct;						// coefficient of thrust at the velocity
		double T_height;					// height of turbine
		double T_radius;					// radius of turbine
		double T_diameter;					// diameter of turbine
		double T_area;						// rotor swept area
		double T_ratedPower;				// rated power of turbine


		
		std::vector<double> GetWakeDiameter(std::vector<double> position);
		std::vector<double> GetWakeVelocity(std::vector<double> position);	

};

void UserDefWakeMergeModel(std::vector<double> &velocitySum,
											 std::vector<double> ambientWind, 
											 std::vector<double> wakeVelocity,
											 double areaRatio,bool reverse);