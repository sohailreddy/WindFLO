/****************************************************************************

	Example 2 on using WindFLO library

  Purpose:
  	
  	This example demonstrate how to use the C++ WindFLO API for optimization of
  	wind farm layout. It uses WindFLO to analyze each layout configuration and
  	a simple genetic algorithm (GA) to perform the optimization. The layout is
  	optimized for maximum power generation and does not incorporate constraints.

  	IMPORTANT: The GA maximizes the fitness function, therefore the positive
  			   farm power is returned by the EvaluateFarm function
  
  Licensing:
  
	This code is distributed under the Apache License 2.0 
	
  Author:
  	Sohail R. Reddy
  	sredd001@fiu.edu
  	
****************************************************************************/


#include <iostream>
#include <vector>


// Include the API.hpp header file
#include "API.hpp"


using namespace std;
// 	Function definitions for the Genetic Algorithm
vector<double> runGA (int maxgens, double pmutation, double pxover );
void initialize ( int npop, vector<double> lbound, vector<double> ubound, double (*fitnessFunc)(vector<double> x), int &seed );


//  WindFLO parameters
#define nTurbines 25	// define the number of turbines
WindFLOAPI windFLO("WindFLO.dat");		// initialize the windFLOAPI class
double EvaluateFarm (vector<double> x);	// Function to evaluate each farm layout


//****************************************************************************80
int main ( ) {

	std::cout << "Running Example 2 \n " ;

	// Two variable per turbines (its x and y coordinates)
	vector<double> lbound(nTurbines*2);		// lower bound of x and y 
	vector<double> ubound(nTurbines*2);		// upper bound of x and y

	// Initialize the upper and lower bounds of each variable
	for(int i = 0; i < nTurbines*2;i++){
		lbound[i] = 0.0;
		ubound[i] = 2000.0;
	}
	
	// Random number generator seed
	int seed = 123456789;  
	
	// Number of population in the GA
	int npop = 50;
	// Max number of generations
	int maxgens = 10;
	// Probability of mutation
	double pmutation =0.15;
	// Probability of crossover
	double pxover = 0.8;
	
	// Initialize the GA... the initialize function takes EvaluateFarm to assign
	// a function pointer to calculate the fitness of each member
	initialize (npop, lbound, ubound, &EvaluateFarm,seed );
	
	// Run the GA ... the optimum solution is returned in xBest
	vector<double> xBest = runGA(maxgens, pmutation, pxover);
	
	// Evaluate the best farm layout
	EvaluateFarm (xBest);
	
	
	std::cout << "Writing optimum \n" ;
	// Write the optimum configuration
	windFLO.write("WindFLO.res");
	
	std::cout << "End execution \n" ;	

	return 0;
}


//****************************************************************************80
// Function to evaluate each farm layout using the WindFLO created above
double EvaluateFarm (vector<double> x)
{
	int nVar = x.size();
	int k = 0;
	for (int i = 0; i < nTurbines; i++ ) {
	   		for (int j = 0; j < 2; j++, k++ ){
	   		// unroll the variable vector 'x' and assign it to turbine positions
			(*windFLO.turbines[i].position)[j] = x[k];	
		}
	}

	// Run WindFLO analysis
	windFLO.run();
			
	// Return the farm power or any other farm output
	return (*windFLO.farmPower);
}
//****************************************************************************80