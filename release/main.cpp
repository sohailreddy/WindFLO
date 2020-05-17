/*

	Copyright 2019, Sohail R. Reddy
	email: sredd001@fiu.edu
	www.sohailreddy.com

*/


#include <vector>
#include <string>
#include <cmath>
#include <iostream>
#include <fstream>


#include "API.hpp"


int main(int argc, char** argv) {
	std::cout << "Copyright 2019, Sohail R. Reddy   (sredd001@fiu.edu)" << std::endl;
	std::string inFile;
	std::string outFile;
	
	if(argc > 1)
		inFile = argv[1];
	if(argc > 2)
		outFile = argv[2];

	WindFLOAPI windFLO(inFile);	
	windFLO.run();	
	windFLO.write(outFile);		
}
