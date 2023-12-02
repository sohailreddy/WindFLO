/*

	Copyright 2019, Sohail R. Reddy
	email: sredd001@fiu.edu
	www.sohailreddy.com

*/

#ifndef CONVEXHULL_HEADER
#define CONVEXHULL_HEADER

#include <stack> 
#include <stdlib.h> 

struct Point 
{ 
    double x, y; 
}; 
  



std::vector<std::vector<double> > computeConvexHull(std::vector<std::vector<double> > &xy);
double ComputeConvexHullArea(std::vector<std::vector<double> > &xy);
std::vector<Point> findConvexHull(Point points[], int n);

#endif 