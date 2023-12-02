/*

	Copyright 2019, Sohail R. Reddy
	email: sredd001@fiu.edu
	www.sohailreddy.com

*/


#include "turbine.hpp"
#include "input.hpp"
#include "convexHull.hpp"


std::vector<std::vector<double> > computeConvexHull(std::vector<std::vector<double> > &xy){

	int nPoints = xy.size();
	Point *points;
	
	points = new Point[nPoints];
	for(int i = 0; i < nPoints ; i++){
		points[i].x = xy[i][0];	points[i].y = xy[i][1];
	}
	
	std::vector<Point> bndPoints = findConvexHull(points, nPoints); 
	delete points;
	
	std::vector<std::vector<double> > boundPoints( bndPoints.size());	
	for(int i = 0; i < bndPoints.size() ; i++){
		boundPoints[i].resize(2);
		boundPoints[i][0] = bndPoints[i].x ;
		boundPoints[i][1] = bndPoints[i].y ;
	}
	
	return boundPoints;
}

double ComputeConvexHullArea(std::vector<std::vector<double> > &xy){

	double area = 0.0;
  	
  	int n = xy.size();
    // Calculate value of shoelace formula 
    int j = n - 1; 
    for (int i = 0; i < n; i++)  { 
        area += (xy[j][0] + xy[i][0]) * (xy[j][1] - xy[i][1]); 
        j = i;  // j is previous vertex to i 
    } 
  
    return abs(area / 2.0); 
}


Point p0; //used to another two points
Point secondTop(std::stack<Point> &stk) {
   Point tempPoint = stk.top(); 
   stk.pop();
   Point res = stk.top();    //get the second top element
   stk.push(tempPoint);      //push previous top again
   return res;
}

double squaredDist(Point p1, Point p2) {
   return ((p1.x-p2.x)*(p1.x-p2.x) + (p1.y-p2.y)*(p1.y-p2.y));
}

int direction(Point a, Point b, Point c) {
   double val = (b.y-a.y)*(c.x-b.x)-(b.x-a.x)*(c.y-b.y);
   if (val == 0.0)
      return 0;    //colinear
   else if(val < 0.0)
      return 2;    //anti-clockwise direction
      return 1;    //clockwise direction
}

int comp(const void *point1, const void*point2) {
   Point *p1 = (Point*)point1;
   Point *p2 = (Point*)point2;
   int dir = direction(p0, *p1, *p2);
   if(dir == 0)
      return (squaredDist(p0, *p2) >= squaredDist(p0, *p1))?-1 : 1;
   return (dir==2)? -1 : 1;
}

std::vector<Point> findConvexHull(Point points[], int n) {
   std::vector<Point> convexHullPoints;
      
   double minY = points[0].y;
   int min = 0;
   for(int i = 1; i<n; i++) {
      double y = points[i].y;
      //find bottom most or left most point
      if((y < minY) || ((minY == y) && points[i].x < points[min].x)) {
         minY = points[i].y;
         min = i;
      }
   } 
   
   std::swap(points[0], points[min]);    //swap min point to 0th location
   p0 = points[0];
   qsort(&points[1], n-1, sizeof(Point), comp);    //sort points from 1 place to end
   int arrSize = 1;    //used to locate items in modified array
   for(int i = 1; i<n; i++) {
      //when the angle of ith and (i+1)th elements are same, remove points
      while(i < n-1 && direction(p0, points[i], points[i+1]) == 0)
         i++;
         points[arrSize] = points[i];
         arrSize++;
   }
   if(arrSize < 3)
      return convexHullPoints;    //there must be at least 3 points, return empty list.
      //create a stack and add first three points in the stack
            
      std::stack<Point> stk;
      stk.push(points[0]); stk.push(points[1]); stk.push(points[2]);

   for(int i = 3; i<arrSize; i++) {    //for remaining vertices
      while(direction(secondTop(stk), stk.top(), points[i]) != 2)
         stk.pop();    //when top, second top and ith point are not making left turn, remove point
         stk.push(points[i]);
   }
   
   while(!stk.empty()) {
      convexHullPoints.push_back(stk.top());    //add points from stack
      stk.pop();
   }
   return convexHullPoints;
}  
  
  
  
  
