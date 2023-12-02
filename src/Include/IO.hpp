/*

	Copyright 2019, Sohail R. Reddy
	email: sredd001@fiu.edu
	www.sohailreddy.com

*/

#ifndef IO_HEADER
#define IO_HEADER

#include <vector>
#include <iostream>

// Print vector
template <typename Type> 
std::ostream& operator<<(std::ostream& os, const std::vector<Type>& A) 
{ 
    for (size_t i = 0; i < A.size(); i++) {
        os << " " << A[i] << " ";
    }
    return os; 
} 


#endif