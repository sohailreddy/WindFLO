#ifndef MATH_HEADER
#define MATH_HEADER

/*
 * Created on Thu Jul 18 2019
 *
 * Copyright (c) 2019 Sohail R. Reddy
 * All rights reserved
 */

#include <iostream>
#include <vector>
#include <math.h>
#include <algorithm>
#include <cmath>





// Vector scaling and dot product operation
template<class Type>
std::vector<Type> operator*(const std::vector<Type> A, const std::vector<Type> B) {
    if(A.size() != B.size()){
        std::cout << " vector mismatch for * operator" << std::endl;
    }
    std::vector<Type> C;
    C.resize(A.size());
    for (size_t i = 0; i < A.size(); i++)
        C[i] = A[i] * B[i];
    return C;
}

template<class Type>
std::vector<Type> operator*(const std::vector<Type> A, const Type B) {
    std::vector<Type> C;
    C.resize(A.size());
    for (size_t i = 0; i < A.size(); i++)
        C[i] = A[i] * B;
    return C;
}

template<class Type>
std::vector<Type> operator*(const Type B, const std::vector<Type> A) {
    std::vector<Type> C;
    C.resize(A.size());
    for (size_t i = 0; i < A.size(); i++)
        C[i] = A[i] * B;
    return C;
} 

// Vector scaling and dot product operation
template<class Type>
std::vector<Type> operator/(const std::vector<Type> A, const std::vector<Type> B) {
    if(A.size() != B.size()){
        std::cout << " vector mismatch for / operator" << std::endl;
    }
    std::vector<Type> C;
    C.resize(A.size());
    for (size_t i = 0; i < A.size(); i++)
        C[i] = A[i] / B[i];
    return C;
}

template<class Type>
std::vector<Type> operator/(const std::vector<Type> A, const Type B) {
    std::vector<Type> C;
    C.resize(A.size());
    for (size_t i = 0; i < A.size(); i++)
        C[i] = A[i] / B;
    return C;
}

template<class Type>
std::vector<Type> operator/(const Type B, const std::vector<Type> A) {
    std::vector<Type> C;
    C.resize(A.size());
    for (size_t i = 0; i < A.size(); i++)
        C[i] = B /A[i];
    return C;
}

// Vector addition operation
template<class Type>
std::vector<Type> operator+(const std::vector<Type> A, const std::vector<Type> B) {
    if(A.size() != B.size()){
        std::cout << " vector mismatch for + operator" << std::endl;
    }
    std::vector<Type> C;
    C.resize(A.size());
    for (size_t i = 0; i < A.size(); i++)
        C[i] = A[i] + B[i];
    return C;
}

template<class Type>
std::vector<Type> operator+(const std::vector<Type> A, const Type B) {
    std::vector<Type> C;
    C.resize(A.size());
    for (size_t i = 0; i < A.size(); i++)
        C[i] = A[i] + B;
    return C;
}


template<class Type>
std::vector<Type> operator+(const Type B, const std::vector<Type> A) {
    std::vector<Type> C;
    C.resize(A.size());
    for (size_t i = 0; i < A.size(); i++)
        C[i] = B + A[i];
    return C;
}

// Vector subtraction operation
template<class Type>
std::vector<Type> operator-(const std::vector<Type> A, const std::vector<Type> B) {
    if(A.size() != B.size()){
        std::cout << " vector mismatch for - operator" << std::endl;
    }
    std::vector<Type> C;
    C.resize(A.size());
    for (size_t i = 0; i < A.size(); i++)
        C[i] = A[i] - B[i];
    return C;
}

template<class Type>
std::vector<Type> operator-(const std::vector<Type> A, const Type B) {
    std::vector<Type> C;
    C.resize(A.size());
    for (size_t i = 0; i < A.size(); i++)
        C[i] = A[i] - B;
    return C;
}

template<class Type>
std::vector<Type> operator-(const Type B, const std::vector<Type> A) {
    std::vector<Type> C;
    C.resize(A.size());
    for (size_t i = 0; i < A.size(); i++)
        C[i] = B - A[i];
    return C;
}


// Vector exponent operation
template<class Type>
std::vector<Type> operator^(const std::vector<Type> A, const std::vector<Type> B) {
    if(A.size() != B.size()){
        std::cout << " vector mismatch for ^ operator" << std::endl;
    }
    std::vector<Type> C;
    C.resize(A.size());
    for (size_t i = 0; i < A.size(); i++)
        C[i] = std::pow(A[i],B[i]);
    return C;
}

template<class Type>
std::vector<Type> operator^(const std::vector<Type> A, const Type B) {
    std::vector<Type> C;
    C.resize(A.size());
    for (size_t i = 0; i < A.size(); i++)
        C[i] = std::pow(A[i],B);
    return C;
}

 template<class Type>
std::vector<Type> operator^(const Type B, const std::vector<Type> A) {
    std::vector<Type> C;
    C.resize(A.size());
    for (size_t i = 0; i < A.size(); i++)
        C[i] = std::pow(B,A[i]);
    return C;
}


template<class Type>
std::vector<Type> abs(const std::vector<Type> A){
    size_t n = A.size();
    std::vector<Type> B(n);
    for (size_t i = 0; i < n; i++)
        B[i] = std::abs(A[i]);
    return B;
}

template<class Type>
Type sum(const std::vector<Type> A){
    Type total = 0.0;
    for (size_t i = 0; i < A.size(); i++)
        total += A[i];
    return total;
}

template<class Type>
Type norm(const std::vector<Type> A, double p=2.0){
    Type total = 0.0;
    for (size_t i = 0; i < A.size(); i++)
        total += std::pow(A[i],p);
    return std::pow(total,1.0/p);
}

template<class Type>
Type dotProduct(const std::vector<Type> A, const std::vector<Type> B){
    Type total = 0.0;
    for (size_t i = 0; i < A.size(); i++)
        total += A[i] * B[i];
    return total;
}

template<class Type>
Type perpDist(const std::vector<Type> A, const std::vector<Type> B, const std::vector<Type> C = {}){

    std::vector<Type> origin(A.size(),0.0);
    if(C.size() != 0)
        std::vector<Type> origin = C;

    std::vector<Type> W = A - origin;
    std::vector<Type> V = B - origin;

    Type AA = dotProduct(V,W);
    Type BB = dotProduct(V,V);
    Type CC = AA/BB;


    std::vector<Type> D = origin + (CC*V);
    std::vector<Type> X = A - D;

    Type dist = std::pow( dotProduct(X,X) ,0.5);
    return dist;
}

template<class Type>
Type calcAngle(const std::vector<Type> A, const std::vector<Type> B){

    Type V = std::pow(dotProduct(A,A),0.5);
    Type W = std::pow(dotProduct(B,B),0.5);
    Type dist = dotProduct(A,B)/ (V*W);

    if(dist > 1.0 || dist < -1.0)
        dist = 1.0;
    Type angle = acos(dist);

    return angle;
}

#endif /* MATH_HEADER */