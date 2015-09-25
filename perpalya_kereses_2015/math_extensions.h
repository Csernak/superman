//============================================================================
// Name        : cycle.cpp
// Author      : Gergely Gyebrószki
// Version     :
// Copyright   : 2014 (C) gyebro.com
// Description : Some extensions to math.h
//============================================================================

#ifndef MATH_EXTENSIONS_H
#define MATH_EXTENSIONS_H

#include <math.h>
#include <vector>
#include <iostream>

using namespace std;

template <class T>
vector<size_t> divisors(T number) {
	size_t max = (size_t)number;
	vector<size_t> divs;
	for (size_t i=1; i<max; i++ ){
		if(number%i == 0) {
			divs.push_back(i);
		}
	}
	return divs;
}

unsigned long long int uintpow(const unsigned long long int num, const size_t pow) {
	unsigned long long int result = num;
	for(size_t i=1; i<pow; i++) {
		result*=num;
	}
	return result;
}

#endif

