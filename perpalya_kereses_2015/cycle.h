//============================================================================
// Name        : cycle.cpp
// Author      : Gergely Gyebroszki
// Version     :
// Copyright   : 2014 (C) gyebro.com
// Description : Prime Cycle Generator
//============================================================================

#ifndef CYCLE_H
#define CYCLE_H

#include <vector>
#include <iostream>
/*#include "math_extensions.h"
#include "armadillo"
using namespace arma;
*/
using namespace std;

namespace MicroChaosTools {

/**
 * Template Cycle class for symbolic dynamics
 *  template parameter should be the smallest type covering the desired symbols
 *   (i.e. byte, unsigned int, unsigned long, etc)
 *  can permute itself
 *  can increase itself
 *  does not check for invalid values
 *  alphabet goes from 0 to N-1 in case of symb = N
 */
template < class T >
class Cycle {
private:
	vector<T> elem; //egy ciklus
	size_t symb; 	/**< Number of different symbols */
	size_t len;		/**< Length of cycle */
public:
	Cycle(size_t alphabet, size_t length) {
		symb = alphabet;
		len = length;
		elem.resize(len); 
	}
	Cycle(size_t alphabet, vector<T> elements) {
		symb = alphabet;
		elem = elements;
		len = elements.size();
	}
	T const& operator() (size_t i) const {
		return elem[i];
	}
	T& operator() (size_t i) {
		return elem[i];
	}
	void fill(T value) {
		for (size_t i=0; i<len; i++) {
			elem[i] = value;
		}
	}
	bool increase(void) {
		size_t depth = len-1;
		bool working = true;
		bool successful = false;
		size_t big = symb-1;
		while(working) {
			if (elem[depth] < big) {
				elem[depth]+=1;
				working = false;
				successful = true;
			}
			else {
				elem[depth] = 0;
				depth-=1;
				if (depth == -1) {
					working = false;
					successful = false;
				}
			}
		}
		return successful;
	}
	Cycle<T> permuteLeft(void) {
		Cycle<T> result(symb,len);
		result(len-1) = elem[0];
		for (size_t i=0; i<len-1; i++) {
			result(i) = elem[i+1];
		}
		return result;
	}
	Cycle<T> permuteRight(void) {
		Cycle<T> result(symb,len);
		result(0) = elem[len-1];
		for (size_t i=1; i<len; i++) {
			result(i) = elem[i-1];
		}
		return result;
	}
	size_t length(void) const {
		return len;
	}
	
	/*vector<T> cvec(void) {
	vector<T> vektor;
	for (size_t i=1; i<len; i++) {
			vektor[i] = elem[i];
	}
	  return vektor;
	}
	*/
	bool isHighestCyclicPermutation(void) {
		bool highest = false;
		bool working = true;
		Cycle<T> perm(symb,len);
		perm = (*this);
		size_t pcount = 0;
		// Generate the other len-1 cyclic permutations
		while(working) {
			// Check permutation count
			if (pcount == len-1) {
				working = false;
				highest = true;
			}
			pcount++;
			perm = perm.permuteLeft();
			// Check digits
			for (size_t i=0; i<len; i++) {
				if(elem[i] > perm(i)) {
					// Stop checking lower permutation, go on next
					break;
				}
				else if(elem[i] < perm(i)){
					// The current permutation is higher than the this
					working = false;
					break;
				}
				else {
					// If last digit is equal then stop
					if (i == len-1) {
						// The current permutation is equal to the origin
						working = false;
						break;
					}
				}
			}
		}
		return highest;
	}
	bool isLowestCyclicPermutation(void) {
		bool lowest = false;
		bool working = true;
		Cycle<T> perm(symb,len);
		perm = (*this);
		size_t pcount = 0;
		// Generate the other len-1 cyclic permutations
		while(working) {
			// Check permutation count
			if (pcount == len-1) {
				working = false;
				lowest = true;
			}
			pcount++;
			perm = perm.permuteLeft();
			// Check digits
			for (size_t i=0; i<len; i++) {
				if(elem[i] < perm(i)) {
					// Stop checking higher permutation, go on next
					break;
				}
				else if(elem[i] > perm(i)){
					// The current permutation is lower than the this
					working = false;
					break;
				}
				else {
					// If last digit is equal then stop
					if (i == len-1) {
						// The current permutation is equal to the origin
						working = false;
						break;
					}
				}
			}
		}
		return lowest;
	}
};

template < class T >
ostream& operator << (ostream& out, Cycle<T> const& c) {
	size_t length = c.length();
	for (size_t i=0; i<length-1; i++) {
		out << (size_t)c(i) << " ";
	}
	out << (size_t)c(length-1) << "\n";
	return out;
}

template < class T >
ostream& operator << (ostream& out, vector<T> const& c) {
	size_t length = c.size();
	for (size_t i=0; i<length-1; i++) {
		out << c[i] << " ";
	}
	out << c[length-1] << "\n";
	return out;
}

} /*end of namespace MicroChaosTools */

#endif
