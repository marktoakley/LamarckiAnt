/*
 * Utils.cpp
 *
 *  Created on: 04/set/2010
 *      Author: lorenzo
 */

#include "Utils.h"

Utils::Utils() {

}

Utils::~Utils() {

}

int Utils::decode_base(char c) {
	c = toupper(c);
	if(c == 'A') return N_A;
	if(c == 'C') return N_C;
	if(c == 'G') return N_G;
	if(c == 'T') return N_T;
	return P_VIRTUAL;
}

char Utils::encode_base(int b) {
	if(b == N_A) return 'A';
	if(b == N_C) return 'C';
	if(b == N_G) return 'G';
	if(b == N_T) return 'T';
	return 'X';
}

template float Utils::gaussian<float>();
template double Utils::gaussian<double>();
