/*
 * ExternalForce.h
 *
 *  Created on: 18/lug/2011
 *      Author: lorenzo
 */

#ifndef EXTERNALFORCE_H_
#define EXTERNALFORCE_H_

#include "defs.h"

template<typename number>
class ExternalForce {
public:
	// we need these members to be public because 
	// we need access in order to copy these numbers
	// to the GPU memory
	number _rate;
	number _F0;
	LR_vector<number> _direction;
	LR_vector<number> _pos0;
	number _stiff;
	
	ExternalForce(number base, number rate, LR_vector<number> dir);
	virtual ~ExternalForce();
	
	virtual LR_vector<number> value(llint step, LR_vector<number> &pos) = 0;
	virtual number potential (llint step, LR_vector<number> &pos) = 0;
};

template<typename number>
class ConstantRateForce : public ExternalForce<number> {
public:
	ConstantRateForce(number base, number rate, LR_vector<number> dir);
	virtual ~ConstantRateForce();
	
	virtual LR_vector<number> value(llint step, LR_vector<number> &pos);
	virtual number potential (llint step, LR_vector<number> &pos);

};

template<typename number>
class MovingTrap : public ExternalForce<number> {
public:
	MovingTrap (number stiff, LR_vector<number> pos0, number rate, LR_vector<number> dir);
	virtual ~MovingTrap() {}
	
	virtual LR_vector<number> value(llint step, LR_vector<number> &pos);
	virtual number potential(llint step, LR_vector<number> &pos);
};

#endif /* EXTERNALFORCE_H_ */
