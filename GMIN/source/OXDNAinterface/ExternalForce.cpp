/*
 * ExternalForce.cpp
 *
 *  Created on: 18/lug/2011
 *      Author: lorenzo
 */

#include "ExternalForce.h"
#include <cstdio>
#include <cstdlib>

template<typename number>
ExternalForce<number>::ExternalForce(number base, number rate, LR_vector<number> dir) {
	_F0 = base;
	_rate = rate;
	_direction = dir;
}

template<typename number>
ExternalForce<number>::~ExternalForce() {

}

template class ExternalForce<double>;
template class ExternalForce<float>;

template<typename number>
ConstantRateForce<number>::ConstantRateForce(number base, number rate, LR_vector<number> dir) : ExternalForce<number>(base, rate, dir) {
}

template<typename number>
ConstantRateForce<number>::~ConstantRateForce() {

}

template<typename number>
LR_vector<number> ConstantRateForce<number>::value(llint step, LR_vector<number> &pos) {
	number x = (this->_F0 + this->_rate * step) * this->_direction.x;
	number y = (this->_F0 + this->_rate * step) * this->_direction.y;
	number z = (this->_F0 + this->_rate * step) * this->_direction.z;
	return LR_vector<number>(x, y, z);
}

template<typename number>
number ConstantRateForce<number>::potential(llint step, LR_vector<number> &pos) {
	return (number) -(this->_F0 + this->_rate * step) * (pos * this->_direction);
}

template class ConstantRateForce<double>;
template class ConstantRateForce<float>;

// added by Flavio a pene di segugio per Debayan Chakraborty
template<typename number>
MovingTrap<number>::MovingTrap(number stiff, LR_vector<number> pos0, number rate, LR_vector<number> dir) : ExternalForce<number>(0,rate,dir) {
    this->_stiff = stiff;
    this->_pos0 = pos0;
}

template<typename number>
LR_vector<number> MovingTrap<number>::value(llint step, LR_vector<number> &pos) {
    LR_vector<number> postrap;
    number x, y, z;

    postrap.x = this->_pos0.x + (this->_rate * step) * this-> _direction.x;
    postrap.y = this->_pos0.y + (this->_rate * step) * this-> _direction.y;
    postrap.z = this->_pos0.z + (this->_rate * step) * this-> _direction.z;

    x = - this->_stiff * (pos.x - postrap.x);
    y = - this->_stiff * (pos.y - postrap.y);
    z = - this->_stiff * (pos.z - postrap.z);
    
    return LR_vector<number>(x, y, z);
}

template<typename number>
number MovingTrap<number>::potential (llint step, LR_vector<number> &pos) {
    LR_vector<number> postrap;

    postrap.x = this->_pos0.x + (this->_rate * step) * this-> _direction.x;
    postrap.y = this->_pos0.y + (this->_rate * step) * this-> _direction.y;
    postrap.z = this->_pos0.z + (this->_rate * step) * this-> _direction.z;

    return (number) (0.5 * this->_stiff * ((pos - postrap) * (pos - postrap)));
}

template class MovingTrap<double>;
template class MovingTrap<float>;

