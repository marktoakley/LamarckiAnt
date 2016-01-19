/*
 * Particle.cpp
 *
 *  Created on: 21/set/2010
 *      Author: lorenzo
 */

#include "Particle.h"

template<typename number>
Particle<number>::Particle() : _N_neigh(0), _verlet_list(NULL), _principal_axis(1,0,0), _stack_axis(0,1,0), _N_ext_forces(0), index(-1), type(P_VIRTUAL) {

}

template<typename number>
void Particle<number>::copy_from(const Particle<number> &p) {
	index = p.index;
	type = p.type;
	_max_neigh = p._max_neigh;
	_N_neigh = p._N_neigh;
	pos = p.pos;
	vel = p.vel;
	orientation = p.orientation;
	pos_list = p.pos_list;
	force = p.force;

	memcpy(_verlet_list, p._verlet_list, _N_neigh * sizeof(int));
}

template<typename number>
Particle<number>::~Particle() {
	if(_verlet_list != NULL) delete[] _verlet_list;
	for(int i = 0; i < _N_ext_forces; i++) delete _ext_forces[i];
}

template<typename number>
bool Particle<number>::add_ext_force(ExternalForce<number> *f) {
	if(_N_ext_forces == MAX_EXT_FORCES) return false;

	_ext_forces[_N_ext_forces] = f;
	_N_ext_forces++;

	return true;
}

template<typename number>
void Particle<number>::init(int max_neigh) {
	_max_neigh = max_neigh;
	_verlet_list = new int[_max_neigh];

	_check();
}

template<typename number>
void Particle<number>::_check() {
	assert(index >= 0);
	assert(type != P_VIRTUAL);
	assert(_verlet_list != NULL);
}

template class Particle<double>;
template class Particle<float>;
