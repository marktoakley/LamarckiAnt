/*
 * Particle.h
 *
 *  Created on: 21/set/2010
 *      Author: lorenzo
 */

#ifndef PARTICLE_H_
#define PARTICLE_H_

#include <cstring>
#include <cstdlib>
#include <cassert>

#include <stdio.h>

#include "defs.h"
#include "ExternalForce.h"

template<typename number>
class Particle {
protected:
	int _N_neigh;
	int *_verlet_list;
	int _max_neigh;
	int _current_neigh_index;
	LR_vector<number> _principal_axis;
	LR_vector<number> _stack_axis;

	void _check();
	void _set_back_position() {pos_back = orientation*_principal_axis*POS_BACK;}
	void _set_stack_position() {pos_stack = orientation*_principal_axis*POS_STACK;}
	void _set_base_position() {pos_base = orientation*_principal_axis*POS_BASE;}

public:
	int _N_ext_forces;
	ExternalForce<number> *_ext_forces[MAX_EXT_FORCES];
	int _next_particle;

	Particle();
	virtual ~Particle();

	void set_positions() {
		_set_back_position();
		_set_stack_position();
		_set_base_position();
	}

	void copy_from(const Particle<number> &);

	void init(int max_neigh);
	void is_neighbour(const Particle &p);
	// inline functions
	void prepare_list() { _current_neigh_index = 0; }
	void reset_lists() {
		_N_neigh = 0;
		pos_list = pos;
	}
	void add_neighbour(int nindex) {
		_verlet_list[_N_neigh++] = nindex;
		assert(_N_neigh < _max_neigh);
	}
	int next_neighbour() { return (_current_neigh_index < _N_neigh) ? _verlet_list[_current_neigh_index++] : P_VIRTUAL; }
	int get_index() const { return index; }

	// returns true if the force was added, false otherwise
	bool add_ext_force(ExternalForce<number> *f);
	void set_initial_forces(llint step) {
		this->force = LR_vector<number>(0, 0, 0);
		for(int i = 0; i < _N_ext_forces; i++) {
			this->force += _ext_forces[i]->value(step, this->pos);
		}
	}
	
	// returns true if the external potential was added, false otherwise
	bool add_ext_potential (ExternalForce<number> *f);
	void set_ext_potential (llint step) {
		this->ext_potential = (number) 0.;
		for(int i = 0; i < _N_ext_forces; i++) {
			this->ext_potential += _ext_forces[i]->potential(step, this->pos);
		}
	}

	int index;
	int type;
	int n3, n5;
	LR_vector<number> pos_list;
	LR_vector<number> torque;
	LR_vector<number> force;
	number ext_potential;
	LR_vector<number> pos;
	// these are the positions of backbone, base and stack sites relative to the center of mass
	LR_vector<number> pos_back;
	LR_vector<number> pos_stack;
	LR_vector<number> pos_base;
	LR_vector<number> L;
	LR_vector<number> vel;
	LR_matrix<number> orientation;
	// transpose (= inverse) orientational matrix
	LR_matrix<number> orientationT;
};

#endif /* PARTICLE_H_ */
