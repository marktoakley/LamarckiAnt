/*
 * SimBackend.cpp
 *
 *  Created on: 03/set/2010
 *      Author: lorenzo
 */

#include "SimBackend.h"
#include "IOManager.h"

template<typename number>
SimBackend<number>::SimBackend(IOManager *IO) :
	_IO(IO), _timer_msgs_number(0), _print_timings(false), _external_forces(false), _are_lists_old(true), _N_updates(0), _confs_to_skip(0)  {
}

template<typename number>
SimBackend<number>::~SimBackend() {
	delete[] _heads;
	delete[] _particles;

	prepare_timer_results(&_timer);

	_IO->log(_IO->LOG_INFO, "Timings informations:\n");
	print_times(&_timer, _IO->get_log_stream());
	_IO->log(_IO->LOG_NOTHING, "");

	if(_print_timings == true) {
		FILE *timings_file = fopen(_timings_filename, "a");
		fprintf(timings_file, "%d %lf\n", _N, _timer.timings[0]);
		fclose(timings_file);
	}

	destroy_timer(&_timer);
}

template<>
void SimBackend<float>::_get_number_settings(input_file &inp) {
	getInputFloat(&inp, "verlet_skin", &_verlet_skin, 1);
}

template<>
void SimBackend<double>::_get_number_settings(input_file &inp) {
	getInputDouble(&inp, "verlet_skin", &_verlet_skin, 1);
}

template<typename number>
void SimBackend<number>::get_settings(input_file &inp) {
	int tmp;

	getInputString(&inp, "topology", _topology_filename, 1);
	getInputInt(&inp, "confs_to_skip", &_confs_to_skip, 0);

	if (getInputInt (&inp, "external_forces", &tmp, 0) == KEY_FOUND) {
		_external_forces = (tmp != 0);
		if (_external_forces) {
			getInputString (&inp, "external_forces_file", _external_filename, 1);
		}
	}

	if(getInputInt(&inp, "print_timings", &tmp, 0) == KEY_FOUND) {
		_print_timings = (tmp != 0);
		if(_print_timings) getInputString(&inp, "timings_filename", _timings_filename, 1);
	}

	char raw_T[256];
	getInputString(&inp, "T", raw_T, 1);
	char deg;
	double tmp_T;
	int res = sscanf(raw_T, "%lf %c", &tmp_T, &deg);
	if(res == 2) {
		deg = tolower(deg);
		switch(deg) {
		case 'c':
			_T = (tmp_T + 273.15) * 0.1 / 300.; // convert to kelvin and then to simulation units
			this->_IO->log(this->_IO->LOG_INFO, "Converting temperature from Celsius (%lf C°) to simulation units (%lf)", tmp_T, _T);
			break;
		case 'k':
			_T = tmp_T * 0.1 / 300.; // convert to simulation units
			this->_IO->log(this->_IO->LOG_INFO, "Converting temperature from Kelvin (%lf K) to simulation units (%lf)", tmp_T, _T);
			break;
		default:
			this->_IO->die("Unrecognizable temperature '%s'", raw_T);
		}
	}
	else _T = tmp_T;

	_get_number_settings(inp);
}

template<typename number>
void SimBackend<number>::_read_external_forces () {
	_IO->log (_IO->LOG_INFO, "Parsing Force file %s", _external_filename);

	//char line[512], typestr[512];
	int open, justopen, a;
	ifstream external (_external_filename);

	if(!external.good ()) _IO->die ("Can't read external_forces_file '%s'", _external_filename);

	justopen = open = 0;
	a = external.get();
	while(external.good()) {
		justopen = 0;
		if (a == '{') {
			open ++;
			justopen = 1;
		}
		if (a == '}') {
			if (justopen) _IO->die ("Syntax error in '%s': nothing between parentheses", _external_filename);
			open --;
		}
		if (open > 1 || open < 0) _IO->die ("Syntax error in '%s': parentheses do not match", _external_filename);
		a = external.get();
	}
	external.clear();
	external.seekg(0, ios::beg);

	FILE * temp;
	//int tempfd;
	//char templ[] = ".extXXXXXX";
	a = external.get();
	while(external.good()) {
		while (a != '{' && external.good()) {
			a = external.get();
		  }

		if(!external.good()) break;

		//tempfd = mkstemp (templ);
		//temp = fdopen (tempfd, "w");
		temp = fopen (".merda", "w");
		_IO->log (_IO->LOG_INFO, "   Using temporary file");
		a = external.get ();
		while (a != '}' && external.good()) {
			fprintf (temp, "%c", a);
			a = external.get ();
		}
		fclose (temp);

		// facciamo il parsing della forza
		temp = fopen (".merda", "r");
		input_file input;
		char type_str[512];
		int type;

		loadInput (&input, temp);
		getInputString (&input, "type", type_str, 1);

		type = -1;
		if (strcmp (type_str, "string") == 0) type = 0;
		if (strcmp (type_str, "twist") == 0) type = 1;
		if (strcmp (type_str, "trap") == 0) type = 2;
		if (type != 0 && type != 2) _IO->die ("force type %s not implemented. Aborting", type_str);

		switch (type) {
		case 0: {
			// forza a constant rate
			double rate, F0, dirx, diry, dirz;
			char strdir[512];
			int npart, tmpi;
			getInputDouble(&input, "F0", &F0, 1);
			getInputDouble(&input, "rate", &rate, 1);
			getInputInt(&input, "particle", &npart, 1);
			getInputString(&input, "dir", strdir, 1);
			tmpi = sscanf(strdir, "%lf,%lf,%lf", &dirx, &diry, &dirz);
			if (tmpi != 3)
				_IO->die("could not parse direction in external_forces_file. Dieing badly");
			_IO->log(_IO->LOG_INFO, "--> adding constant rate force f = (%g + %g * t) [%g,%g,%g] on part %i", F0, rate, dirx, diry, dirz, npart);
			// TODO: passare al valore il tempo in "secondi" e non quello
			// in passi, cosi' se uno cambia passo non muore gonfio
			ConstantRateForce<number> *extF = new ConstantRateForce<number>((number) F0, (number) rate, LR_vector<number>((number)dirx, (number)diry, (number)dirz));
			_particles[npart].add_ext_force(extF);
			break;
		}
		case 2: {
			double rate, pos0x, pos0y, pos0z, stiff, dirx, diry, dirz;
			char strdir[512], posdir[512];
			int npart, tmpi;
			getInputDouble(&input, "stiff", &stiff, 1);
			getInputDouble(&input, "rate", &rate, 1);
			getInputInt(&input, "particle", &npart, 1);
			getInputString(&input, "dir", strdir, 1);
			tmpi = sscanf(strdir, "%lf,%lf,%lf", &dirx, &diry, &dirz);
			if (tmpi != 3)
				_IO->die("could not parse dir in external_forces_file. Dieing badly");
			getInputString(&input, "pos0", posdir, 1);
			tmpi = sscanf(posdir, "%lf,%lf,%lf", &pos0x, &pos0y, &pos0z);
			if (tmpi != 3)
				_IO->die("could not parse pos0 in external_forces_file. Dieing badly");
			_IO->log(_IO->LOG_INFO, "--> adding MovingTrap with stiffnes %lf and pos=[%g,%g,%g] + (%g * t) [%g,%g,%g] on part %i", stiff, pos0x, pos0y, pos0z, rate, dirx, diry, dirz, npart);

			MovingTrap<number> *extF = new MovingTrap<number>((number) stiff, LR_vector<number> ((number)pos0x, (number) pos0y, (number) pos0z), (number) rate, LR_vector<number>((number)dirx, (number)diry, (number)dirz));
			_particles[npart].add_ext_force(extF);
			break;
		}
		default:
			_IO->log(_IO->LOG_INFO,	"Probably should't reach this point. Hoping for the best");
			break;
		}
		fclose (temp);
	}

	_IO->log (_IO->LOG_INFO, "   Force file parsed", _external_filename);
}

template<typename number>
void SimBackend<number>::_check_input_sanity() {
	for(int i = 0; i < _N; i++) {
		Particle<number> *p = &_particles[i];
		if(p->n3 >= _N) _IO->die("Wrong topology for particle %d (n3 neighbor is %d, should be < N = %d)", i, p->n3, _N);
		if(p->n5 >= _N) _IO->die("Wrong topology for particle %d (n5 neighbor is %d, should be < N = %d)", i, p->n5, _N);

		// check that the distance between bonded neighbor doesn't exceed a reasonable threshold
		number mind =  FENE_R0 - FENE_DELTA;
		number maxd =  FENE_R0 + FENE_DELTA;
		p->set_positions();
		if(_particles[i].n3 != P_VIRTUAL) {
			Particle<number> *q = &_particles[p->n3];
			q->set_positions();
			LR_vector<number> rv = p->pos + p->pos_back - (q->pos + q->pos_back);
			number r = sqrt(rv*rv);
			if(r > maxd || r < mind)
				_IO->die("Distance between bonded neighbors %d and %d exceeds acceptable values (d = %lf)", i, p->n3, r);
		}

		if(_particles[i].n5 != P_VIRTUAL) {
			Particle<number> *q = &_particles[p->n5];
			q->set_positions();
			LR_vector<number> rv = p->pos + p->pos_back - (q->pos + q->pos_back);
			number r = sqrt(rv*rv);
			if(r > maxd || r < mind)
				_IO->die("Distance between bonded neighbors %d and %d exceeds acceptable values (d = %lf)", i, p->n5, r);
		}
	}
}

template<typename number>
void SimBackend<number>::_read_topology() {
	char line[512];
	ifstream topology(_topology_filename);

	if(!topology.good()) _IO->die("Can't read topology file '%s'", _topology_filename);

	topology.getline(line, 512);
	sscanf(line, "%d %d\n", &_N, &_N_strands);

	_particles = new Particle<number>[_N];

	char base;
	int strand, i = 0;
	while(topology.good()) {
		topology.getline(line, 512);
		if(strlen(line) == 0 || line[0] == '#') continue;
		if(i == _N) _IO->die("Too many particles found in the topology file (should be %d), aborting", _N);

		int res = sscanf(line, "%d %c %d %d", &strand, &base, &_particles[i].n3, &_particles[i].n5);
		if(res < 4) _IO->die("Line %d of the topology file has an invalid syntax", i+2);

		if(_particles[i].n3 < 0) _particles[i].n3 = P_VIRTUAL;
		if(_particles[i].n5 < 0) _particles[i].n5 = P_VIRTUAL;

		_particles[i].type = Utils::decode_base(base);
		if(_particles[i].type == P_VIRTUAL) _IO->die("Particle #%d in strand #%d contains a non valid base '%c'", i, strand, base);
		_particles[i].index = i;
		i++;
	}

	if(i < _N) _IO->die("Not enough particles found in the topology file (should be %d), aborting", _N);

	topology.close();
}

template<typename number>
void SimBackend<number>::init(ifstream &conf_input) {
	init_timer(&_timer, _timer_msgs_number, _timer_msgs);
	_interaction.init(_T);

	_read_topology();

	if (_external_forces) _read_external_forces ();

	conf_input.seekg(0);
	char line[512];
	conf_input.getline(line, 512);
	sscanf(line, "t = %*d");
	conf_input.getline(line, 512);
	double box_side;
	sscanf(line, "b = %lf %*f %*f", &box_side);
	_box_side = (number) box_side;
	conf_input.getline(line, 512);
	sscanf(line, "E = %*f %*f %*f");

	_sqr_verlet_skin = SQR(_verlet_skin);
	// we choose rcut as the max of the range interaction of excluded volume between backbones and
	// hydrogen bonding
	number rcutback = 2 * fabs(POS_BACK) + 2 * EXCL_RC1;
	number rcutbase = 2 * fabs(POS_BASE) + 2 * HYDR_RCHIGH;
	_rcut = fmax(rcutback, rcutbase);
	_sqr_rcut = SQR(_rcut);
	_sqr_rverlet = SQR(_rcut + _verlet_skin * 2);

	_N_cells_side = (int) floor(this->_box_side / sqrt(_sqr_rverlet));
	while(_N_cells_side > ceil(pow(2*_N, 1/3.)) && _N_cells_side > 3) {
		_N_cells_side--;
	}

	if(_N_cells_side < 3) _IO->die("N_cells_side (%d) must be > 2", _N_cells_side);

	_N_cells = _N_cells_side*_N_cells_side*_N_cells_side;
	_heads = new int[_N_cells];

	_IO->log(_IO->LOG_INFO, "N_cells: %d, N_cells_side: %d", _N_cells, _N_cells_side);

	// if we are performing a CUDA simulation then we don't need CPU verlet lists
	if(_is_CUDA_sim == false) {
		// volume of a sphere whose radius is ceil(rverlet / (smallest sigma)) times maximum hard sphere density (sqrt(2)).
		_max_neigh = (int) ceil((4 * M_PI * pow(ceil(sqrt(_sqr_rverlet)), 3) / 3.) * sqrt(2.));
		if(_max_neigh > _N) _max_neigh = _N-1;
	}
	else _max_neigh = 0;

	// we need to skip a certain number of lines, depending on how many particles we have and how many
	// configurations we want to be skipped
	if(_confs_to_skip > 0) {
		char trash[5000];
		int rows_to_skip = _confs_to_skip * (_N+3);
		for(int i = 0; i < rows_to_skip; i++) conf_input.getline(trash, 5000);
	}

	int i = 0;
	while(!conf_input.eof() && i < _N) {
		conf_input >> _particles[i].pos.x >> _particles[i].pos.y >> _particles[i].pos.z;
		conf_input >> _particles[i].orientation.v1.x >> _particles[i].orientation.v1.y >> _particles[i].orientation.v1.z;
		conf_input >> _particles[i].orientation.v3.x >> _particles[i].orientation.v3.y >> _particles[i].orientation.v3.z;
		_particles[i].orientation.v1.normalize();
		_particles[i].orientation.v3.normalize();
		// orthonormalization
		_particles[i].orientation.v1 -= _particles[i].orientation.v3 * (_particles[i].orientation.v1*_particles[i].orientation.v3);
		_particles[i].orientation.v1.normalize();
		_particles[i].orientation.v2 = _particles[i].orientation.v3.cross(_particles[i].orientation.v1);
		_particles[i].orientation.transpone();

		conf_input >> _particles[i].vel.x >> _particles[i].vel.y >> _particles[i].vel.z;
		conf_input >> _particles[i].L.x >> _particles[i].L.y >> _particles[i].L.z;

		_particles[i].init(_max_neigh);

		i++;
	}

	if(i != this->_N) {
		if(_confs_to_skip > 0) _IO->die("Wrong number of particles (%d) found in configuration. Maybe you skipped too many configurations?", i);
		else _IO->die("The number of lines found in configuration file (%d) doesn't match the parsed number of particles (%d)", i, this->_N);
	}

	if(_is_CUDA_sim == false) _IO->log(_IO->LOG_INFO, "N: %d, max_neigh: %d", _N, _max_neigh);
	else _IO->log(_IO->LOG_INFO, "N: %d", _N);

	this->_check_input_sanity();
}

template<>
void SimBackend<float>::_fill_cell_index(const LR_vector<float> &pos, int cell_index[3]) {
	cell_index[0] = (int) ((pos.x / _box_side - floor(pos.x / _box_side)) * (1.f - FLT_EPSILON) * _N_cells_side);
	cell_index[1] = (int) ((pos.y / _box_side - floor(pos.y / _box_side)) * (1.f - FLT_EPSILON) * _N_cells_side);
	cell_index[2] = (int) ((pos.z / _box_side - floor(pos.z / _box_side)) * (1.f - FLT_EPSILON) * _N_cells_side);
}

template<>
void SimBackend<double>::_fill_cell_index(const LR_vector<double> &pos, int cell_index[3]) {
	cell_index[0] = (int) ((pos.x / _box_side - floor(pos.x / _box_side)) * (1. - DBL_EPSILON) * _N_cells_side);
	cell_index[1] = (int) ((pos.y / _box_side - floor(pos.y / _box_side)) * (1. - DBL_EPSILON) * _N_cells_side);
	cell_index[2] = (int) ((pos.z / _box_side - floor(pos.z / _box_side)) * (1. - DBL_EPSILON) * _N_cells_side);
}

template<typename number>
void SimBackend<number>::_create_cells() {
	for(int i = 0; i < _N_cells; i++) _heads[i] = P_VIRTUAL;

	int ind[3];
	for(int i = 0; i < this->_N; i++) {
		Particle<number> *p = &this->_particles[i];
		_fill_cell_index(p->pos, ind);

		int cell_index = (ind[0] * _N_cells_side + ind[1]) * _N_cells_side + ind[2];
		int old_head = _heads[cell_index];
		_heads[cell_index] = i;
		p->_next_particle = old_head;
	}
}

template<typename number>
void SimBackend<number>::_update_lists() {
	Particle<number> *p, *q;

	_create_cells();

	int ind[3], loop_ind[3];
	for(int i = 0; i < this->_N; i++) {
		p = &this->_particles[i];
		p->reset_lists();

		_fill_cell_index(p->pos, ind);
		for(int j = -1; j < 2; j++) {
			loop_ind[0] = (ind[0] + j + _N_cells_side) % _N_cells_side;
			for(int k = -1; k < 2; k++) {
				loop_ind[1] = (ind[1] + k + _N_cells_side) % _N_cells_side;
				for(int l = -1; l < 2; l++) {
					loop_ind[2] = (ind[2] + l + _N_cells_side) % _N_cells_side;
					int loop_index = (loop_ind[0] * _N_cells_side + loop_ind[1]) * _N_cells_side + loop_ind[2];

					int qindex = _heads[loop_index];
					while(qindex != P_VIRTUAL) {
						q = &this->_particles[qindex];
						// if this is an MC simulation then we need full lists, otherwise i-th particle will have neighbours with index > i
						if((i < qindex || (_sim_type == SIM_MC && i != qindex)) && p->pos.minimum_image(q->pos, this->_box_side).norm() < _sqr_rverlet)
							if(p->n3 != qindex && p->n5 != qindex) p->add_neighbour(qindex);
						qindex = q->_next_particle;
					}
				}
			}
		}
	}

	_are_lists_old = false;
}

template class SimBackend<float>;
template class SimBackend<double>;

