/*
 * IOManager.h
 *
 *  Created on: 03/set/2010
 *      Author: lorenzo
 */

#ifndef IOMANAGER_H_
#define IOMANAGER_H_

#include <iostream>
#include <fstream>
// TODO: c-style, da trasformare in c++, se possibile
#include <cstdarg>
#include <cstdio>
#include <cstdlib>

extern "C" {
#include "parse_input.h"
}

#include "defs.h"
#include "MDBackend.h"
#include "MCBackend.h"

using namespace std;

class IOManager {
protected:
	bool _debug;

	// log part
	FILE *_log_stream;
	char _log_level_strings[5][256];
	bool _log_open;

	void _log(int log_level, const char *format, va_list &ap);
	void _set_stream(const char *filename);

	// simulation's I/O part
	FILE *_energy_stream;
	bool _no_stdout_energy;
	bool _restart_step_counter;
	char _trajectory_file[256];
	char _energy_file[256];
	char _reduced_conf_output_dir[256];
	char _last_conf_file[256];
	// prints out backend's current configuration.
	template<typename number> void _print_conf(SimBackend<number> &backend, llint curr_step, ofstream &conf_output, bool reduced=false);

public:
	enum { LOG_INFO = 0, LOG_WARNING = 1, LOG_ERROR = 2, LOG_CRITICAL = 3, LOG_DEBUG = 4, LOG_NOTHING = 5 };

	IOManager();
	virtual ~IOManager();

	void init();
	void get_settings(input_file &inp);
	void log(int log_level, const char *format, ...);
	void debug(const char *format, ...);
	void die(const char *format, ...);
	void print_version();

	FILE *get_log_stream() { return _log_stream; }

	template<typename number> void print_conf(SimBackend<number> &backend, llint curr_step, bool only_last);
	template<typename number> void print_reduced_conf(SimBackend<number> &backend, llint curr_step);

	template<typename number> void print_energy(MDBackend<number> &backend, llint curr_step);
	template<typename number> void print_energy(MCBackend<number> &backend, llint curr_step);
};

#endif /* IOMANAGER_H_ */
