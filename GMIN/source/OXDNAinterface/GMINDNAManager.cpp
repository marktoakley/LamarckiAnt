#include "GMINDNAManager.h"

bool PrintDataManager::stop = false;

PrintDataManager::PrintDataManager(IOManager &IO, const char *input_file, int confs_to_skip):
  _IO(IO), _print_reduced_conf_every(0), _print_energy_every(1000), _restart_step_counter(false) {
	_start_step = _cur_step = _steps = 0;
	_IO = IO;
	loadInputFile(&_input, input_file);
	if(_input.state == ERROR)
		_IO.die("Caught an error while opening the input file");
	_backend = new ProcessData_Backend<double>(&_IO,confs_to_skip);
}

ProcessData_Backend<double> *  PrintDataManager::get_backend_ptr () {
	return _backend;
}

PrintDataManager::~PrintDataManager() {
	cleanTimeScale(&_time_scale_manager);
	delete _backend;
	_IO.log(_IO.LOG_INFO, "FATTO...");
}

void PrintDataManager::_get_options() {
	getInputString(&_input, "conf_file", _conf_file, 1);
}

void PrintDataManager::load_options() {
	_get_options();
	_IO.get_settings(_input);
	_backend->get_settings(_input);
}

void PrintDataManager::init(const char* input_fname) {
	conf_input.open(input_fname);
	if(conf_input.good() == false)
	  _IO.die("Can't read configuration file '%s'", _conf_file);
	_backend->init(conf_input);
}

void PrintDataManager::run() {
	return;
}

void PrintDataManager::clean() {
	char unread[1024];
	unread[0] = '\0';
	// print unread (i.e. unused) keys
	setUnreadKeys(&_input);
	if(_input.N_unread_keys > 0) {
		for(int i = 0; i < _input.N_unread_keys; i++) {
			if(strlen(unread) < 1024)
				sprintf(unread+strlen(unread), "\n\t%s", _input.unread_keys[i]);
		}
		_IO.debug("The following keys found in the input file were not used: %s", unread);
	}
	cleanInputFile(&_input);
}

ostream& PrintDataManager::OutputState(ostream& out, int output_t)
{
	_backend->OutputBonds(out,_iteration);
	return out;
}
using namespace std;

#ifdef __cplusplus
extern "C" {
#endif

//static PrintDataManager * mysimptr;

void userpot_init_()
{
	//cout << "### initialise () " << endl;;
	// here we initialise the system from standard variables
    static IOManager IO;
    
	static PrintDataManager mysim(IO, "gmindnainput", 0);

    mysim.load_options();
    
    mysim.init("gmindnainitconf");

	mysimptr = &mysim;

	return;
}

void
userpot_get_natoms_ (int * natoms)
{
	ProcessData_Backend<double> * mp;
	mp = mysimptr->get_backend_ptr ();
	//* dof = (3 + 3) * (mp->get_N ());
 	* natoms = 2 * (mp->get_N ());
	return;
}

// sets coordinates from the array
void userpot_initialize_gmin_(int *ndof, double * coords)
{
// TODDO: vr274 check ndof == mp->get_N
	int i, N;
	ProcessData_Backend<double> * mp;
	mp = (ProcessData_Backend<double> *) mysimptr->get_backend_ptr ();
		
	N = mp->get_N ();
	Particle <double> * mols;
	mols = mp->get_particles_array();
	
	// coordinates
	for (i = 0; i < N; i ++ ) {
		coords[3 * i] = mols[i].pos.x;
		coords[3 * i + 1] = mols[i].pos.y;
		coords[3 * i + 2] = mols[i].pos.z;
	}
	
	// angles; we have to set things to angle-axis
	for (i = 0; i < N; i ++ ) {
		double ux, uy, uz, alpha, r, t;
		ux = mols[i].orientation.v3.y - mols[i].orientation.v2.z;
		uy = mols[i].orientation.v1.z - mols[i].orientation.v3.x;
		uz = mols[i].orientation.v2.x - mols[i].orientation.v1.y;
		r = sqrt (ux * ux + uy * uy + uz * uz);
		t = mols[i].orientation.v1.x + mols[i].orientation.v2.y + mols[i].orientation.v3.z;
		alpha = atan2 (r, t - 1.);
		if (alpha < 0) {
			// should not get here, but just in case...
			alpha = fabs (alpha); 
			ux = - ux;
			uy = - uy;
			uz = - uz;
		}
			
		coords[3 * N + 3 * i] = ux * alpha / r;
		coords[3 * N + 3 * i + 1] = uy * alpha / r;
		coords[3 * N + 3 * i + 2] = uz * alpha / r;
	}
	
	return;
}

// gets coordinates from the array
void calc_potential_(double * coords, double * grad, double * energy)
{
	int i, N;
	Particle<double> * mols;
	ProcessData_Backend<double> * mp;
	
	mp = mysimptr->get_backend_ptr ();
	N = mp->get_N ();
	mols = mp->get_particles_array ();
	
	for (i = 0; i < N; i ++ ) {
		mols[i].pos.x = coords[3 * i];
		mols[i].pos.y = coords[3 * i + 1];
		mols[i].pos.z = coords[3 * i + 2];
	}
	
	for (i = 0; i < N; i ++ ) {
		double p[3], theta, st, ct, ux, uy, uz;
		p[0] = coords[3 * N + 3 * i];
		p[1] = coords[3 * N + 3 * i + 1];
		p[2] = coords[3 * N + 3 * i + 2];
		theta = sqrt (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]); // angle
		
		// now translate those to the orientation matrix
		// http://en.wikipedia.org/wiki/Rotation_representation_(mathematics)
		ux = p[0] / theta; uy = p[1] / theta; uz = p[2] / theta;
		st = sin (theta); ct = cos(theta);
		LR_matrix<double> Ris;
		if (theta > 1.e-12) {
			// theta not extremely small
			Ris = LR_matrix<double> (ct + ux * ux * (1. - ct), ux * uy * (1. - ct) - uz * st, ux * uz * (1. - ct) + uy * st,
								     uy * ux * (1. - ct) + uz * st, ct + uy * uy * (1. - ct), uy * uz * (1. - ct) - ux * st,
                                     uz * ux * (1. - ct) - uy * st, uz * uy * (1. - ct) + ux * st, ct + uz * uz * (1. - ct));
		}
		else {
			// small theta... remember that p[i] = theta * u[i];
			Ris = LR_matrix<double> (1., -p[2], p[1],
									 p[2], 1., -p[0],
									 -p[1], p[0], 1.);
		}
		
		mols[i].orientation = Ris;
	}
	
	// set the forces and torques to 0. on all particles 
	for (i = 0; i < N; i ++) {
		mols[i].force.x = mols[i].force.y = mols[i].force.z = 0.;
		mols[i].torque.x = mols[i].torque.y = mols[i].torque.z = 0.;
	}
		
	// from the coordinates, it should fill gradients and energy;
	mp->OutputBonds(cout, 0);
	
	// now fill the arrays grad with the gradients and int with energy;
	* energy = mp->get_U ();
	for (i = 0; i < N; i ++) {
		grad[3 * i] =  - mols[i].force.x;
		grad[3 * i + 1] = - mols[i].force.y;
		grad[3 * i + 2] = - mols[i].force.z;
		
		// following Victor and Dwaipaian
		// we do not handle small theta yet...
		LR_matrix<double> R, dRx, dRy, dRz, pgn, I3, dpgndx, dpgndy, dpgndz;
		double px, py, pz, a, sa, ca;
		px = coords[3 * N + 3 * i];
		py = coords[3 * N + 3 * i + 1];
		pz = coords[3 * N + 3 * i + 2];
		
		a = sqrt (px * px + py * py + pz * pz);
		sa = sin(a); ca = cos(a);
		
		// we have to explicitely handle the case theta << 1
		if (a > 1.e-12) { // general case
			I3 = LR_matrix<double> (1, 0, 0, 0, 1, 0, 0, 0, 1);
		
			pgn = LR_matrix<double> (0, -pz/a,  py/a,
									 pz/a,  0, -px/a,
        	                        -py/a,  px/a, 0.);
		
			R = I3 + (pgn * pgn) * (1. - ca) + pgn * sa;
			
			dpgndx = LR_matrix<double> (0, px*pz/(a*a*a), -px*py/(a*a*a),
	                                   -px*pz/(a*a*a), 0, -1./a+px*px/(a*a*a),
										px*py/(a*a*a), 1./a-px*px/(a*a*a), 0);
			
			dpgndy = LR_matrix<double> (0, py*pz/(a*a*a), -py*py/(a*a*a)+1./a,
   		                                -py*pz/(a*a*a), 0, +px*py/(a*a*a),
										 py*py/(a*a*a)-1/a, -px*py/(a*a*a), 0);
			
			dpgndz = LR_matrix<double> (0, pz*pz/(a*a*a)-1./a, -pz*py/(a*a*a),
            	                        -pz*pz/(a*a*a)+1./a, 0, +px*pz/(a*a*a),
										+py*pz/(a*a*a), -px*pz/(a*a*a), 0);
		
			dRx = pgn * pgn * (px * sa / a)	+ (dpgndx * pgn + pgn * dpgndx) * (1. - ca) + pgn * (px * ca / a) + dpgndx * sa;
			dRy = pgn * pgn * (py * sa / a)	+ (dpgndy * pgn + pgn * dpgndy) * (1. - ca) + pgn * (py * ca / a) + dpgndy * sa;
			dRz = pgn * pgn * (pz * sa / a)	+ (dpgndz * pgn + pgn * dpgndz) * (1. - ca) + pgn * (pz * ca / a) + dpgndz * sa;
		
			// calcoliamo la matrice Rprime che avremmo per un pkprime
			/* TEST FOR DERIVATIVES OF MATICES
			double dc = 1.e-4, aprime;
			pz += dc;
			aprime = sqrt (px * px + py * py + pz * pz);
			LR_matrix<double> pgnprime = LR_matrix<double> (0, -pz/aprime,  py/aprime,
								 pz/aprime,  0, -px/aprime,
                                -py/aprime,  px/aprime, 0.);
		
			LR_matrix<double> Rprime = I3 + (pgnprime * pgnprime) * (1. - cos(aprime)) + pgnprime * sin(aprime);
			LR_matrix<double> dRR = (Rprime - R)/dc;
			printf("###############################\n");
			printf("## % 10.5lf  % 10.5lf  % 10.5lf\n", dRR.v1.x, dRR.v1.y, dRR.v1.z);
			printf("## % 10.5lf  % 10.5lf  % 10.5lf\n", dRR.v2.x, dRR.v2.y, dRR.v2.z);
			printf("## % 10.5lf  % 10.5lf  % 10.5lf\n\n", dRR.v3.x, dRR.v3.y, dRR.v3.z);
			dRR = dRz;
			printf("## % 10.5lf  % 10.5lf  % 10.5lf\n", dRR.v1.x, dRR.v1.y, dRR.v1.z);
			printf("## % 10.5lf  % 10.5lf  % 10.5lf\n", dRR.v2.x, dRR.v2.y, dRR.v2.z);
			printf("## % 10.5lf  % 10.5lf  % 10.5lf\n", dRR.v3.x, dRR.v3.y, dRR.v3.z);
			printf("###############################\n");
			pz -= dc;
			*/
			
			// now we find the matrix corresponding to -p
			px = -px; py = -py; pz = -pz;
			I3 = LR_matrix<double> (1, 0, 0, 0, 1, 0, 0, 0, 1);
			pgn = LR_matrix<double> (0, -pz/a,  py/a,
									 pz/a,  0, -px/a,
	                                -py/a,  px/a, 0.);
			R = I3 + pgn * pgn * (1. - ca) + pgn * sa; // rotation matrix
		}
		else {
			// case where the rotation angle is very small;
			// we set the rotation matrix and its derivatives explicitely
			R = LR_matrix<double> (1, -pz, py, pz, 1, -px, -py, px, 1);
			dRx = LR_matrix<double> (0, 0, 0, 0, 0,-1, 0, 1, 0);
			dRy = LR_matrix<double> (0, 0, 1, 0, 0, 0,-1, 0, 0);
			dRz = LR_matrix<double> (0,-1, 0, 1, 0, 0, 0, 0, 0);
		}
		
		// multiply
		LR_matrix<double> DRMI1 = dRx * R;
		LR_matrix<double> DRMI2 = dRy * R;
		LR_matrix<double> DRMI3 = dRz * R;
		
		// all true in working version
		/*
		assert(fabs(DRMI1.v3.y + DRMI1.v2.z) < 1.e-6);
		assert(fabs(DRMI1.v1.z + DRMI1.v3.x) < 1.e-6);
		assert(fabs(DRMI1.v2.x + DRMI1.v1.y) < 1.e-6);
		
		assert(fabs(DRMI2.v3.y + DRMI2.v2.z) < 1.e-6);
		assert(fabs(DRMI2.v1.z + DRMI2.v3.x) < 1.e-6);
		assert(fabs(DRMI2.v2.x + DRMI2.v1.y) < 1.e-6);
		
		assert(fabs(DRMI3.v3.y + DRMI3.v2.z) < 1.e-6);
		assert(fabs(DRMI3.v1.z + DRMI3.v3.x) < 1.e-6);
		assert(fabs(DRMI3.v2.x + DRMI3.v1.y) < 1.e-6);
		*/
		
		// we put back the torque in the lab frame; the -1 can come from 
		// anywhere, possibly the definition of the torque as -grad(U)
		LR_vector<double> labtorque = mols[i].orientation * mols[i].torque * (-1.);
		
		grad[3 * N + 3 * i] = 0.5 * (labtorque.x * (DRMI1.v3.y - DRMI1.v2.z) +  
								     labtorque.y * (DRMI1.v1.z - DRMI1.v3.x) +  
								     labtorque.z * (DRMI1.v2.x - DRMI1.v1.y));  

		grad[3 * N + 3 * i + 1] = 0.5 * (labtorque.x * (DRMI2.v3.y - DRMI2.v2.z) +  
								         labtorque.y * (DRMI2.v1.z - DRMI2.v3.x) +  
								         labtorque.z * (DRMI2.v2.x - DRMI2.v1.y));  
		
		grad[3 * N + 3 * i + 2] = 0.5 * (labtorque.x * (DRMI3.v3.y - DRMI3.v2.z) +  
								         labtorque.y * (DRMI3.v1.z - DRMI3.v3.x) +  
								         labtorque.z * (DRMI3.v2.x - DRMI3.v1.y));  
	}
	
	return;
}

// gets coordinates from the array
void userpot_dump_configuration_ (char *str, double * coords, long int strlen)
{
    int i, N;
    Particle<double> * mols;
    ProcessData_Backend<double> * mp;
    ofstream conf_output;
    double U=0; // todo add potential here
    double K=0; // kinetic energy

    char fname[512];
    strncpy(fname, str, strlen);
    fname[strlen] = 0;
    conf_output.open(fname);
    
    mp = mysimptr->get_backend_ptr ();
    N = mp->get_N ();
    mols = mp->get_particles_array ();
    
    for (i = 0; i < N; i ++ ) {
        mols[i].pos.x = coords[3 * i];
        mols[i].pos.y = coords[3 * i + 1];
        mols[i].pos.z = coords[3 * i + 2];
    }
    
    for (i = 0; i < N; i ++ ) {
        double p[3], theta, st, ct, ux, uy, uz;
        p[0] = coords[3 * N + 3 * i];
        p[1] = coords[3 * N + 3 * i + 1];
        p[2] = coords[3 * N + 3 * i + 2];
        theta = sqrt (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]); // angle
        
        // now translate those to the orientation matrix
        // http://en.wikipedia.org/wiki/Rotation_representation_(mathematics)
        ux = p[0] / theta; uy = p[1] / theta; uz = p[2] / theta;
        st = sin (theta); ct = cos(theta);
        LR_matrix<double> Ris;
        if (theta > 1.e-12) {
            // theta not extremely small
            Ris = LR_matrix<double> (ct + ux * ux * (1. - ct), ux * uy * (1. - ct) - uz * st, ux * uz * (1. - ct) + uy * st,
                                     uy * ux * (1. - ct) + uz * st, ct + uy * uy * (1. - ct), uy * uz * (1. - ct) - ux * st,
                                     uz * ux * (1. - ct) - uy * st, uz * uy * (1. - ct) + ux * st, ct + uz * uz * (1. - ct));
        }
        else {
            // small theta... remember that p[i] = theta * u[i];
            Ris = LR_matrix<double> (1., -p[2], p[1],
                                     p[2], 1., -p[0],
                                     -p[1], p[0], 1.);
        }
        
        mols[i].orientation = Ris;
    }
    
   
      conf_output.precision(15);

      conf_output << "t = " << 0 << endl; // potentially this can be used to label configurations
      conf_output << "b = " << mp->get_box_side() << " " << mp->get_box_side() << " " << mp->get_box_side() << endl;
      conf_output << "E = " << K/N + U/N << " " << U/N << " " << K/N  << endl;

      for(int i = 0; i < N; i++) {
              LR_matrix<double> oT = mols[i].orientation.get_transpose();
              conf_output << mols[i].pos.x << " " << mols[i].pos.y << " " << mols[i].pos.z << " ";
              conf_output << oT.v1.x << " " << oT.v1.y << " " << oT.v1.z << " ";
              conf_output << oT.v3.x << " " << oT.v3.y << " " << oT.v3.z << " ";
              conf_output << mols[i].vel.x << " " << mols[i].vel.y << " " << mols[i].vel.z << " ";
              conf_output << mols[i].L.x << " " << mols[i].L.y << " " << mols[i].L.z << endl;
      }
      conf_output.close();
}
#ifdef __cplusplus
}
#endif

