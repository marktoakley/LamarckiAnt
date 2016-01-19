#include "GMINDNAManager.h"

// fake main function...
// as it is set now, it will
int main(int argc, char *argv[]) {
	
	int ndof, i, j;
	double * gino, * grad, energy;
   	
	initialise_(); // initialisation
	
	get_dof_(&ndof);
	//cout << "(main) ndof: " << ndof << endl ;;
	
	gino = (double *) malloc (ndof* sizeof (double));
	copy_coords_ (gino);
	
	grad = (double *) malloc (ndof* sizeof (double));
	
	potential_(gino, grad, &energy);
	
	// TEST: see if it works; it does after a long series of rhyming prophanities.	
	double * newc, e2, dc;
	newc = (double *) malloc (ndof * sizeof (double));
	if (argc < 2) {
		dc = 1.e-6;
	}
	else {
		dc = atof (argv[1]);
	}
	printf ("## testing the gradients with dc = %lf\n", dc);
	printf ("#d.o.f index - num. grad - grad - pert. energy - energy, dc\n");
	for (i = 0; i < ndof; i ++) {
		for (j = 0; j < ndof; j ++)
			newc[j] = gino[j];
		newc[i] += dc; // perturbation
		potential_(newc, grad, &e2);
		printf ("%5i %+10.7lf %+10.7lf %+10.7lf %+10.7lf %+10.7lf\n", i, (e2-energy)/double(dc), grad[i], e2, energy, dc);
	}
	// end of test
	
	free (newc);
	free (grad);
	free (gino);
		
    return 0;
}

