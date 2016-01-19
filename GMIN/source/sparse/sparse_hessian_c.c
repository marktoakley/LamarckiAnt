#include <cholmod.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

int quench_num;

void get_pointer(double ** hessian);
void free_pointer(double ** hessian);

void errors_to_file(int status, const char *file, int line, const char *message) {
   // Prints errors and warnings to an output file, with default name "cholmod_errors".

   // Open the file
   FILE * error_file;
   if (quench_num == 0) {
      error_file = fopen("cholmod_error", "w");
   }
   else {
      error_file = fopen("cholmod_error", "a");
   }

   // Write the error message to the file.
   fprintf(error_file, "Quench: %d\n", quench_num);
   fprintf(error_file, "CHOLMOD warning/error: %s\n", message);
   fprintf(error_file, "Encountered in file %s on line %d\n", file, line);

   // Close the file
   fclose(error_file);
}

cholmod_sparse *cholmod_hessian(double *hessian, int dim, cholmod_common *common) {
    // This function assigns the Hessian matrix from a dense matrix for CHOLMOD to use.
    // Allocate a cholmod triplet matrix of appropriate size
    cholmod_triplet *triplet_hessian;
    triplet_hessian = cholmod_allocate_triplet(dim, dim, dim*dim, 0, CHOLMOD_REAL, common);
    // Loop through values of hessian and assign their row/column index and values to triplet_hessian.
    int loop;
    for (loop = 0; loop < (dim * dim); loop++) {
        // If the hessian contains zeros, skip to the next item
        if (hessian[loop] == 0) continue;
        // Only add lower triangular entries
        if ((loop / dim) < (loop % dim)) continue;
        ((int*)triplet_hessian->i)[triplet_hessian->nnz] = loop / dim;
        ((int*)triplet_hessian->j)[triplet_hessian->nnz] = loop % dim;
        ((double*)triplet_hessian->x)[triplet_hessian->nnz] = hessian[loop];
        triplet_hessian->nnz++;
    }
    // Convert the triplet to a sparse matrix and return.
    cholmod_sparse *sparse_hessian;
    sparse_hessian = cholmod_triplet_to_sparse(triplet_hessian, (dim * dim), common);
    // Make sure that CHOLMOD knows that the hessian is lower triangular (i.e. stype = -1, other way around
    // from triplet stype)
    sparse_hessian->stype = -1;
    // Free the triplet
    cholmod_free_triplet(&triplet_hessian, common);

    return sparse_hessian;
}

double determinant(cholmod_sparse *sparse_hessian, int dim, cholmod_common *common) {
    // This function factorizes the sparse hessian matrix and returns the determinant.
    double det = 0.0;
    int i, j;
    cholmod_factor *factor;

    // Analyze and then factorize the sparse_hessian, storing the result in factor.
    factor = cholmod_analyze(sparse_hessian, common);
    cholmod_factorize(sparse_hessian, factor, common);

    // Change it to a simplicial LDL^T factorisation. This can also be done by setting
    // common->final_super = FALSE, common->final_ll = FALSE. I don't think it makes a
    // difference.
    cholmod_change_factor(CHOLMOD_REAL, 0, 0, 0, 0, factor, common);

    // Cast the column pointers and value pointers to the right types
    int *col_ptr = (int *)factor->p;
    int *row_idx = (int *)factor->i;
    int *non_zeros = (int *)factor->nz;
    double *values = (double *)factor->x;

    // Accumulate the log determinant
    for (i = 0; i < factor->n; i++) {
        det += log(fabs(values[col_ptr[i]]));
    }

    // Free used arrays
    cholmod_free_factor(&factor, common);

    return det;
}

void get_determinant(int *dim, double *det, int *quench) {
    double *hessian;
    cholmod_sparse *sparse_hess;
    cholmod_common common;
    
    // Start CHOLMOD
    cholmod_start(&common);

    // Turn off printing and use our own error handler
    (&common)->print = 0;
    (&common)->error_handler = &errors_to_file;
    quench_num = *quench;

    // Get the hessian from GMIN. The first element is at *hessian.
    get_pointer(&hessian);

    // Get a sparse copy of the hessian
    sparse_hess = cholmod_hessian(hessian, *dim, &common);

    // Calculate the determinant using that sparse hessian
    *det = determinant(sparse_hess, *dim, &common);

    // Free the hessian and point at the null pointer, just to be sure.
    free_pointer(&hessian); 
    hessian = 0;

    // Free CHOLMOD memory and finish
    cholmod_free_sparse(&sparse_hess, &common);
    cholmod_finish(&common);
}
