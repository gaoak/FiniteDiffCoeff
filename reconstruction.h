// fortran to c arguments
// integer    int*
// integer array, int *
extern "C" void recon_coef_u_(int *r, int *k, double *coeff);
extern "C" void recon_coef_a_(int *r, int *k,double * coef, double *stencil);
