#ifndef PRTIMESLHOOD_I_INCLUDED
#define PRIMESTLHOOD_I_INCLUDED


#include <gsl/gsl_vector.h>

void abs_double(double x);
double pow_double(double x, double y);

void wrapper(double *volumeVect, double *ang0Vect,
                double *basFricAngVect, double *intAngVect, 
                double *volBasVect, double *y,
                double *circterms, double *parest, int npts,
                double *combo);
                

double prTimesLhood(double *volumeVect, double *ang0Vect,
                double *basFricAngVect, double *intAngVect,
                double *volBasVect, double *y, double *thetas, 
                double *circterms, double *parest, int npts);


double wrprTimesLhood(double *thetas, double *wrparams);
double gsl_wrprTimesLhood(const gsl_vector *thetasvect, void * wrparams);
void posteriorMode(double *thetas, double *wrparams);

#endif 
