#ifndef SYLVESTER_H_
#define SYLVESTER_H_
#include "GenericMatrix.h"

void createSylvester(GenericMatrix<double> * pol1,GenericMatrix<double> * pol2,GenericMatrix<double *> * syl);

void createPolynomials(GenericMatrix<double *> * syl,SylvPol<GenericMatrix <double> > * sylp);

#endif
