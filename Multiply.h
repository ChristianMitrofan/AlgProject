#ifndef MULTIPLY_H_
#define MULTIPLY_H_

#include "GenericMatrix.h"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

void multiply(GenericMatrix<double *> * Sylv,GenericMatrix<double *> * res);

void CreateZpol(SylvPol<GenericMatrix <double> > * sylp,SylvPol<GenericMatrix <double> > * sylpZ,int t1,int t2,int t3,int t4);

double findBinCoeff(int n,int k);

double * factorialarr(int t1,int t2,int d);

double * polmult(double * pol1,double * pol2,int d1,int d2);



#endif
