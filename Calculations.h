#ifndef CALCULATIONS_H_
#define CALCULATIONS_H_

#include "GenericMatrix.h"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <lapacke.h>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::VectorXcd;
using Eigen::MatrixXcd;

int detChecks(SylvPol <GenericMatrix<double> > * pol);

bool GeneralizedEigenProblem(MatrixXd & A,MatrixXd & B,MatrixXd & v,MatrixXd & lambda);

VectorXcd StandardEigenProblemVals(MatrixXd & A);

MatrixXcd StandardEigenProblemVecs(MatrixXd & A);

int computeK(SylvPol <GenericMatrix <double > > * sylp,int B,int height,double * K);

void PolynomialRoots(GenericMatrix<double> *pol1,GenericMatrix<double> *pol2,double * vx,double * vy,double * vm,int c,int * tz,char hidden);

#endif
