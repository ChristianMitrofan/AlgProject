#ifndef UNIVARIATESOLVING_H_
#define UNIVARIATESOLVING_H_

#include <cmath>
#include <iostream>
#include <cstdlib>
#include <Eigen/Core>
#include <Eigen/Dense>
#include "GenericMatrix.h"
#include "Calculations.h"

using Eigen::VectorXd;

VectorXd Vectorize(GenericMatrix<double> * pol,double root,char hidden,int * s);

int FindCommonRoots(MatrixXd Cm1, MatrixXd Cm2, double root, char hvar);

#endif
