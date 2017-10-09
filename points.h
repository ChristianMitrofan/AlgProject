#ifndef POINTS_H_
#define POINTS_H_


#include <Eigen/Core>

using Eigen::MatrixXd;
using Eigen::VectorXd;

void checkPst(int * error,int rank ,int kpoints);

MatrixXd createInterpolation(double ** points,int d,int kpoints);

int computeRank(MatrixXd matrix);

VectorXd computeKernel(MatrixXd matrix);

#endif
