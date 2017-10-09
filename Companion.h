#ifndef COMPANION_H
#define COMPANION_H

#include "GenericMatrix.h"
#include <iostream>
#include <cstdlib>
#include <Eigen/Core>
#include <Eigen/Dense>

using Eigen::MatrixXd; 

class Companion
{
	private:
		int size;
		MatrixXd CM;
	
	public:
		Companion(int m,SylvPol<GenericMatrix<double> > * pol,int d);
		//~Companion();
		void PrintMatrix();
		void SetMatrix(int r,int c,double element);
		double GetMatrix(int r,int c);
		MatrixXd GetMatrix();
		int Size();
};

#endif
