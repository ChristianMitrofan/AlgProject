#ifndef Lmatrices_H
#define Lmatrices_H

#include "GenericMatrix.h"
#include <iostream>
#include <cstdlib>
#include <Eigen/Core>
#include <Eigen/Dense>

using Eigen::MatrixXd; 

class Lmatrices
{
	private:
		int size;
		MatrixXd L0;
		MatrixXd L1;
	
	public:
		Lmatrices(int m,SylvPol<GenericMatrix<double> > * pol,int d);
		//~Lmatrices();
		void PrintL0Matrix();
		void PrintL1Matrix();
		void SetL0Matrix(int r,int c,double element);
		void SetL1Matrix(int r,int c,double element);
		double GetL0Matrix(int r,int c);
		double GetL1Matrix(int r,int c);
		MatrixXd GetL0Matrix();
		MatrixXd GetL1Matrix();
		int Size();
};

#endif
