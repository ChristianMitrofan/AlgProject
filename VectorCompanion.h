#ifndef VECTORCOMPANION_H
#define VECTORCOMPANION_H

#include "GenericMatrix.h"
#include <iostream>
#include <cstdlib>
#include <Eigen/Core>
#include <Eigen/Dense>

using Eigen::MatrixXd;
using Eigen::VectorXd; 

class VectorCompanion
{
	private:
		int size;
		MatrixXd VCM;
	
	public:
		VectorCompanion(VectorXd Mononomial,int s);
		//~VectorCompanion();
		void PrintMatrix();
		void SetMatrix(int r,int c,double element);
		double GetMatrix(int r,int c);
		MatrixXd GetMatrix();
		int Size();
};

#endif
