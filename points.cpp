#include "points.h"
#include <stdio.h>
#include <iostream>
#include <cstdlib>
#include <Eigen/SVD>
#include <Eigen/LU>


using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::JacobiSVD;
//using Eigen::FullPivLU;

using namespace std;

void checkPst(int * error,int rank ,int kpoints){
	if(rank==kpoints){
		*error=0;
		cout << "\nthe system of equations is well-constrained \n" <<endl;
		return;
	}else if(rank<kpoints){
		*error=1;
		cout << "\nthe problem is infeasible or the solution is numerically unstable\n" <<endl;
		return;
	}/*else{
		*error=2;
		cout << "\nunderdefined, usually infinite number of solutions \n" <<endl;
		return;
	}*/	
}

MatrixXd createInterpolation(double ** points, int d, int kpoints){ // Interpolation matrix
    MatrixXd interp;
    interp=MatrixXd::Zero(kpoints,kpoints+1);
	for(int i=0;i<kpoints;i++){
		int maxpow=1;
		int powx=1;
		int powy=0;
		for(int j=0;j<=kpoints;j++){
			if(j==0){
				interp(i,0)=1;
			}else{
                //cout <<"maxpower = " << maxpow << " powx = " << powx << " powy = " << powy << endl;
                //cout <<"x is " << points[i][1] <<" ^ " << powx << " y is " << points[i][0]<<" ^ " << powy <<endl;
				if(powx==0)
					interp(i,j)=pow(points[i][1],powy);
				else if(powy==0)		
					interp(i,j)=pow(points[i][0],powx);
				else
					interp(i,j)=pow(points[i][0],powx)*pow(points[i][1],powy);
				//cout << temp(i,j) << endl;
				if(powy==maxpow){
					maxpow++;
					powx=maxpow;
					powy=0;
				}else{
					powx--;
					powy++;
				}				
			}						
		}
	}

    //cout << interp << endl;
    return interp;
}

int computeRank(MatrixXd  matrix){ // rank of a matrix
	int rank=0;
	int s =matrix.rows();
	JacobiSVD<MatrixXd> svd(matrix);							//Singular Values Calculation
	VectorXd sings = svd.singularValues();
	for(int i=0;i<s;i++)
		if (sings(i)>(pow(10,-7)))
			rank++;
	cout << rank << endl;
	return rank;

}

VectorXd computeKernel(MatrixXd matrix){ // kernel of a given matrix
    VectorXd cos;
    cos = VectorXd::Zero(matrix.cols());
	cos=matrix.fullPivLu().kernel();
    //cout << matrix <<endl;
    //cout << cos << endl;
    return cos;
}

