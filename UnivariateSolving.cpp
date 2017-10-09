#include "UnivariateSolving.h"
#include "GenericMatrix.h"
#include "GenericMatrix.cpp"
#include <iostream>
#include <fstream>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::VectorXcd;
using Eigen::EigenSolver;
using Eigen::JacobiSVD;

using namespace std;

VectorXd Vectorize(GenericMatrix<double> * pol,double root,char hidden,int * s)	// Given a polynomial and the root that must be substituted it creates a mononomial
{
	int i,j,r,c;
	r=pol->GetRows();
	c=pol->GetColumns();
	MatrixXd Polynomial(r,c);
	VectorXd Mononomial;
	for(i=0; i<r; i++)
	{
		for(j=0; j<c; j++)
		{
			Polynomial(i,j) = pol->GetMatrix(i,j);
		}
	}
	cout << hidden << endl;
	cout << "The root is " << root << " and the polynomial is : " << endl;
	cout << Polynomial << endl;
	if(hidden=='y')
	{
		*s=r;
		Mononomial=VectorXd::Zero(c);		// All zeros
		for(i=1; i<r; i++)
		{
			for(j=0; j<c; j++)
			{
				Mononomial[j]=Polynomial(0,j)+pow(root,i)*Polynomial(i,j); // We substitute the y and put the result in the mononomial(Mononomial[0]=fixed term)
			}
		}
		for(i=0; i<c-1; i++)
		{
			Mononomial[i]=Mononomial[i]/Mononomial[c-1];	// So that the first term of the mononomial will always be one
		}
	}
	else if(hidden=='x')					//Same as y
	{
		*s=c;
		Mononomial=VectorXd::Zero(r);
		for(i=1; i<c; i++)
		{
			for(j=0; j<r; j++)
			{
				Mononomial[j]=Polynomial(0,j)+pow(root,i)*Polynomial(i,j);
			}
		}
		for(i=0; i<r-1; i++)
		{
			Mononomial[i]=Mononomial[i]/Mononomial[r-1];
		}
	}
	cout << "The mononomial is : " << endl;
	cout << Mononomial << endl;
	return Mononomial;
}

int FindCommonRoots(MatrixXd Cm1,MatrixXd Cm2,double root,char hvar)
{
	int c=0;
	int maxroots;
	VectorXcd vals1;
	VectorXcd vals2;
	vals1=StandardEigenProblemVals(Cm1);		//We find the eigenvalues of both mononomials
	vals2=StandardEigenProblemVals(Cm2);
	maxroots=vals1.size();						
	if(maxroots<vals2.size())
		maxroots=vals2.size();
	VectorXcd roots(maxroots);
	for(int i=0; i<vals1.size(); i++)			//Check if the roots match
	{
		for(int j=0; j<vals2.size(); j++)
		{
			if(vals1[i]==vals2[j])
			{
				roots[c]=vals1[i];
				c++;
			}
		}
	}
   /* ofstream myfile;
    myfile.open("roots.txt");
    for(int i=0;i<c;i++){
        if(hvar=='x'){
            myfile<< root << "  " << roots[i] << endl;
        }else{
            myfile<< roots[i] << "  " << root << endl;
        }
    }
    myfile.close();*/
	cout << "Eigenvalues : " << endl;
	cout << vals1 << endl;
	cout << "More eigenvalues :" << endl;
	cout << vals2 << endl;
	cout << "Common roots :" << endl;
	cout << roots << endl;
	return 1;
}


