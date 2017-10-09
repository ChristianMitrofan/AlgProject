#include "Lmatrices.h"
#include "GenericMatrix.h"
#include "GenericMatrix.cpp"

using Eigen::MatrixXd;

using namespace std;

Lmatrices::Lmatrices(int m,SylvPol<GenericMatrix<double> > * pol,int d)
{
	size=m;
	int i,j;
	//Create the d coefficient of the Sylvester polynomial
	MatrixXd Md(m,m);
	for(i=0; i<m; i++)
	{
		for(j=0; j<m; j++)
		{
			Md(i,j) = pol->GetMatrixI(d)->GetMatrix(i,j);
		}
	}
	//Create L0
	L0 = MatrixXd::Zero(d*m,d*m);
	MatrixXd Diag(m*(d-1),m*(d-1));
	Diag << MatrixXd::Identity((d-1)*m,(d-1)*m);
	Diag=(-1)*Diag;
	L0.topRightCorner((d-1)*m,(d-1)*m)=Diag;
	MatrixXd Bot(m,m*d);						//Matrix that has all the Sylvester Coefficients together
	for(int c=0; c<d; c++)	
	{
		for(i=0; i<m; i++)
		{
			for(j=0; j<m; j++)
			{
				Bot(i,j+c*m) = pol->GetMatrixI(c)->GetMatrix(i,j);
			}
		}
	}
	L0.bottomRows(m)=Bot;
	//Create L1
	L1 = MatrixXd::Identity(d*m,d*m);
	L1.bottomRightCorner(m,m)=Md;
	L1=(-1)*L1;
	return;
}

void Lmatrices :: SetL0Matrix(int i,int j,double elem)
{
	L0(i,j)=elem;
}

void Lmatrices :: SetL1Matrix(int i,int j,double elem)
{
	L1(i,j)=elem;
}

void Lmatrices :: PrintL0Matrix()
{
	cout << L0 <<endl;
} 

void Lmatrices :: PrintL1Matrix()
{
	cout << L1 <<endl;
}

double Lmatrices :: GetL0Matrix(int i,int j)
{
	return L0(i,j);
}

double Lmatrices :: GetL1Matrix(int i,int j)
{
	return L1(i,j);
}

MatrixXd Lmatrices :: GetL0Matrix()
{
	return L0;
}

MatrixXd Lmatrices :: GetL1Matrix()
{
	return L1;
}

int Lmatrices :: Size() 
{
	return size;
}