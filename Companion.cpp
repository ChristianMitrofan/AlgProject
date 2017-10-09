#include "Companion.h"
#include "GenericMatrix.h"
#include "GenericMatrix.cpp"

using Eigen::MatrixXd;

using namespace std;

Companion::Companion(int m,SylvPol<GenericMatrix<double> > * pol,int d)
{
	size=m;
	int i,j;
	cout << m << d << endl;
	//Create the d coefficient of the Sylvester polynomial
	MatrixXd Md(m,m);
	for(i=0;i<m;i++)
	{
		for(j=0;j<m;j++)
		{
			Md(i,j) = pol->GetMatrixI(d)->GetMatrix(i,j);
		}
	}
	MatrixXd Mdinv(m,m);
	Mdinv=Md.inverse();
	Mdinv=(-1)*Mdinv;
	cout<< Mdinv << endl;
	//Create the Companion matrix
	
	CM = MatrixXd::Zero(d*m,d*m);				//All zeros
	MatrixXd Diag(m*(d-1),m*(d-1));
	Diag << MatrixXd::Identity((d-1)*m,(d-1)*m);	
	CM.topRightCorner((d-1)*m,(d-1)*m)=Diag; 	//Create all the m*m identity matrices in the companion matrix
	MatrixXd Bot(m,m*d);						//Matrix that has all the Sylvester Coefficients together
	MatrixXd Mi(m,m);	
	for(int c=0; c<d; c++)					
	{
		for(i=0; i<m; i++)
		{
			for(j=0; j<m; j++)
			{
				Mi(i,j)=pol->GetMatrixI(c)->GetMatrix(i,j);
			}
		}
		cout << Mi << endl;
		Mi=Mdinv*Mi;
		cout << Mi << endl;
		for(i=0; i<m; i++)
		{
			for(j=0; j<m; j++)
			{
				Bot(i,j+c*m) = Mi(i,j);
			}
		}
	}
	CM.bottomRows(m)=Bot;
	cout << CM << endl;
	return;
}

void Companion :: SetMatrix(int i,int j,double elem)
{
	CM(i,j)=elem;
}

void Companion :: PrintMatrix()
{
	cout << CM <<endl;
} 

double Companion :: GetMatrix(int i,int j)
{
	return CM(i,j);
}

MatrixXd Companion :: GetMatrix()
{
	return CM;
}

int Companion :: Size() 
{
	return size;
}
