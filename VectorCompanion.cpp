#include "VectorCompanion.h"
#include "GenericMatrix.h"
#include "GenericMatrix.cpp"

using Eigen::MatrixXd;
using Eigen::VectorXd;

using namespace std;

VectorCompanion::VectorCompanion(VectorXd Mononomial,int s)
{
	size=s;
	cout << "The size of the companion is : " << s << endl;
	VCM = MatrixXd::Zero(s,s);				//All zeros
	MatrixXd Diag(s-1,s-1);
	Diag << MatrixXd::Identity(s-1,s-1);	
	VCM.topRightCorner(s-1,s-1)=Diag;
	for(int i=0; i<s; i++)
	{
		VCM(s-1,i)=Mononomial[i];
	}
	cout << "The companion matrix is : " << endl << VCM << endl;
}

void VectorCompanion :: SetMatrix(int i,int j,double elem)
{
	VCM(i,j)=elem;
}

void VectorCompanion :: PrintMatrix()
{
	cout << VCM <<endl;
} 

double VectorCompanion :: GetMatrix(int i,int j)
{
	return VCM(i,j);
}

MatrixXd VectorCompanion :: GetMatrix()
{
	return VCM;
}

int VectorCompanion :: Size() 
{
	return size;
}
