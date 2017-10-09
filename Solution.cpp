#include "Solution.h"
#include "UnivariateSolving.h"
#include "VectorCompanion.h"
	
using namespace std;
using Eigen::EigenSolver;
using Eigen::MatrixXd;
using Eigen::MatrixXcd;
using Eigen::VectorXd;
using Eigen::VectorXcd;		
	
bool StandardProblem(SylvPol <GenericMatrix<double> > * pol,int size,int d,GenericMatrix<double> * pol1,GenericMatrix<double> * pol2,int * tz,char hvar)
{
	int c=0;
	int i,j;
	VectorXcd v;
	double deval;
	MatrixXcd Vecs;
	VectorXcd Vals;
	double  vx[size*d];
	double  vy[size*d];
	double  vm[size*d];
	MatrixXd Cm(size*d,size*d);
	
	//Creating Companion matrix;
	
	Companion * comp = new Companion(size,pol,d);
	Cm=comp->GetMatrix();
	Vals=StandardEigenProblemVals(Cm);
	Vecs=StandardEigenProblemVecs(Cm);
	//cout << "The eigenvalues of A are:" << endl << Vals << endl;
	//cout << "The matrix of eigenvectors, V, is:" << endl << Vecs << endl << endl;
	for(i=0; i<size*d; i++)
	{
		vm[i]=1;
	}
	for(i=0; i<size*d; i++)
	{
		v=Vecs.col(i);
		if(i!=0)
		{
			deval=Vals[i].real();
			for(j=0; j<i; j++)
			{
				if(deval==vy[j])
				{
					vm[j]++;
				}
					
			}
		}
        if((v[0].imag()==0) && (fabs(Vals[i].real())>pow(10,-5)))	//If the first value of the eigenvecture is not complex and the real part is not zero then we have a solution
		{
			vx[c]=v[1].real()/v[0].real();
			vy[c]=Vals[i].real();
			c++;
		}
        else if(Vals[i].imag()!=0){
            cout << Vals[i] <<"  "<<v[1]/v[0] <<" : Imaginary "<<endl;
        }
	}
	cout << "Roots" << endl << "-----" << endl;
	//Variables for multiple roots
	VectorXd Mononomial1;
	VectorXd Mononomial2;
	MatrixXd Cmatrix1;
	MatrixXd Cmatrix2;
	int csize;
	for(i=0; i<c; i++)
	{
		if(hvar=='y')
		{
			if(vm[i]>1)
			{
				cout << "y = " << vy[i] << "," << " multiplicity = " << vm[i] << endl;
				Mononomial1=Vectorize(pol1,vy[i],hvar,&csize);					// Create the first mononomial (from the first polynomial) and its companion matrix
				VectorCompanion * vcm1=new VectorCompanion(Mononomial1,csize);
				Cmatrix1=vcm1->GetMatrix();
				Mononomial2=Vectorize(pol1,vy[i],hvar,&csize);					// Create the second mononomial and its companion matrix
				VectorCompanion * vcm2=new VectorCompanion(Mononomial2,csize);
				Cmatrix2=vcm2->GetMatrix();
                FindCommonRoots(Cmatrix1,Cmatrix2,vy[i],'y');					// Find their common roots
			}
		}
		if(hvar=='x')
		{
			if(vm[i]>1)
			{
				cout << "x = " << vx[i] << "," << "multiplicity = " << vm[i] << endl;
				Mononomial1=Vectorize(pol1,vx[i],hvar,&csize);
				VectorCompanion * vcm1=new VectorCompanion(Mononomial1,csize);
				Cmatrix1=vcm1->GetMatrix();
				Mononomial2=Vectorize(pol1,vx[i],hvar,&csize);
				VectorCompanion * vcm2=new VectorCompanion(Mononomial2,csize);
				Cmatrix2=vcm2->GetMatrix();
                FindCommonRoots(Cmatrix1,Cmatrix2,vx[i],'x');
			}
		}
	}
   if(hvar == 'y')
        PolynomialRoots(pol1,pol2,vx,vy,vm,c,tz,hvar);
   else
        PolynomialRoots(pol1,pol2,vy,vx,vm,c,tz,hvar);
}
		
bool GeneralizedProblem(SylvPol <GenericMatrix<double> > * pol,int size,int d,GenericMatrix<double> * pol1,GenericMatrix<double> * pol2,int * tz,char hvar)
{	
	//cout << "Max power is : " << d << "(Over 9000)" << endl;
	//cout << "Size : " << size << endl;
	int c=0;
	int i,j;
	double deval;
	VectorXd evec;
	VectorXd eval;
	double vx[size*d];
	double vy[size*d];
	double vm[size*d];
	MatrixXd gl(size*d,3);
	MatrixXd gv(size*d,size*d);
	MatrixXd L0(size*d,size*d);
	MatrixXd L1(size*d,size*d);

	//Creating L0 and L1;
	
	Lmatrices * lm = new Lmatrices(size,pol,d);
	L0=lm->GetL0Matrix();
	L1=lm->GetL1Matrix();
	//cout << "Printing L0 : " << endl << L0 << endl << endl;
	//cout << "Printing L1 : " << endl << L1 << endl << endl;
	GeneralizedEigenProblem(L1,L0,gv,gl);
	//cout << "The eigenvalues of A are:" << endl << gl << endl << endl;
	//cout << "The matrix of eigenvectors, V, is:" << endl << gv << endl << endl;
	for(i=0; i<size*d; i++)
	{
		vm[i]=1;
	}
	for(i=0; i<size*d; i++)
	{
		eval=gl.row(i);
		evec=gv.col(i);
		if(i!=0)
		{
			deval=eval[2]/eval[0];
			for(j=0; j<i; j++)
			{
				if(deval==vy[j])
				{
					vm[j]++;
				}
					
			}
		}	
		if((eval[1]==0) && (fabs(eval[0])>pow(10,-10)) && (fabs(evec[0])>pow(10,-10)))	//If the first value of the eigenvecture is not complex and the real part is not zero then we have a solution	
		{
			vx[c]=evec[1]/evec[0];
			vy[c]=eval[2]/eval[0];
			c++;
		}
        else{
            //cout <<endl<< evec[1]/evec[0] <<" "<<eval[2]/eval[0] << ": Imaginary"<<endl;
        }
	}
	cout << "Roots" << endl << "-----" << endl;
	//Variables for multiple roots
	VectorXd Mononomial1;
	VectorXd Mononomial2;
	MatrixXd Cmatrix1;
	MatrixXd Cmatrix2;
	int csize;
	for(i=0; i<c; i++)
	{
		if(hvar=='y')
		{
			if(vm[i]>1)
			{
				cout << "y = " << vy[i] << "," << " multiplicity = " << vm[i] << endl;
				Mononomial1=Vectorize(pol1,vy[i],hvar,&csize);
				VectorCompanion * vcm1=new VectorCompanion(Mononomial1,csize);
				Cmatrix1=vcm1->GetMatrix();
				Mononomial2=Vectorize(pol1,vy[i],hvar,&csize);
				VectorCompanion * vcm2=new VectorCompanion(Mononomial2,csize);
				Cmatrix2=vcm2->GetMatrix();
                FindCommonRoots(Cmatrix1,Cmatrix2,vy[i],'y');
			}
		}
		if(hvar=='x')
		{
			if(vm[i]>1)
			{
				cout << "x = " << vx[i] << "," << "multiplicity = " << vm[i] << endl;
				Mononomial1=Vectorize(pol1,vx[i],hvar,&csize);
				VectorCompanion * vcm1=new VectorCompanion(Mononomial1,csize);
				Cmatrix1=vcm1->GetMatrix();
				Mononomial2=Vectorize(pol1,vx[i],hvar,&csize);
				VectorCompanion * vcm2=new VectorCompanion(Mononomial2,csize);
				Cmatrix2=vcm2->GetMatrix();
                FindCommonRoots(Cmatrix1,Cmatrix2,vx[i],'x');
			}
		}
    }
   if(hvar == 'y')
        PolynomialRoots(pol1,pol2,vx,vy,vm,c,tz,hvar);
   else
        PolynomialRoots(pol1,pol2,vy,vx,vm,c,tz,hvar);
}

