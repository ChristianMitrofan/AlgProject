#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "GenericMatrix.cpp"
#include "GenericMatrix.h"
#include "Multiply.h"


using Eigen::MatrixXd;

void multiply(GenericMatrix<double *> * Sylv,GenericMatrix<double *> * res){		//(given matrix,calculated matrix)
	int vectorsize=Sylv->GetColumns();
	GenericMatrix<double> * vector = new GenericMatrix<double>(vectorsize,1);		//the vector to be given
	double el;
	printf("please give a vector with %d elements \n",Sylv->GetColumns());	//inserted by the user
	for(int i=0;i<vectorsize;i++){
		scanf("%lf",&el);
		vector->SetMatrix(i,0,el);						//set the vector parts
	}
	int h=Sylv->GetHigh();
	double result[h];
	double sum[h];
	double element[h];
	for(int i=0;i<vectorsize;i++){				//matrix multiplication
		for(int k1=0;k1<h;k1++)	
			sum[k1]=0;
		for(int j=0;j<vectorsize;j++){
			for(int kh=0;kh<h;kh++){
				element [kh] = Sylv->GetMatrix(i,j,kh);
				result[kh]=0;
			}
			for(int h1=0;h1<h;h1++){
				result[h1] = element[h1]*vector->GetMatrix(j,0);			
				sum[h1]+=result[h1];					//calculations between same powers
				}
			}
		res->SetMatrix(i,0,sum);
	}
	delete vector;
}


double findFactorial(int n){ //factorial of an int number
	if ((n==1)|| (n==0))
		return 1;
	else
		return n*findFactorial(n-1);
}

double findBinCoeff(int n,int k){	//Calculation of Binomial Coeff
	double n1=findFactorial(n);
	double k1=findFactorial(k);
	double nk=findFactorial(n-k);
	//printf("n1 %d k1 %d ,nk %d \n",n1,k1,nk);
	return n1/(nk*k1);
}

double * factorialarr(int t1,int t2,int d){ // returns polynomial^power in open form
	int coeff,c,c1,c2;
	double * farr = (double *)malloc(sizeof(double)*(d+1));
	for(int k=0;k<=d;k++)
		farr[k]=0;
	if(d==0){
		farr[0]=1;
		return farr;	
	}
	for(int k=0;k<=d;k++){
		c=findBinCoeff(d,k);
		c1=pow(t1,k);
		c2=pow(t2,d-k);
		coeff=c1*c2*c;
		//printf("k %d,c %d,c1 %d,c2 %d,coeff %d\n",k,c,c1,c2,coeff);
		farr[k]+=coeff;
	}
	return farr;
}

double * polmult(double * pol1,double * pol2,int d1,int d2){	//Calculates the multiplication of two polynomials in ooen form
	int coeff,power;
	double * result = (double *)malloc(sizeof(double)*(d1+d2+1));
	for(int i=0;i<=d1+d2;i++)
		result[i]=0.0;
	for(int i=0;i<=d1;i++)
		for(int j=0;j<=d2;j++){
			coeff=pol1[i]*pol2[j];
			power=i+j;
			result[power]+=coeff;
		}
	//printf("HEREEEEE\n");
return result;
	
}

 GenericMatrix<double> * MandCmult(GenericMatrix<double> * pol,int c){	//Calculates and returns the int*Matrix
	int s=pol->GetRows();
	GenericMatrix<double> * polres = new GenericMatrix<double>(s,s);
	for(int i=0;i<s;i++)
		for(int j=0;j<s;j++)
			polres->SetMatrix(i,j,c*pol->GetMatrix(i,j));
	return polres;
}

 GenericMatrix<double> * MandMsum(GenericMatrix<double> * m1,GenericMatrix<double> * m2){ //Calculates and returns the Matrix1+Matrix2
	int s = m1->GetRows();
	GenericMatrix<double> * polres = new GenericMatrix<double>(s,s);
	for(int i=0;i<s;i++)
		for(int j=0;j<s;j++)
			polres->SetMatrix(i,j,m1->GetMatrix(i,j)+m2->GetMatrix(i,j));
	//printf("polres\n");
	//polres->PrintMatrix();
	return polres;
}
void CreateZpol(SylvPol<GenericMatrix <double> > * sylp,SylvPol<GenericMatrix <double> > * sylpZ,int t1,int t2,int t3,int t4){ //Calculates the Z polynomial
	cout << "Calculating Z pol with "<<" t1 :"<<t1<<" t2 :"<<t2<<" t3 :"<<t3<<" t4 :"<<t4<<endl;
	double * pol,*polz2,*polz1,*mult;
	int d,npower,k;
	int dmax=sylp->GetMatrixI(0)->GetRows();
 	MatrixXd Md(dmax,dmax);
 	MatrixXd Md2(dmax,dmax);
	GenericMatrix <double> * Mmult;
	d=sylp->GetSizeI();
	SylvPol<GenericMatrix <double> > ** nsylps = new SylvPol<GenericMatrix <double> > *[d];
	for(int i=0;i<d;i++){				//each matrix Md,d=0...d is Multiplicated with polynomials (t1*x+t2)^d*(t3*x+t4)^d-i
		polz1 = factorialarr(t1,t2,d-i-1);
		polz2 = factorialarr(t3,t4,i);
		mult = polmult(polz1,polz2,d-i-1,i);
		nsylps[i] = new SylvPol<GenericMatrix <double> > (d);
		Mmult =sylp->GetMatrixI(d-1-i);//h i
		for(int j=d-1;j>=0;j--){
			for(int i2=0;i2<dmax;i2++)
 			{
 	 			for(int j2=0;j2<dmax;j2++)
 				{
  				 Md(i2,j2) = Mmult->GetMatrix(i2,j2);
  				}
 			}
			Md=Md*mult[j];
			GenericMatrix <double> * Ms=new GenericMatrix <double> (dmax,dmax);
			for(int i1=0;i1<dmax;i1++)
 			{
 	 			for(int j1=0;j1<dmax;j1++)
 				{
  					Ms->SetMatrix(i1,j1,Md(i1,j1));
  				}
 			}
			nsylps[i]->SetMatrixI(Ms,j);
		}
		delete polz1;
		delete polz2;
		delete mult;
	}
	for(int l=0;l<d;l++){
		for(int i=0;i<dmax;i++)
 		{
 	 		for(int j=0;j<dmax;j++)
 			{
  				 Md(i,j) = nsylps[0]->GetMatrixI(l)->GetMatrix(i,j);
  			}
 		}
		//sums = nsylps[0]->GetMatrixI(l);
		for(int m=1;m<d;m++){
			for(int i=0;i<dmax;i++)
 			{
 	 			for(int j=0;j<dmax;j++)
 				{
  					 Md2(i,j) = nsylps[m]->GetMatrixI(l)->GetMatrix(i,j);
  				}
 			}
			Md+=Md2;
			//sums=MandMsum(sums,nsylps[m]->GetMatrixI(l));
		}
		GenericMatrix<double> * sums = new GenericMatrix<double>(dmax,dmax);
		for(int i=0;i<dmax;i++)
 		{
 	 		for(int j=0;j<dmax;j++)
 			{
  				sums->SetMatrix(i,j,Md(i,j));
  			}
 		}
		sylpZ->SetMatrixI(sums,l);
	}
	//delete sums;
	for(int i=0;i<d;i++)
		delete nsylps[i];  
	delete nsylps;
	delete Mmult;
	//sylpZ->PrintMatrix();
}
