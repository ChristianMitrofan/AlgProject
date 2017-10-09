#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctime>
#include <cmath>
#include "Calculations.h"
#include "Multiply.h"
#include "Polynomials.h"
#include "Sylvester.h"
#include "Solution.h"
#include "GenericMatrix.cpp"
#include "GenericMatrix.h"
#include "Companion.h"
#include "Lmatrices.h"
#include "points.h"


#define BUFFSIZE 200					//size to read(each line read must not be greater than BUFFSIZE bytes)

using namespace std;
using Eigen::EigenSolver;
using Eigen::MatrixXd;
using Eigen::MatrixXcd;
using Eigen::VectorXd;
using Eigen::VectorXcd;

int main (int argc,char * argv[])
{
	/*FILE *fpb;	
	char * buf1,*buf2,*tmp_buffer;	//Used in reading from the file
	int d1;							//Max rank of y in the first polynomial
	int d2;							//Max rank of x in the first polynomial
	int d3;							//Max rank of y in the second polynomial
	int d4;							//Max rank of x in the second polynomial
	int dmax;						//Size of Sylvester
	int height;						//Size of the polynomials in Sylvester
	int i,j,k,el,x,y,sign;			//Used in iteration and insertingx,y,sign,*write;
	char hidden;					//Stores the hidden variable
	char answer;				
	double* write;
	srand(time(NULL));
	GenericMatrix<double> * pol1;
	GenericMatrix<double> * pol2;
	//Create Polynomials
	if (argc == 1) 					//Random Polynomials
	{ 
		int power1,power2;
		printf("Insert the two max polynomial powers (>1) \n");
		if(scanf("%d %d",&power1,&power2)!=2)
		{
			printf("Error \n");
			return 1;
		}
		d2=rand()%power1;
		d4=rand()%power2;
		if(d2==0)
			d2=1;
		if(d4==0)
			d4=1;
		d1=power1-d2;
		d3=power2-d4;
		pol1 = new GenericMatrix<double>(d1+1,d2+1);
		pol2 = new GenericMatrix<double>(d3+1,d4+1);
		for(i = 0 ; i<=d1 ; i++)
		{
			for(j = 0 ; j<=d2 ; j++)
			{
				sign = rand()%10;
				if(sign>4)
					sign=-1;
				else
					sign=1;
				pol1->SetMatrix(i,j,sign*rand()%10);
			}
		}
		pol1->SetMatrix(d1,d2,sign*rand()%10+1);			//x^maxpowerx * y^maxpowery must exist
		for(i = 0 ; i<=d3 ; i++)
		{
			for(j = 0 ; j<=d4 ; j++)
			{
				sign = rand()%10;
				if(sign>4)
					sign=1;
				else
					sign=-1;
				pol2->SetMatrix(i,j,sign*rand()%10);
			}
		}
		pol2->SetMatrix(d3,d4,sign*rand()%10+1);			//x^maxpowerx * y^maxpowery must exist
	}
	else	//Polynomials read from file
	{
		if (argc!=2) 
		{
      		printf("Correct syntax is: %s File\n", argv[0]);
      		return(1);
   		}
   		fpb = fopen (argv[1],"rb");
   		if (fpb==NULL) 
		{
      		printf("Cannot open file\n");
      		return 1;
   		}
		buf1 = (char *)malloc(sizeof(char)*BUFFSIZE);		//each polynomial is stored in a buffer
		buf2 = (char *)malloc(sizeof(char)*BUFFSIZE);
		fgets(buf1,BUFFSIZE,fpb);
		fgets(buf2,BUFFSIZE,fpb);
		findmaxpowers(buf1,&x,&y);							//max x,y powers of the polynomials
		d2=x;
		d1=y;
		pol1 = new GenericMatrix<double>(d1+1,d2+1);			//Allocated according to the powers
		findmaxpowers(buf2,&x,&y);
		d4=x;
		d3=y;
		pol2 = new GenericMatrix<double>(d3+1,d4+1);
		for(i=0;i<=d1;i++)									//Initialization
			for(j=0;j<=d2;j++)
				pol1->SetMatrix(i,j,0);
		for(i=0;i<=d3;i++)
			for(j=0;j<=d4;j++)
				pol2->SetMatrix(i,j,0);
		cout << buf1 << endl << endl;
		cout << buf2 << endl;
		tmp_buffer=buf1;
		createpolmatrix(buf1,tmp_buffer,pol1);				//1st polynomial in matrix form
  		free (tmp_buffer);
		tmp_buffer=buf2;
		createpolmatrix(buf2,tmp_buffer,pol2);				//2nd polynomial in matrix form
  		free (tmp_buffer);
		fclose (fpb);
	}	
	//Set Sylvester Matrix
	GenericMatrix<double *> * syl;
	if((d2>=d1&&d4>=d3)||(d2>=d1)&&(d2>=d3)||(d4>=d1)&&(d4>=d3))	//If the rank of x is greater or equal to the rank of y
	{
		height=max(d1,d3)+1;			//Height will be the max rank of y
		dmax=d2+d4;
		syl=new GenericMatrix<double *>(dmax,dmax,height);
		hidden='y';
	}
	else															//The rank of y is greater than the rank of x
	{
		height=max(d4,d2)+1;			//Height will be the max rank of x
		dmax=d1+d3;
		syl=new GenericMatrix<double *>(dmax,dmax,height);
		hidden='x';
	}
	cout << "The hidden variable is " << hidden << endl;
	write=new double[height];				//Array that will be written in Sylvester
	for(i=0; i<height; i++)				//Sylvester is all zeroes
		write[i]=0;
	for(i=0; i<dmax; i++) 
	{  
		for(j=0; j<dmax; j++) 
		{
			syl->SetMatrix(i,j,write);
		}
	}
	pol1->PrintMatrix();
	pol2->PrintMatrix();
	createSylvester(pol1,pol2,syl);
	cout << "There are " << d1 << " rows and " << d2 << " columns in the first matrix" << endl;
	cout << "There are " << d3 << " rows and " << d4 << " columns in the second matrix" << endl;
	height=syl->GetHigh();
	dmax=syl->GetRows();
	cout << "There are " << dmax << " rows and " << dmax << " columns in the Sylvester matrix " ; 
	cout << "and the rank of y is " << height-1 << endl;
	cout << "Printing Sylvester" << endl;
	syl->PrintMatrix();
	//Set Sylvester Polynomial
	
	SylvPol<GenericMatrix <double> > * sylp=new SylvPol<GenericMatrix <double> >(height);
	createPolynomials(syl,sylp);
	//sylp->PrintnMatrix(0);
	if(detChecks(sylp))
		cout << "Determinant check is true" << endl;
	else
	{
		cout << "Determinant check is false" << endl;
	}
	/*
	while(42)
	{
		cout << "Which coefficient of the sylvester polynomial do you want to be printed?" << endl;
		scanf("%d",&k);
		if(k<0||k>height-1)
			break;
		sylp->PrintnMatrix(k);
	}	
	cout << "Do you wish to perform a multiplication?(SylvesterMatrix * vector) - (y/n)\n" << endl;
	cin >> answer;
	if(answer == 'y')
	{
		GenericMatrix<double *> * multres = new GenericMatrix<double *>(syl->GetColumns(),1,syl->GetHigh());
		multiply(syl,multres);
		multres->PrintMatrix();
		delete multres;
	}*/
	
	/*  													/
	/			  Here starts the code of the second part of the project			/
	/												       */					
	/*int B=7;		//Variables used in the second part
	double bestKfactor,Kfactor;
	int c=0;
	int ts[5];
	ts[0]=rand()%10;
	ts[1]=rand()%10;
	ts[2]=rand()%10;
	ts[3]=rand()%10;
	ts[4]=0;
	//compute K
	if(computeK(sylp,B,height,&Kfactor))				//Calculate K and solve accordingly
	{
		StandardProblem(sylp,dmax,height-1,pol1,pol2,ts,hidden);
	}
	else
	{
		GeneralizedProblem(sylp,dmax,height-1,pol1,pol2,ts,hidden);	
	}
	bestKfactor=Kfactor;
	SylvPol<GenericMatrix <double> > * sylpZ=new SylvPol<GenericMatrix <double> >(height);	//Z polynomial
	ts[4]=1;
	for(int tries=0;tries<3;tries++)
	{	
		SylvPol<GenericMatrix <double> > * sylpZ=new SylvPol<GenericMatrix <double> >(height);	//Z polynomial
		ts[0]=rand()%10;
		ts[1]=rand()%10;
		ts[2]=rand()%10;
		ts[3]=rand()%10;
		CreateZpol(sylp,sylpZ,ts[0],ts[1],ts[2],ts[3]);	//Create Z polynomial with random ts y=(t1*z+t2)/(t3*z+t4)
		if(computeK(sylpZ,B,height,&Kfactor))		//Compute new K and solve accordingly if K(z)<K(y)
		{
			if(Kfactor<bestKfactor){		//K(z)<K(y)
				cout<<"yep"<<endl;
				StandardProblem(sylpZ,dmax,height-1,pol1,pol2,ts,hidden);
				break;
			}
				cout<<"nope"<<endl;		//K(z)>K(y)
		}
		else
		{
			if(Kfactor<bestKfactor){		//K(z)<K(y)
				GeneralizedProblem(sylpZ,dmax,height-1,pol1,pol2,ts,hidden);	
				cout<<"yep"<<endl;
				break;
			}
				cout<<"nope"<<endl;		//K(z)>K(y)
		}
	delete sylpZ;
	}
	delete pol1;
	delete pol2;
	//delete sylp;
	delete syl;
	delete [] write;
	return 42;*/
	double ** points=new double * [5];
	for(int i=0;i<5 ;i++) 
		points[i]=new double[2];
	points[0][0]=-1.0;
	points[0][1]=0.0;
	points[1][0]=4.0;
	points[1][1]=-1.0;
	points[2][0]=-1.0;
	points[2][1]=-5.0;
	points[3][0]=-4.0;
	points[3][1]=2.0;
	points[4][0]=4.0;
	points[4][1]=-4.0;
	int kpoints=5;
	int d=2;
	cout<<"Hey\n"<<endl;
	MatrixXd matrix(kpoints,kpoints+1);	
	createInterpolation(points,matrix,d,kpoints);
	int error;
	checkPst(&error,computeRank(matrix),kpoints);
	if(error==0)
		computeKernel(matrix);
	//cout<<"Hey\n"<<endl;

}

