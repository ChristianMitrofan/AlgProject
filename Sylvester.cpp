#include <stdio.h>
#include <stdlib.h>
#include "GenericMatrix.cpp"
#include "GenericMatrix.h"

void createSylvester(GenericMatrix<double> * pol1,GenericMatrix<double> * pol2,GenericMatrix<double *> * syl)
{
	int d1,d2,d3,d4,dmax,height,i,j,k;
	double * write;
	d1=pol1->GetRows()-1;
	d2=pol1->GetColumns()-1;
	d3=pol2->GetRows()-1;
	d4=pol2->GetColumns()-1;
	height=syl->GetHigh();
	write=new double[height];
	if((d2>=d1&&d4>=d3)||(d2>=d1)&&(d2>=d3)||(d4>=d1)&&(d4>=d3))	//If x has the greater rank then the y polynomials are the rows
	{	
		for(i=0; i<=d2; i++)	//For Columns
		{
			for(j=0; j<=d1; j++)	//For Rows
			{
				write[j]=pol1->GetMatrix(j,i);	//The element of the array is the coefficient while the position is the power			
			}
			for(k=0; k<d4; k++)				
			{
				syl->SetMatrix(k,k+i,write);	
			}
		}
		for(i=0; i<height; i++)
			write[i]=0;
		for(i=0; i<=d4; i++)	//For Columns
		{
			for(j=0; j<=d3; j++)	//For Rows
			{
				write[j]=pol2->GetMatrix(j,i);	//The element of the array is the coefficient while the position is the power			
			}
			for(k=0; k<d2; k++)				
			{
				syl->SetMatrix(k+d4,k+i,write); //The second polynomial will be written d4 rows lower
			}
		}
	}
	else															//If y has the greater rank then the x polynomials are the columns
	{
		for(i=0; i<=d1; i++)	//For Rows
		{
			for(j=0; j<=d2; j++)	//For Cols
			{
				write[j]=pol1->GetMatrix(i,j);	//The element of the array is the coefficient while the position is the power			
			}
			for(k=0; k<d3; k++)				
			{
				syl->SetMatrix(k,k+i,write);	//Den eimai sigouros oti einai k,k
			}
		}
		for(i=0; i<height; i++)
			write[i]=0;
		for(i=0; i<=d3; i++)	//For Columns
		{
			for(j=0; j<=d4; j++)	//For Rows
			{
				write[j]=pol2->GetMatrix(i,j);	//The element of the array is the coefficient while the position is the power			
			}
			for(k=0; k<d1; k++)				
			{
				syl->SetMatrix(k+d3,k+i,write); //The second polynomial will be written d3 rows lower
			}
		}
	}
	delete [] write;
}

void createPolynomials(GenericMatrix<double *> * syl,SylvPol<GenericMatrix <double> > * sylp)
{
	double el;
	int height=syl->GetHigh();
	int dmax=syl->GetColumns();
	for(int k=0; k<height; k++)				//The first coefficient of the sylvester polynomial is in position 0 of the array
	{										//and has in it all the coefficient of y that are of zero power
		GenericMatrix <double>  * newM = new GenericMatrix <double> (dmax,dmax);
		for(int i=0;i<dmax;i++)
			for(int j=0;j<dmax;j++)
			{
				el=syl->GetMatrix(i,j,k);
		 		newM->SetMatrix(i,j,el);
		 	}
		sylp->SetMatrixI(newM,k);
	}
}
