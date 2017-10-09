#include "Calculations.h"
#include "GenericMatrix.cpp"
#include "GenericMatrix.h"
#include <iostream>
#include <fstream>


using Eigen::MatrixXd;
using Eigen::MatrixXcd;
using Eigen::VectorXd;
using Eigen::VectorXcd;
using Eigen::EigenSolver;
using Eigen::JacobiSVD;

double fRand(double fMin, double fMax)	//Random double numbers
{
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}
int detChecks(SylvPol <GenericMatrix<double> > * pol){		//Bonus no1 checking with Random y if dets are 0
	int flag = 1;						//return 1 if not 0 if yes
	pol->PrintMatrix();
	int dmax = pol->GetMatrixI(0)->GetRows();
 	MatrixXd Md(dmax,dmax);
	MatrixXd ress(dmax,dmax);
	MatrixXd sums=MatrixXd::Zero(dmax,dmax);
 	int Size = pol->GetSizeI();
	double check=fRand(-5,5);
	for(int times=0;times<4;times++){
		for(int pols=0;pols<Size;pols++){	
			for(int i=0;i<dmax;i++)
 			{
 	 			for(int j=0;j<dmax;j++)
 				{
  				 Md(i,j) = pol->GetMatrixI(pols)->GetMatrix(i,j);
  				}
 			}
			sums+=ress;
		}
		if (sums.determinant()!=0){
			flag=1;
			break;
		}
		check=fRand(-5,5);
	}	
	if (flag == 0)
		return 0;
	else 
		return 1;
}	


int computeK(SylvPol <GenericMatrix <double > > * sylp,int B,int height,double * K){	//returns 0 if it is not invertible or 1 if it is	
	int dmax = sylp->GetMatrixI(0)->GetRows();					//k is stored in double * K
 	MatrixXd Md(dmax,dmax);
 	for(int i=0;i<dmax;i++)								//Md Matrix
 	{
 	 	for(int j=0;j<dmax;j++)
 		{
  		 Md(i,j) = sylp->GetMatrixI(height-1)->GetMatrix(i,j);
  		}
 	}
	cout <<"Printing Md Matrix : "<< endl << Md << endl;
	int s =Md.rows();
	JacobiSVD<MatrixXd> svd(Md);							//Singular Values Calculation
	VectorXd sings = svd.singularValues();
	if (sings(0)==0){
		cout << "k = inf  Bound: ill-conditioned Md, generalized eigenproblem" << endl;
		return 0;
	}
	cout<<sings(0)<<endl;
	cout<<sings(s-1)<<endl;
	double result = (sings(0)/sings(s-1));						//K calculation(division with 0 is excluded)
	*K=result;	
	if(result<pow(10,B)){								//print and return according to K
		cout << "k = " << result << " Bound: non-singular  Md, standard eigenproblem" << endl;
		return 1;
	}else{
		cout << "k = " << result << " Bound: ill-conditioned Md, generalized eigenproblem" << endl;
		return 0;
	}
}

bool GeneralizedEigenProblem(MatrixXd & A,MatrixXd & B,MatrixXd & v,MatrixXd & lambda){ //Eigenvectors and Eigenvalues if Matrix is not invertible
	int N = A.cols();
	if((B.cols()!=N) || (A.rows()!=N) || (B.rows()!=N)	)
		return -1;
	v.resize(N,N);
	lambda.resize(N,3);
	int LDA = A.outerStride();
	int LDB = B.outerStride();
	int LDV = v.outerStride();
	int INFO = 0;

	double * alphar = lambda.col(0).data();
	double * alphai = lambda.col(1).data();
	double * beta = lambda.col(2).data();

	INFO = LAPACKE_dggev(LAPACK_COL_MAJOR,'N','V',N,A.data(),LDA,B.data(),LDB,alphar,alphai,beta,0,LDV,v.data(),LDV);
	return INFO==0;
}


VectorXcd StandardEigenProblemVals(MatrixXd & A){		//Eigenvalues if Matrix is invertible

	EigenSolver<MatrixXd> es(A); 
	VectorXcd evals = es.eigenvalues();
	return evals;
}

MatrixXcd StandardEigenProblemVecs(MatrixXd & A){		//Eigenvectors if Matrix is invertible

	EigenSolver<MatrixXd> es(A); 
	MatrixXcd evecs = es.eigenvectors();
	return evecs;
}

void PolynomialRoots(GenericMatrix<double> *pol1, GenericMatrix<double> *pol2, double * vx, double * vy, double * vm, int c, int * tz, char hidden)
{
	double x,y;
	bool flag=false;
	double sum=0;
	int pow1y=pol1->GetRows();
	int pow1x=pol1->GetColumns();
	int pow2y=pol2->GetRows();
	int pow2x=pol2->GetColumns();
   /* for(int i=0 ; i<c ; i++){
        vx[i]=(tz[1]-tz[3]*vx[i])/(tz[2]*vx[i]-tz[0]);
        vy[i]=(tz[1]-tz[3]*vy[i])/(tz[2]*vy[i]-tz[0]);
    }*/

    ofstream myfile ("roots.txt");

    for (int i=0; i<c; i++)
	{
		if(vm[i]==1)
		{
			sum=0;
			flag = false;
			for(int rows=0;rows<pow1y;rows++)
				for(int cols=0;cols<pow1x;cols++)
					sum+=pol1->GetMatrix(rows,cols)*pow(vx[i],cols)*pow(vy[i],rows);
	
			if(sum<pow(10,-5))
				flag = true ;

			sum = 0;
			for(int rows=0;rows<pow2y;rows++)
				for(int cols=0;cols<pow2x;cols++)
					sum+=pol2->GetMatrix(rows,cols)*pow(vx[i],cols)*pow(vy[i],rows);

            //if(sum>pow(10,-5))
                //flag = false ;

			if(hidden=='y')
			{
                if((sum<pow(10,-5)) && (flag==true)){
                    cout << "y = " << vy[i] << " x = " << vx[i] <<  endl;
                    if (myfile.is_open())
                    {
                        myfile << vx[i] << "    " << vy[i] <<endl;
                    }
                       else cout << "Unable to open file";
                }
			}
			else
			{
                if((sum<pow(10,-5)) && (flag==true)){
					cout << "y = " << vx[i] << " x = " << vy[i] <<  endl;
                    if (myfile.is_open())
                    {
                        myfile << vy[i] << "    " << vx[i] <<endl;
                    }
                       else cout << "Unable to open file";
                }
			}	
		}
    }
   myfile.close();
}

