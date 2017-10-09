#ifndef POLYNOMIALS_H_
#define POLYNOMIALS_H_
#include "GenericMatrix.h"

int findmaxpowers(char * input,int * x,int * y);

int findnext(char * buffer,char const * n);
			
int findpower(char * buffer,int * pass);

double findcoeff(char * buffer,int * pass);

void createpolmatrix(char * buffer,char * tmp_buffer,GenericMatrix<double> * pol);

void RandomizePolynomial(GenericMatrix<double> * pol,int r,int c);

#endif
