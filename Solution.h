#ifndef SOLUTION_H_
#define SOLUTION_H_

#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include "Companion.h"
#include "Lmatrices.h"
#include "Calculations.h"
#include "GenericMatrix.h"


bool StandardProblem(SylvPol <GenericMatrix<double> > * pol,int size,int d,GenericMatrix<double> * pol1,GenericMatrix<double> * pol2,int * tz,char hvar);

bool GeneralizedProblem(SylvPol <GenericMatrix<double> > * pol,int size,int d,GenericMatrix<double> * pol1,GenericMatrix<double> * pol2,int * tz,char hvar);

#endif

