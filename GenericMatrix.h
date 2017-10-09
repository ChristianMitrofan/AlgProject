#ifndef GENERIC_MATRIX_H
#define GENERIC_MATRIX_H

#include <iostream>
#include <cstdlib>

using namespace std;

template <typename T>
class GenericMatrix
{
	private:
		int rows;
		int cols;
		T **matrix;
	
	public:
		GenericMatrix(int r,int c);
		~GenericMatrix();
		void PrintMatrix();
		void SetMatrix(int r,int c,T element);
		T GetMatrix(int r,int c);
		int GetRows();
		int GetColumns();
};

template <typename T>
class GenericMatrix<T *>
{
	private:
		int rows;
		int cols;
		int height;
		T ***matrix;
	
	public:
		GenericMatrix(int r,int c,int h);
		~GenericMatrix();
		void PrintMatrix();
		void SetMatrix(int r,int c,T element[]);
		T* GetMatrix(int r,int c);
		T GetMatrix(int r,int c,int h);
		int GetRows();
		int GetColumns();
		int GetHigh();
};

template <typename T>
class SylvPol
{
	private:
		int size;
		T ** arr;
	
	public:
		SylvPol(int p);
		~SylvPol();
		void PrintnMatrix(int n);
		void PrintMatrix();
		void SetMatrixI(T * N,int power);
		T * GetMatrixI(int n);
		int GetSizeI();
};
#endif
