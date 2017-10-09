#include <string>
#include "GenericMatrix.h"

using namespace std;

template <typename T>
GenericMatrix<T>::GenericMatrix(int r,int c)
{
	rows=r;
	cols=c;
	matrix=new T*[rows]; 
    for (int i=0; i<rows; ++i)
		matrix[i]=new T[cols]; 
    for(int i=0;i<r;i++)
        for(int j=0;j<c;j++)
            matrix[i][j]=0;
}

template <typename T>
GenericMatrix<T>::~GenericMatrix()
{
	for(int i=0; i<rows; i++)
		delete [] matrix[i];
	delete [] matrix;	
}

template <typename T>
int GenericMatrix<T>::GetRows()
{
	return rows;
}

template <typename T>
int GenericMatrix<T>::GetColumns()
{
	return cols;
}


template <typename T>
void GenericMatrix<T>::PrintMatrix() 
{ 
    cout<<"Current Matrix"<<endl;
    for(int i=0; i<rows; i++) 
    {  
		for(int j=0; j<cols; j++) 
		{ 
			cout << matrix[i][j] << " " ; 
		}
		cout<<endl;
    } 
 }
 
template <typename T>
void GenericMatrix<T>::SetMatrix(int r,int c,T element)
{
	matrix[r][c]=element;
}

template <typename T>
T GenericMatrix<T>::GetMatrix(int r,int c)
{
	return matrix[r][c];
}

template <typename T>
GenericMatrix<T *>::GenericMatrix(int r,int c,int h)
{
	{
		rows=r;
		cols=c;
		height=h;
		matrix=new T**[rows]; 
		for (int i=0; i<rows; i++)
		{
			matrix[i]=new T*[cols];
			for(int j=0; j<cols; j++)
				matrix[i][j]=new T[h];
		}
	}
}

template <typename T>
GenericMatrix<T *>::~GenericMatrix()
{
	for(int i=0; i<rows; i++)
	{
		for(int j=0; j<cols; j++)
			delete matrix[i][j];
		delete [] matrix[i];
	}
	delete [] matrix;	
}

template <typename T>
void GenericMatrix<T *>::PrintMatrix() 
{
	cout<<"Current Matrix"<<endl;
    for(int i=0; i<rows; i++) 
    {  
		for(int j=0; j<cols; j++) 
		{ 
			for(int k=0; k<height; k++)
			{
				if(matrix[i][j][k]==0)
					cout<< 0 ;
				else
					cout << matrix[i][j][k] << "y^" << k ;
				if(k!=height-1)
					cout << "+" ;
			}
			cout << " " ;
			//cout << matrix[i][j][0] << "y^" << 0 << "+" << matrix [i][j][1] << "y^" << 1 << "+"<< matrix [i][j][2] << "y^" << 2 << " " ; 
		}
		cout<<endl;
    } 
}

template <typename T>
void GenericMatrix<T *>::SetMatrix(int r,int c,T element[])
{
	for(int i=0; i<height; i++)
		matrix[r][c][i]=element[i];
}

template <typename T>
T* GenericMatrix<T *>::GetMatrix(int r,int c)
{
	return matrix[r][c];
}

template <typename T>
T GenericMatrix<T *>::GetMatrix(int r,int c,int h)
{
	return matrix[r][c][h];
}

template <typename T>
int GenericMatrix<T *>::GetRows()
{
	return rows;
}

template <typename T>
int GenericMatrix<T*>::GetColumns()
{
	return cols;
}

template <typename T>
int GenericMatrix<T*>::GetHigh()
{
	return height;
}
template <typename T>
SylvPol<T>::SylvPol(int p)
{
	size=p;
	arr=new T * [size] ; 
}

template <typename T>
SylvPol<T>::~SylvPol()
{
	for(int i=0; i<size; i++)
		arr[i]->~T ();
	delete [] arr;	
}

template <typename T>
void SylvPol<T>::PrintnMatrix(int n)
{
	arr[n]->PrintMatrix();  
}

template <typename T>
void SylvPol<T>::PrintMatrix()
{
	for(int i=0;i<size;i++)
		PrintnMatrix(i);
}

template <typename T>
void SylvPol<T>::SetMatrixI(T * newM,int power)
{
	arr[power]=newM;  
}

template <typename T>
T * SylvPol<T>::GetMatrixI(int power)
{
	return arr[power];  
}

template <typename T>
int SylvPol<T>::GetSizeI()
{
	return size;  
}
