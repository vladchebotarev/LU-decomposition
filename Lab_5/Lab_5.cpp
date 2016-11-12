#include "stdafx.h"
#include <iostream>
#include <Windows.h>
#include<stdio.h>
#include<math.h>
#include <float.h>      
#include <limits>  

using namespace std;
const double eps = 1e-7;

void MatrixPrint(double **matrix, int indeksy[4])
{

	for (int i = 0; i<4; i++)
	{
		cout << "| ";
		for (int j = 0; j<4; j++)
		{
			cout.width(15);
			cout << matrix[indeksy[i]][j];
		}
		cout << " |" << endl;
	}
}

void VectorPrint(double *vector, int indeksy[4])
{
	for (int i = 0; i<4; i++)
	{

		cout << "| ";
		cout.width(8);
		cout << vector[indeksy[i]];
		cout << "| ";
		cout << endl;
	}
}

void DataInsert(double **A, double *b)
{
	A[0][0] = 1.0;   A[0][1] = -20.0;  A[0][2] = 30.0;   A[0][3] = -4.0;
	A[1][0] = 2.0;   A[1][1] = -40.0;  A[1][2] = -6.0;   A[1][3] = 50.0;
	//A[1][0] = -16.0;   A[1][1] = 15.0;  A[1][2] = -140.0;   A[1][3] = 13.0;
	A[2][0] = 9.0;   A[2][1] = -180.0; A[2][2] = 11.0;   A[2][3] = -12.0;
	A[3][0] = -16.0; A[3][1] = 15.0;   A[3][2] = -140.0; A[3][3] = 13.0;
	//A[3][0] = 2.0; A[3][1] = -40.0;   A[3][2] = -6.0; A[3][3] = 50.0;

	b[0] = 35.0; b[1] = 104.0; b[2] = -366.0; b[3] = -354.0;
}

bool PartSelection(double **A, int indeks[4], int x)
{
	double max = 0;
	int indeks_max = x, temp;
	bool search = false;

	for (int i = x; i < 4; ++i)
	{
		if (fabs(A[indeks[i]][x]) > max)
		{
			max = fabs(A[indeks[i]][x]);
			indeks_max = i;
			search = true;
		}
	}

	if (search)
	{
		temp = indeks[indeks_max];
		indeks[indeks_max] = indeks[x];
		indeks[x] = temp;
	}

	return search;

}
void DecompositionLU()
{
	double **A, **L;
	double *b, *y, *x, wspl, temp;
	int indeksy[4];
	int n = 4, i, j, k;


	A = new double *[n];
	L = new double *[n];
	b = new double[n];
	y = new double[n];
	x = new double[n];

	for (i = 0; i < n; ++i)
		indeksy[i] = i;

	for (i = 0; i < n; ++i)
		A[i] = new double[n];

	for (i = 0; i < n; ++i)
		L[i] = new double[n];

	for (i = 0; i < n; ++i)
	{
		for (j = 0; j < n; ++j)
			L[i][j] = 0;
	}

	DataInsert(A, b);
	cout << "Macierz A:" << endl;
	MatrixPrint(A, indeksy);
	cout << endl;
	cout << " Wektor b:" << endl;
	VectorPrint(b, indeksy);
	cout << endl;

	for (i = 0; i < 3; i++)
	{
		for (j = i + 1; j < 4; j++)
		{
			if (A[indeksy[i]][indeksy[i]] == 0.0)
			{
				if (!PartSelection(A, indeksy, i))
				{


					break;
				}
			}

			L[indeksy[j]][i] = A[indeksy[j]][i] / A[indeksy[i]][i];
			/*MatrixPrint( L, indeksy );
			cout << L[j][i] << endl;*/
			wspl = A[indeksy[j]][i] / A[indeksy[i]][i];
			for (k = i; k < 4; k++)
			{
				A[indeksy[j]][k] = A[indeksy[j]][k] - A[indeksy[i]][k] * wspl; //U


			}
			cout << " Macierz U1" << endl;
			MatrixPrint(A, indeksy);
			cout << endl;
		}
	}

	for (i = 0; i < 4; ++i)
	{
		L[indeksy[i]][i] = 1;
	}

	cout << " matrixierz U:" << endl;
	MatrixPrint(A, indeksy);
	cout << endl;
	cout << " matrixierz L:" << endl;
	MatrixPrint(L, indeksy);
	cout << endl;
	//cout << " Wektor b:" << endl;
	//VectorPrint( b, indeksy );
	//cout << endl;

	for (i = 0; i < 4; i++)
	{
		temp = 0;
		for (j = 0; j < i; j++)
		{
			temp += L[indeksy[i]][j] * y[indeksy[j]];
		}
		y[indeksy[i]] = (b[indeksy[i]] - temp);
		//cout << "y: " << y[indeksy[i]] << endl;

	}

	cout << " Wektor b: " << endl;
	VectorPrint(b, indeksy);
	cout << endl;

	cout << " Wektor y: " << endl;
	VectorPrint(y, indeksy);
	cout << endl;

	for (i = n - 1; i >= 0; i--)
	{
		temp = 0.0;
		for (j = i + 1; j < 4; j++)
		{
			temp += A[indeksy[i]][j] * x[indeksy[j]];
		}

		if (A[indeksy[i]][i] == 0)
			A[indeksy[i]][i] = DBL_MIN;

		x[indeksy[i]] = (1 / A[indeksy[i]][i]) *(y[indeksy[i]] - temp);

	}

	cout << " Wektor x: " << endl;
	VectorPrint(x, indeksy);

	for (i = 0; i < n; i++)
	{
		delete[] A[i];
		delete[] L[i];
	}

	delete[] A;
	delete[] L;
	delete[] b;
	delete[] x;
	delete[] y;

}

int main()
{
	DecompositionLU();

	system("pause");
	return 0;
}