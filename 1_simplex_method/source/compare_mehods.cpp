#pragma once
#include <vector>
#include <cstdio>
#include <iostream>
#include <tuple>
#include <cmath>

using namespace std;

#define N 15

#pragma warning(disable:4996)

// вектор вещественных чисел
typedef vector<double> column_t;

// матрица вещественных чисел
typedef vector<column_t> matrix_t;

/*считывание матриц*/
matrix_t ArrayRead(FILE* file)
{

	column_t str;
	str.resize(N);
	matrix_t arr;
	arr.resize(N);
	for (int i = 0; i < N; i++)
	{
		arr[i] = str;
	}

	for (int j = 0;j < N;j++)
	{
		for (int i = 0;i < N;i++)
		{
			fscanf(file, "%lf", &arr[j][i]);
			
		}
		
	}
	return arr;
}

/*считывание столбцов*/
column_t ColumnRead(FILE* file)
{
	column_t arr;
	arr.resize(N);

	for (int i = 0;i < N;i++)
	{
		fscanf(file, "%lf", &arr[i]);
	}
	return arr;
}

/*Метод вращений решения СЛАУ*/
column_t RotationMethod(matrix_t A, column_t b)
{
	long double c, s;
	int sz = size(b);
	int k = 0;
	column_t X;
	X.resize(sz);
	for (int i = 0;i < sz - 1;i++)
	{
		for (int j = i + 1;j < sz;j++)
		{
			c = A[i][k] / sqrt(A[i][k] * A[i][k] + A[j][k] * A[j][k]);
			s = A[j][k] / sqrt(A[i][k] * A[i][k] + A[j][k] * A[j][k]);
			for (int m = k;m < sz;m++)
			{
				long double a1 = c * A[i][m] + s * A[j][m];
				long double a2 = c * A[j][m] - s * A[i][m];
				A[i][m] = a1;
				A[j][m] = a2;
			}

			long double b1 = c * b[i] + s * b[j];
			long double b2 = c * b[j] - s * b[i];
			b[i] = b1;
			b[j] = b2;
		}
		k++;
	}
	for (int i = sz - 1;i >= 0;i--)
	{
		long double tmp = b[i];
		for (int j = sz - 1;j > i;j--)
		{
			tmp -= A[i][j] * X[j];
		}
		X[i] = tmp / A[i][i];
	}
	return X;
}

column_t GaussMethod(matrix_t A, column_t b)
{
    int n = size(b);
    double temp;
    column_t X;
    X.resize(n);
    int k = 0, index;
	//прямой ход
	for (int i = 0; i < n-1;i++)
	{
		temp = A[i][i];
		for (int j = i + 1; i < n; i++)
		{

			for (int k = 0; k < n; k++)
			{
				A[k][j] = A[k][j] - A[i][j] / temp;
			}
			b[j] = b[j] - b[i] / temp;
		}
	
	}
   
    // обратная подстановка
    for (k = n - 1; k >= 0; k--)
    {
        X[k] = b[k];
        for (int i = 0; i < k; i++)
            b[i] = b[i] - A[i][k] * X[k];
    }
    return X;
}

int main()
{
	FILE* rotation = fopen("rotationresult.txt", "w");
	FILE* gauss = fopen("gaussresult.txt", "w");

	FILE* af = fopen("A.txt", "r");
	FILE* bf = fopen("b.txt", "r");
	for (int i = 0;i < 10;i++)
	{
		matrix_t A;
		column_t b;
		column_t Xrotation;
		column_t Xgauss;
		A = ArrayRead(af);
		b = ColumnRead(bf);
		Xrotation = RotationMethod(A, b);
		Xgauss = GaussMethod(A, b);
		for (int i = 0; i < N; i++)
		{
			fprintf(rotation, "%.16lf\n ", Xrotation[i]);
			fprintf(gauss, "%.16lf\n ", Xgauss[i]);
		}

	}
	fclose(af);
	fclose(bf);


	fclose(rotation);
	fclose(gauss);
	return 0;
}