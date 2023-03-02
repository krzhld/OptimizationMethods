#include <vector>
#include <cstdio>
#include <iostream>
#include <tuple>
#include <cmath>

using namespace std;

typedef vector<double> column_t;
typedef vector<column_t> matrix_t;
typedef tuple<matrix_t, column_t, column_t> problem_t;
typedef pair<column_t, double> solving_t;
typedef vector<int> comb_t;
typedef vector<comb_t> combinations_t;

/*Транспонирвоание и умножение на -1 матрицы для сведения канонической задачи к двойственной*/
matrix_t TransformationMatrix(matrix_t matrix)
{
	matrix_t transform;
	int n = size(matrix);
	int m = size(matrix[1]);
	transform.resize(m);
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			transform[j][i] = -matrix[i][j];
		}
	}
	return transform;
}

/*Получение двойственной задачи к канонической задаче*/
problem_t getDualLinearProblemForCanonicalProblem(problem_t problem) 
{
	matrix_t A; 
	column_t b, c;
	tie(A, b, c) = problem;
	matrix_t dualA; 
	column_t dualb, dualc;
	dualb.resize(size(c));
	dualc.resize(size(b));
	dualA = TransformationMatrix(A);
	for (int i = 0; i < size(c); i++)
	{
		dualb[i] = -c[i];
	}
	for (int j = 0; j < size(b); j++)
	{
		dualc[j] = -b[j];
	}
	problem_t dualProblem = make_tuple(dualA, dualb, dualc);
	return dualProblem;
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

/*Умножение транспонированного вектора на вектор такой же размерности*/
double MultipliedVectors(column_t v1, column_t v2)
{
	double result = 0;
	for (int i = 0; i < size(v1); i++)
	{
		result += v1[i] * v2[i];
	}
	return result;
}

/*Следующая комбинация в поиске сочетаний без повторений*/
bool NextSet(vector<int>& a, int n, int m)
{
	int k = m;
	for (int i = k - 1; i >= 0; --i)
		if (a[i] < n - k + i + 1)
		{
			a[i] += 1;
			for (int j = i + 1; j < k; ++j)
				a[j] = a[j - 1] + 1;
			return true;
		}
	return false;
}

/*Факториал числа*/
int Factorial(int n)
{
	int res = 1;
	for (int i = 1; i <= n; i++) {
		res = res * i;
	}
	return res;
}

/*Поиск всех сочетаний из n по m без повторений*/
combinations_t CombinationsWithoutRepetitions(int m, int n)
{
	comb_t temp, comb;
	int num = 0;
	combinations_t matrixCombinations;
	int numberCombinations = (int)(Factorial(n) / (Factorial(m) * Factorial(n)));
	temp.resize(n);
	comb.resize(m);
	matrixCombinations.resize(numberCombinations);
	for (int i = 0; i < n; i++)
		temp[i] = i + 1;
	for (int j = 0; j < m;j++)
	{
		comb[j] = temp[j];
	}
	matrixCombinations[num] = comb;
	while (NextSet(temp, n, m))
	{
		num++;
		for (int j = 0; j < m;j++)
		{
			comb[j] = temp[j];
		}
		matrixCombinations[num] = comb;
    } 
	return matrixCombinations;
}

/*Проверка совпадения числа и одного из значений в сочетании*/
bool IsNumberInCombination(int number, comb_t comb)
{
	for (int i = 0; i < size(comb); i++)
	{
		if (number == comb[i])
			return true;
	}
	return false;
}

// Получение матрицы без i-й строки и j-го столбца
void GetMatr(matrix_t& mas, matrix_t &p, int i, int j, int m)
{
	int ki, kj, di, dj;
	di = 0;
	for (ki = 0; ki < m - 1; ki++) 
	{ // проверка индекса строки
		if (ki == i) di = 1;
		dj = 0;
		for (kj = 0; kj < m - 1; kj++) 
		{ // проверка индекса столбца
			if (kj == j) dj = 1;
			p[ki][kj] = mas[ki + di][kj + dj];
		}
	}
}

// Рекурсивное вычисление определителя
double Determinant(matrix_t &mas, int m) 
{
	int i, j, k, n;
	double det;
	matrix_t p; 
	column_t tmp;
	p.resize(m);
	tmp.resize(m);
	for (i = 0; i < m; i++)
	{
		p[i] = tmp;
	}
	j = 0; det = 0;
	k = 1; //(-1) в степени i
	n = m - 1;
	if (m == 1) 
	{
		det = mas[0][0];
		return(det);
	}
	if (m == 2) 
	{
		det = mas[0][0] * mas[1][1] - (mas[1][0] * mas[0][1]);
		return(det);
	}
	if (m > 2) 
	{
		for (i = 0; i < m; i++) {
			GetMatr(mas, p, i, 0, m);
			det = det + k * mas[i][0] * Determinant(p, n);
			k = -k;
		}
	}
	return(det);
}

/*Решение задачи методом перебора крайних точек*/
solving_t IteratingThroughExtremePoints(problem_t problem)
{
	matrix_t A;
	column_t b, c;
	tie(A, b, c) = problem;
	int n = size(c),  m = size(b);
	if (n <= m) 
	{
		cout << "Incorred problem";
		column_t e1;
		e1[0] = 0;
		double e2 = 0;
		solving_t error = make_pair(e1, e2);
		return error;
	}
	combinations_t combinations = CombinationsWithoutRepetitions(m, n);
	matrix_t A1; 
	A1.resize(m);
	int numb=0;
	column_t X, solveSystem, minPoint;
	X.resize(n);
	solveSystem.resize(m);
	minPoint.resize(n);
	double min = 1e20;
	double valueInPoint; 
	for (auto& comb : combinations)
	{
		for (int i = 0; i < n; i++)
		{
			if (IsNumberInCombination(i, comb))
			{
				for (int j = 0; j < m; j++)
				{

					A1[numb][j] = A[i][j];
				}
				numb++;
			}
		}

		for (int i = 0; i < n;i++)
		{
			if (!IsNumberInCombination(i, comb))
			{
				X[i] = 0;
			}
		}
		if (Determinant(A1, m) != 0)
		{
			solveSystem = RotationMethod(A1, b);
			int numbInSolv = 0;
			for (int k = 0; k < n; k++)
			{
				if (IsNumberInCombination(k, comb))
				{
					X[k] = solveSystem[numbInSolv];
					numbInSolv++;
				}
			}
			valueInPoint = MultipliedVectors(c, X);
			if (valueInPoint < min)
			{
				minPoint = X;
				min = valueInPoint;
			}
		}
		
	}
	solving_t solving = make_pair(minPoint, min);
	return solving; 
}



int main()
{
	return 0;
}