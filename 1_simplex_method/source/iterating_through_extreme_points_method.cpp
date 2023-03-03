#include "iterating_through_extreme_points_method.h"

using namespace std;

/*Нахождение ранга матрицы*/
int RankMatrix(matrix_t A, int N, int M)
{
	double temp;
	int rank = M;
	double sum=0;
	for (int i = 0; i < N - M; i++)
	{
		temp = A[i][i];
		for (int j = i + 1; j < M; j++)
		{
			for (int k = 0; k < N; k++)
			{
				A[k][j] = A[k][j] - A[i][j] / temp;
			}
		}
		
	}
	for (int i = 0; i < M; i++)
	{
		double sum = 0;
		for (int j = 0; j < N; j++)
		{
			sum += A[j][i];
		}
		if (sum == 0)
		{
			rank--;
		}
	}
	return rank;
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
bool NextSet(comb_t& a, int n, int m)
{
	int k = m;
	for (int i = k - 1; i >= 0; --i)
		if (a[i] < n - k + i + 1)
		{
			a[i] += 1;
			for (int j = i + 1; j < k; ++j)
			{
				a[j] = a[j - 1] + 1;
			}
			return true;
		}
	return false;
}

/*Факториал числа*/
long long int Factorial(int n)
{
	long long int res = 1;
	for (int i = 1; i <= n; i++) {
		res = res * i;
	}
	return res;
}

/*Поиск всех сочетаний из n по m без повторений*/
combinations_t CombinationsWithoutRepetitions(int m, int n)
{
	comb_t temp, comb;
	temp.resize(n); //промежуточный вектор
	comb.resize(m); //вектор целых чисел - один из вариантов сочетаний размерности m 


	int num = 0; //переменная - счетчик

	/*число всех сочетаний из n по m*/
	//int numberCombinations = (int)(Factorial(n) / (Factorial(m) * Factorial(n - m)));

	/*инициализация матрицы сочетаний*/
	combinations_t matrixCombinations;
	//matrixCombinations.resize(numberCombinations);

	matrixCombinations.resize(1);

	for (int i = 0; i < n; i++)
		temp[i] = i + 1;

	/*получение первого сочетания*/
	for (int j = 0; j < m;j++)
	{
		comb[j] = temp[j];
	}

	matrixCombinations[num] = comb;

	/*получение последующих сочетаний*/
	while (NextSet(temp, n, m))
	{
		num++;
		matrixCombinations.resize(num + 1);
		for (int j = 0; j < m;j++)
		{
			comb[j] = temp[j];
		}
		matrixCombinations[num] = comb;
	}

	return matrixCombinations;
}

/*Проверка совпадения числа и одного из значений в сочетании*/
bool IsNumberInCombination(int number, comb_t comb, int combSize)
{
	for (int i = 0; i < combSize; i++)
	{
		if (number == comb[i])
			return true;
	}
	return false;
}

/*Определитель матрицы*/
double Determinant(matrix_t matrix)
{
	int index; 
	double det = 1, total = 1;
	double num1, num2;
	int n = matrix.size();
	column_t temp; // временный массив для хранения столбца
	temp.resize(n + 1);

	/* цикл для обхода диагональных элементов*/ 
	for (int i = 0; i < n; i++)
	{
		index = i; // инициализация индекса

		// поиск индекса, который имеет ненулевое значение
		while (index < n && matrix[index][i] == 0)
		{
			index++;
		}
		if (index == n) // если существует ненулевой элемент
		{
			continue; // определитель матрицы равен нулю
		}
		if (index != i) 
		{
			/* цикл для замены строки диагонального элемента и элемента индексной строки*/
			for (int j = 0; j < n; j++)
			{
				swap(matrix[index][j], matrix[i][j]);
			}
			// знак определителя меняется, когда мы сдвигаем строки
			det = det * pow(-1, index - i);
		}

		// сохранение значений элементов диагональной строки
		for (int j = 0; j < n; j++)
		{
			temp[j] = matrix[i][j];
		}
		
		// обход каждой строки под диагональным элементом
		for (int j = i + 1; j < n; j++)
		{
			num1 = temp[i]; // значение диагонального элемента
			num2 = matrix[j][i]; // значение элемента следующей строки

			/* обход каждого элемента строки и умножение на каждую строку */
			for (int k = 0; k < n; k++)
			{
				matrix[j][k] = (num1 * matrix[j][k]) - (num2 * temp[k]);
			}
			total = total * num1; // Det(kA)=kDet(A);
		}
	}

	// умножение диагональных элементов для получения определителя
	for (int i = 0; i < n; i++)
	{
		det = det * matrix[i][i];
	}

	// Det(kA)/k=Det(A);
	return (det / total);
}

/*Проверка неотрицательности компонент вектора*/
bool NonNegativityOfVector(column_t v)
{
	for (int i = 0; i < size(v); i++)
	{
		if (v[i] < 0)
			return false;
	}
	return true;
}

/*Решение задачи методом перебора крайних точек*/
solving_t IteratingThroughExtremePoints(canon_problem_t &problem)
{
	matrix_t A;
	column_t b, c;
	tie(A, b, c) = problem;
	int n = size(c), m = size(b);

	/*проверка корректности задачи*/
	if (n <= m)
	{
		cout << "Incorred problem";
		exit(-1);
	}

	/*получение всех комбинаций из n (количество переменных) по m (количество ограничений)*/
	combinations_t combinations = CombinationsWithoutRepetitions(m, n);

	column_t maxPoint; // вектор, в котором достигается оптимальное решение
	maxPoint.resize(n);
	double max = -1e20; // максимальное значение функции цели
	comb_t optComb; // базисная комбинация, при которой достигается оптимальное решение
	optComb.resize(m);


	/*перебор всех возможных комбинаций сочетаний*/
	for (auto& comb : combinations)
	{

		double valueInPoint; // значение функции цели при данной комбинации
		column_t X; // опорный вектор при данной комбинации 
		X.resize(n);
		column_t columnA1; //матрица для решения системы уравнений, из которой отброшены столбцы , номера которых не входят в данную комбинацию 

		/*размерность матрицы A1 m на m*/
		columnA1.resize(m); 
		matrix_t A1;
		A1.resize(m); 
		for (int k = 0; k < m;k++)
		{
			A1[k] = columnA1;
		}

		column_t solveSystem; //решение СЛАУ с матрицей A1 и вектором b, имеет размерность 
		solveSystem.resize(m);


		int numb = 0; // счетсчик для задания компонент матрицы A1

		/*задание компонент матрицы A1, берутся те столбцы, номера которых есть в данной комбинации*/
		for (int i = 0; i < n; i++)
		{
			if (IsNumberInCombination(i+1, comb, m))
			{
				for (int j = 0; j < m; j++)
				{

					A1[j][numb] = A[i][j];
				}
				numb++;
			}
		}

		/*те компоненты опорного вектора, номера которых не входят в данную комбинацию принимают нулевые значения*/
		for (int i = 0; i < n;i++)
		{
			if (!IsNumberInCombination(i+1, comb, m))
			{
				X[i] = 0;
			}
		}

		if (Determinant(A1) != 0) // проверка невырожденности матрицы
		{
			/*решение СЛАУ методом вращений*/
			solveSystem = RotationMethod(A1, b);

			if (NonNegativityOfVector(solveSystem)) // проверка на неотрицательность компонент решения СЛАУ (иначе вектор не будет допустимым)
			{
				int numbInSolv = 0; //счетчик 

			/*те компоненты опорного вектора, номера которых входят в данную комбинацию принимают соответствующие значения решения СЛАУ*/
				for (int k = 0; k < n; k++)
				{
					if (IsNumberInCombination(k + 1, comb, m))
					{
						X[k] = solveSystem[numbInSolv];
						numbInSolv++;
					}
				}

				/*вычисление значения целевой функции при данном опорном векторе*/
				valueInPoint = MultipliedVectors(c, X);

				/*если полученное значение функции цели меньше минимума, то оно становится минимумом, 
				вектору, в котором достигается минимум присваивается значение X*/
				if (valueInPoint > max)
				{
					maxPoint = X;
					max = valueInPoint;
					optComb = comb;
				}
					
				
			}
			

		}
		
	}
	/*создание кортежа из вектора, в котором достигается оптимальное решение и оптимального решения*/
	solving_t solving = make_tuple(maxPoint, max, optComb);
	return solving;
}
