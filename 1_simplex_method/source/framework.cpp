#include "framework.h"
using namespace std;

#pragma warning(disable:4996)


bool NextSet(comb_t& a, int n, int m) {
	int k = m;
	for (int i = k - 1; i >= 0; --i)
		if (a[i] < n - k + i + 1) {
			a[i] += 1;
			for (int j = i + 1; j < k; ++j)
				a[j] = a[j - 1] + 1;
			return true;
		}
	return false;
}

combinations_t CombinationsWithoutRepetitions(int m, int n) {
	comb_t temp, comb;
	temp.resize(n); //промежуточный вектор
	comb.resize(m); //вектор целых чисел - один из вариантов сочетаний размерности m 

	int num = 0; //переменная - счетчик

	/*инициализация матрицы сочетаний*/
	combinations_t matrixCombinations;

	matrixCombinations.resize(1);

	for (int i = 0; i < n; i++)
		temp[i] = i + 1;

	/*получение первого сочетания*/
	for (int j = 0; j < m; j++)
		comb[j] = temp[j];

	matrixCombinations[num] = comb;

	/*получение последующих сочетаний*/
	while (NextSet(temp, n, m)) {
		num++;
		matrixCombinations.resize(num + 1);
		for (int j = 0; j < m; j++)
			comb[j] = temp[j];

		matrixCombinations[num] = comb;
	}

	return matrixCombinations;
}

bool IsNumberInCombination(int number, comb_t comb) {
	for (int i = 0; i < size(comb); i++) {
		if (number == comb[i])
			return true;
	}
	return false;
}

bool NonNegativityOfVector(column_t v) {
	for (int i = 0; i < size(v); i++) {
		if (v[i] < 0)
			return false;
	}

	return true;
}

bool NonPositivityOfVector(column_t v) {
	for (int i = 0; i < size(v); i++) {
		if (v[i] > 0)
			return false;
	}

	return true;
}

column_t MultipliedMatrixAndColumn(matrix_t M, column_t c) {
	int n = size(c);

	column_t X(n); ///результат перемножения

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++)
			X[i] += M[i][j] * c[j];
	}

	return X;
}

int RankMatrix(matrix_t A) {
	int M = A.size();
	int N = A[0].size();
	double temp;
	int rank = M;
	double sum = 0;
	for (int i = 0; i < N - M; i++) {
		temp = A[i][i];
		for (int j = i + 1; j < M; j++) {
			for (int k = 0; k < N; k++)
				A[j][k] = A[j][k] - A[j][i] / temp;
		}
	}
	for (int i = 0; i < M; i++) {
		double sum = 0;
		for (int j = 0; j < N; j++)
			sum += A[i][j];
		if (sum == 0)
			rank--;
	}
	return rank;
}

matrix_t InverseMatrix(matrix_t A) {
	int n = size(A);

	column_t Ecolumn; // строка единичной матрицы
	Ecolumn.resize(n);

	matrix_t E; // инициализация единичной матрицы
	E.resize(n);
	for (int i = 0; i < n; i++)
		E[i] = Ecolumn;

	double temp;
	for (int i = 0; i < n; i++)
		E[i][i] = 1;

	for (int k = 0; k < n; k++) {
		temp = A[k][k];

		for (int j = 0; j < n; j++) {
			A[k][j] /= temp;
			E[k][j] /= temp;
		}

		for (int i = k + 1; i < n; i++) {
			temp = A[i][k];

			for (int j = 0; j < n; j++) {
				A[i][j] -= A[k][j] * temp;
				E[i][j] -= E[k][j] * temp;
			}
		}
	}

	for (int k = n - 1; k > 0; k--) {
		for (int i = k - 1; i >= 0; i--) {
			temp = A[i][k];

			for (int j = 0; j < n; j++) {
				A[i][j] -= A[k][j] * temp;
				E[i][j] -= E[k][j] * temp;
			}
		}
	}

	// преобразованная единичная матрица является обратной матрицей к матрице А 
	return E;
}

column_t RotationMethod(matrix_t A, column_t b) {
	long double c, s;
	int sz = size(b);
	int k = 0;
	column_t X;
	X.resize(sz);
	for (int i = 0; i < sz - 1; i++) {
		for (int j = i + 1; j < sz; j++) {
			c = A[i][k] / sqrt(A[i][k] * A[i][k] + A[j][k] * A[j][k]);
			s = A[j][k] / sqrt(A[i][k] * A[i][k] + A[j][k] * A[j][k]);
			for (int m = k; m < sz; m++) {
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
	for (int i = sz - 1; i >= 0; i--)
	{
		long double tmp = b[i];
		for (int j = sz - 1; j > i; j--)
		{
			tmp -= A[i][j] * X[j];
		}
		X[i] = tmp / A[i][i];
	}
	return X;
}

column_t GaussMethod(matrix_t A, column_t b) {
	int n = size(b);
	double temp;
	column_t X;
	X.resize(n);
	int k = 0;
	// прямой ход
	for (int i = 0; i < n - 1; i++) {
		temp = A[i][i];
		for (int j = i + 1; i < n; i++) {

			for (int k = 0; k < n; k++)
				A[k][j] = A[k][j] - A[i][j] / temp;

			b[j] = b[j] - b[i] / temp;
		}
	}

	// обратная подстановка
	for (k = n - 1; k >= 0; k--) {
		X[k] = b[k];
		for (int i = 0; i < k; i++)
			b[i] = b[i] - A[i][k] * X[k];
	}
	return X;
}

double MultipliedVectors(column_t v1, column_t v2) {
	double result = 0;
	for (int i = 0; i < size(v1); i++)
		result += v1[i] * v2[i];

	return result;
}

double Determinant(matrix_t matrix) {
	int index;
	double det = 1, total = 1;
	double num1, num2;
	int n = matrix.size();
	column_t temp; // временный массив для хранения столбца
	temp.resize(n + 1);

	// цикл для обхода диагональных элементов
	for (int i = 0; i < n; i++) {
		index = i; // инициализация индекса

		// поиск индекса, который имеет ненулевое значение
		while (index < n && matrix[index][i] == 0)
			index++;

		if (index == n) // если существует ненулевой элемент
			continue; // определитель матрицы равен нулю

		if (index != i) {
			// цикл для замены строки диагонального элемента и элемента индексной строки
			for (int j = 0; j < n; j++)
				swap(matrix[index][j], matrix[i][j]);

			// знак определителя меняется, когда мы сдвигаем строки
			det = det * pow(-1, index - i);
		}

		// сохранение значений элементов диагональной строки
		for (int j = 0; j < n; j++)
			temp[j] = matrix[i][j];

		// обход каждой строки под диагональным элементом
		for (int j = i + 1; j < n; j++) {
			num1 = temp[i]; // значение диагонального элемента
			num2 = matrix[j][i]; // значение элемента следующей строки

			// обход каждого элемента строки и умножение на каждую строку
			for (int k = 0; k < n; k++)
				matrix[j][k] = (num1 * matrix[j][k]) - (num2 * temp[k]);

			total = total * num1; // Det(kA)=kDet(A);
		}
	}

	// умножение диагональных элементов для получения определителя
	for (int i = 0; i < n; i++)
		det = det * matrix[i][i];

	// Det(kA)/k=Det(A);
	return (det / total);
}

general_problem_t FromFileConvertToGeneral(string filename) {
	// открываем поток чтения из файла
	ifstream input_file_stream(filename);

	string cur_string;
	int N = 0; // кол-во переменных
	int N1 = 0; // кол-во полож. переменных
	int M1 = 0; // сколько ограничений в виде неравенств <= или >=
	double cur_number = 0; // буферное число
	column_t column_C; // вектор функции цели
	column_t column_B; // вектор правой части ограничений в виде равенств и неравенств
	matrix_t matrix_A; // матрица ограничений
	vector<int> N1_vec; // индексы полож. переменных 

	getline(input_file_stream, cur_string);
	// если файл не начинается с min_func: то завершаем программу
	if (cur_string != "min_func:") {
		cout << "Error reading file!" << endl;
		exit(-1);
	}
	else {
		getline(input_file_stream, cur_string);
		// открываем поток чтения из прочитанной строки
		stringstream str_stream(cur_string);

		// считываем функцию цели
		while (str_stream >> cur_number)
			column_C.push_back(cur_number);

		// определяем кол-во переменных
		N = column_C.size();

		//сколько будет столбцов у матрицы A
		matrix_A.resize(N);
	}

	getline(input_file_stream, cur_string);

	if (cur_string == ">=") {
		int cur_n = 0;
		getline(input_file_stream, cur_string);

		// считываем неравенства >= в матрицу A
		while (cur_string != "<=") {
			stringstream str_stream(cur_string);
			cur_n = 0;
			while (str_stream >> cur_number && cur_n != N) {
				matrix_A[cur_n].push_back(cur_number);
				++cur_n;
			}
			column_B.push_back(cur_number);
			++M1;
			getline(input_file_stream, cur_string);
		}
		getline(input_file_stream, cur_string);

		// считываем неравенства <= в матрицу A, превращая их в >=
		while (cur_string != "==") {
			stringstream str_stream(cur_string);
			cur_n = 0;
			while (str_stream >> cur_number && cur_n != N) {
				matrix_A[cur_n].push_back(-cur_number);
				++cur_n;
			}
			column_B.push_back(-cur_number);
			++M1;
			getline(input_file_stream, cur_string);
		}
		getline(input_file_stream, cur_string);

		// считываем равенства в матрицу A
		while (cur_string != "x_i >= 0") {
			stringstream str_stream(cur_string);
			cur_n = 0;
			while (str_stream >> cur_number && cur_n != N) {
				matrix_A[cur_n].push_back(cur_number);
				++cur_n;
			}
			column_B.push_back(cur_number);
			getline(input_file_stream, cur_string);
		}
		getline(input_file_stream, cur_string);
		stringstream str_stream(cur_string);

		// записываем индексы полож. переменных 
		int cur_index = 0;
		while (str_stream >> cur_index)
			N1_vec.push_back(cur_index);
		N1 = N1_vec.size();
	}
	// закрываем поток файла
	input_file_stream.close();

	/*
	здесь мы переставляем в матрице A и столбце C столбцы/координаты так,
	чтобы первые N1 столцов/координат соответствовали по порядку x_i, i in N1_vec
	*/
	matrix_t matrix_A_mod(N);
	column_t column_C_mod(N);
	vector<int> iterator(N);
	for (int i = 0; i < N; i++)
		iterator[i] = i;
	int i = 0;
	for (auto index : N1_vec) {
		matrix_A_mod[i] = matrix_A[index - 1];
		column_C_mod[i] = column_C[index - 1];
		iterator[index - 1] = -1;
		++i;
	}
	for (auto index : iterator) {
		if (index != -1) {
			matrix_A_mod[i] = matrix_A[index];
			column_C_mod[i] = column_C[index];
			++i;
		}
	}
	// очищаем вспомогательные векторы/матрицы
	N1_vec.clear();
	iterator.clear();
	matrix_A.clear();
	column_C.clear();

	general_problem_t problem = make_tuple(matrix_A_mod, column_B, column_C_mod, M1, N1);
	return problem;
}

canon_problem_t ConvertGeneralToCanon(general_problem_t& problem) {
	column_t column_C;
	column_t column_B;
	matrix_t matrix_A;
	int N1;
	int M1;
	tie(matrix_A, column_B, column_C, M1, N1) = problem;
	int N = column_C.size();
	int N2 = N - N1;

	/*
	расширяем матрицу A, чтобы в нее поместить ограничения,
	связанные с нововведенными полож. переменными;
	заполняем её (см. явные формулы в пособии)
	*/
	matrix_A.resize(N + N2 + M1);
	int M = column_B.size();
	int cur_M;
	for (int i = N; i < N + N2; ++i) {
		cur_M = 0;
		while (cur_M != M) {
			matrix_A[i].push_back(-matrix_A[i - N2][cur_M]);
			++cur_M;
		}
	}
	cur_M = 0;
	for (int i = N + N2; i < N + N2 + M1; ++i) {
		matrix_A[i].resize(M);
		if (cur_M != M1) {
			matrix_A[i][cur_M] = -1;
			++cur_M;
		}
	}

	/*
	преобразовываем функцию цели, учитывая нововведенные переменные
	(см. явные формулы в пособии)
	*/
	column_C.resize(N + N2 + M1);
	for (int i = N; i < N + N2; ++i)
		column_C[i] = -column_C[i - N2];

	/*
	преобразываем систему ограничений: b >= 0
	*/
	int columns_A = matrix_A.size();
	for (int i = 0; i < M; ++i) {
		if (column_B[i] < 0) {
			column_B[i] *= -1;
			for (int j = 0; j < columns_A; ++j)
				matrix_A[j][i] *= -1;
		}
	}

	canon_problem_t canon_problem = make_tuple(matrix_A, column_B, column_C);
	return canon_problem;
}

void TransposeQuadMatrix(matrix_t& matrix) {
	int M = matrix.size();
	for (int i = 0; i < M; ++i) {
		for (int j = 0; j < i; ++j)
			swap(matrix[i][j], matrix[j][i]);
	}
}

matrix_t TransformationMatrix(matrix_t matrix) {
	column_t columnTransform;
	matrix_t transform;
	int n = size(matrix);
	int m = size(matrix[1]);
	columnTransform.resize(n);
	transform.resize(m);
	for (int k = 0; k < m; k++)
		transform[k] = columnTransform;

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++)
			transform[j][i] = -matrix[i][j];
	}

	return transform;
}

general_problem_t GetDualLinearProblem(general_problem_t& problem) {
	matrix_t A;
	column_t b, c;
	int M1, N1;
	tie(A, b, c, M1, N1) = problem;
	// инциализация двойственной матрицы и двойственных столбцов 
	matrix_t dualA;
	column_t dualb, dualc;

	// получение матрицы A двойственной задачи 
	dualA = TransformationMatrix(A);

	//столбец b двойственной задачи  - это столбец c прямой задачи, где все значения берутся со знаком - 
	/* все элементы столбца нужно умножить на - 1, чтобы свести полученную двойственную задаче к общей задаче ЛП
	(получить нужный знак неравенства(>= )) */
	dualb.resize(size(c));
	for (int i = 0; i < size(c); i++)
		dualb[i] = -c[i];

	// столбец c двойственной задачи - это столбец b прямой задачи
	dualc.resize(size(b));
	for (int j = 0; j < size(b); j++)
		dualc[j] = b[j];

	/* числа N1 и M1 - это количество переменных с ограничениями на знак
	и количество неравенств в системе ограничений соответственно */
	int dualN1 = M1;
	int dualM1 = N1;

	general_problem_t dualProblem = make_tuple(dualA, dualb, dualc, dualM1, dualN1);

	return dualProblem;
}

solution_t IteratingThroughExtremePoints(canon_problem_t& problem) {
	matrix_t A;
	column_t b, c;
	tie(A, b, c) = problem;
	int n = size(c), m = size(b);

	/* проверка корректности задачи */
	if (n <= m) {
		cout << "Incorred problem";
		exit(-1);
	}

	/* получение всех комбинаций из n (количество переменных) по m (количество ограничений) */
	combinations_t combinations = CombinationsWithoutRepetitions(m, n);

	column_t maxPoint; // вектор, в котором достигается оптимальное решение
	maxPoint.resize(n);
	double max = -1e20; // максимальное значение функции цели
	comb_t optComb; // базисная комбинация, при которой достигается оптимальное решение
	optComb.resize(m);


	/* перебор всех возможных комбинаций сочетаний */
	for (auto& comb : combinations) {

		double valueInPoint; // значение функции цели при данной комбинации
		column_t X; // опорный вектор при данной комбинации 
		X.resize(n);
		column_t columnA1; //матрица для решения системы уравнений, из которой отброшены столбцы , номера которых не входят в данную комбинацию 

		/*размерность матрицы A1 m на m*/
		columnA1.resize(m);
		matrix_t A1;
		A1.resize(m);

		for (int k = 0; k < m; k++)
			A1[k] = columnA1;

		column_t solveSystem; //решение СЛАУ с матрицей A1 и вектором b, имеет размерность 
		solveSystem.resize(m);


		int numb = 0; // счетсчик для задания компонент матрицы A1

		/* задание компонент матрицы A1, берутся те столбцы, номера которых есть в данной комбинации */
		for (int i = 0; i < n; i++) {
			if (IsNumberInCombination(i + 1, comb)) {
				for (int j = 0; j < m; j++)
					A1[j][numb] = A[i][j];

				numb++;
			}
		}

		/* те компоненты опорного вектора, номера которых не входят в данную комбинацию принимают нулевые значения */
		for (int i = 0; i < n; i++) {
			if (!IsNumberInCombination(i + 1, comb))
				X[i] = 0;
		}

		// проверка невырожденности матрицы 
		if (Determinant(A1) != 0) { 
			/*решение СЛАУ методом вращений*/
			solveSystem = RotationMethod(A1, b);

			// проверка на неотрицательность компонент решения СЛАУ (иначе вектор не будет допустимым)
			if (NonNegativityOfVector(solveSystem)) {
				int numbInSolv = 0; //счетчик 

			/*те компоненты опорного вектора, номера которых входят в данную комбинацию принимают соответствующие значения решения СЛАУ*/
				for (int k = 0; k < n; k++) {
					if (IsNumberInCombination(k + 1, comb)) {
						X[k] = solveSystem[numbInSolv];
						numbInSolv++;
					}
				}

				/*вычисление значения целевой функции при данном опорном векторе*/
				valueInPoint = MultipliedVectors(c, X);

				/*если полученное значение функции цели меньше минимума, то оно становится минимумом,
				вектору, в котором достигается минимум присваивается значение X*/
				if (valueInPoint > max) {
					maxPoint = X;
					max = valueInPoint;
					optComb = comb;
				}
			}
		}
	}
	/*создание кортежа из вектора, в котором достигается оптимальное решение и оптимального решения*/
	solution_t solving = make_tuple(maxPoint, max, optComb);
	return solving;
}

column_t SolvingDualProblem(canon_problem_t problem, column_t X, comb_t optBasis) {
	matrix_t A;
	column_t b, c;
	tie(A, b, c) = problem;

	int n = size(c);
	int m = size(b);

	column_t Y; // решение двойственной задачи 

	column_t columnBasis; // столбец базиса
	columnBasis.resize(m);

	matrix_t basis; //базис оптимального решения 
	basis.resize(m);

	for (int i = 0; i < m; i++)
		basis[i] = columnBasis;

	column_t cB; // вектор компонент целевой функцию, соответствующий базису оптимального решения
	cB.resize(m);

	int numb = 0; //счетчик

	for (int i = 0; i < n; i++) {
		if (IsNumberInCombination(i + 1, optBasis)) {
			for (int j = 0; j < m; j++)
				basis[numb][j] = A[i][j];

			cB[numb] = c[i];
			numb++;
		}
	}

	matrix_t invBasis = InverseMatrix(basis); // инвертируем матрицу базиса

	Y = MultipliedMatrixAndColumn(invBasis, cB);

	double result = 0;
	for (int i = 0; i < size(Y); i++)
		result += Y[i] * b[i];

	cout << endl << " " << result << " " << endl;

	return Y;
}

void GetCurMatrix(matrix_t& A, matrix_t& cur_matrix, comb_t cur_columns) {
	int j = 0;
	int M = cur_matrix.size();
	for (auto cur_index : cur_columns) {
		for (int i = 0; i < M; ++i)
			cur_matrix[i][j] = A[i][cur_index - 1];

		++j;
	}
}

void CalculateNextSupportVector(comb_t& Nk_indexes, column_t& u_k, column_t& cur_X, int N_number) {
	double theta = DBL_MAX;
	double potential_theta = 0;
	int i_k = -1;
	for (auto temp : Nk_indexes) {
		if (u_k[temp - 1] <= 0)
			continue;
		else {
			potential_theta = cur_X[temp - 1] / u_k[temp - 1];
			if (potential_theta < theta) {
				theta = potential_theta;
				i_k = temp;
			}
		}
	}
	// считаем x_{k+1} = x_k - \theta * u_k
	for (int i = 0; i < N_number; ++i) {
		cur_X[i] = cur_X[i] - theta * u_k[i];
	}
	// cur_X[i_k - 1] = 0;
}

bool IsNeedChangeBasis(comb_t& Nk_indexes, comb_t& N_plus_indexes, column_t& u_k) {
	for (auto temp : Nk_indexes) {
		if (IsNumberInCombination(temp, N_plus_indexes))
			continue;
		else {
			if (u_k[temp - 1] > 0)
				return true;
		}
	}
	return false;
}

SimplexState IterSimplex(matrix_t& main_A, int M_number, comb_t& basis, column_t& main_C, column_t& cur_X, column_t& cur_Y) {
	// инициализация буферной матрицы для нахождения опорного вектора
	matrix_t A_M_Nk(M_number);
	for (int i = 0; i < M_number; ++i)
		A_M_Nk[i].resize(M_number);

	int N_number = main_C.size();

	// формируем набор индексов N+ 
	vector<int> N_plus_indexes;
	for (int i = 0; i < N_number; ++i) {
		if (cur_X[i] > 0)
			N_plus_indexes.push_back(i + 1);
	}

	comb_t Nk_indexes(M_number);
	bool is_cur_X_singular;

	// если опорный вектор невырожден (|N+| = |M|), то формируем квадратную матрицу по индексам столбцов из N+
	if (M_number == N_plus_indexes.size()) {
		GetCurMatrix(main_A, A_M_Nk, N_plus_indexes);
		Nk_indexes = N_plus_indexes;
		is_cur_X_singular = false;
	}
	else { 
		/* если же вырожден, то нужно к текущим столбцам N+ добавлять столбцы из N0, чтобы текущая матрица была невырождена */
		is_cur_X_singular = true;

		// инициализируем первое сочетание
		for (int i = 0; i < M_number; ++i)
			Nk_indexes[i] = i + 1;

		bool trigger = false;

		while (true) {
			// проверяем: все ли столбцы из A[][N+] на месте
			for (int i = 0; i < N_plus_indexes.size(); ++i) {
				if (!IsNumberInCombination(N_plus_indexes[i], Nk_indexes)) {
					// не нашли какой-то столбец - генерим следующее сочетание, вызывая триггер и выходя из перебора индексов
					trigger = NextSet(Nk_indexes, N_number, M_number);
					break;
				}
			}
			// триггер сработан - выключаем его и идем дальше проверять только что измененное сочетание
			if (trigger == true) {
				trigger = false;
				continue;
			}
			// если все же все столбцы на месте, то формируем квадратную матрицу из исходной
			GetCurMatrix(main_A, A_M_Nk, Nk_indexes);
			// и проверяем ее определитель
			if (Determinant(A_M_Nk) != 0) {
				break;
			}
			else { // если вырожденная - генерим следующее сочетание и идем дальше проверять только что измененное сочетание
				NextSet(Nk_indexes, N_number, M_number);
				continue;
			}
		}
	}

	// получим текущий вектор C[N_k]
	column_t C_Nk(M_number);
	int i = 0;
	for (auto temp : Nk_indexes) {
		C_Nk[i] = main_C[temp - 1];
		++i;
	}

	/* найдем y^T [M] * A[M][Nk] = c^T[M] */
	/* для этого транспонируем матрицу 
	   A[Nk][M] * y[M] = c[M]	*/
	TransposeQuadMatrix(A_M_Nk);
	// и решим СЛАУ методом вращения
	cur_Y = RotationMethod(A_M_Nk, C_Nk);

	// Найдем обратную к A[M][Nk] матрицу B[Nk][M]
	matrix_t B_Nk_M = InverseMatrix(A_M_Nk);

	// вернем матрицу к исходному состоянию
	TransposeQuadMatrix(A_M_Nk);

	// сформируем массив индексов L_k 
	comb_t L_k;
	int cur_l = 1;
	for (auto temp : Nk_indexes) {
		while (cur_l < temp) {
			L_k.push_back(cur_l);
			++cur_l;
		}
	}

	// сделаем cur_X[L_k] = 0
	for (auto temp : L_k) {
		cur_X[temp - 1] = 0.0;
	}

	// получим d_k[L_k]
	column_t d_k(L_k.size());
	i = 0;
	double cur_sum = 0;
	for (auto temp : L_k) {
		d_k[i] = main_C[temp - 1];
		for (int j = 0; j < M_number; ++j)
			cur_sum += cur_Y[j] * main_A[j][temp - 1];
		d_k[i] -= cur_sum;
		cur_sum = 0;
		++i;
	}

	// если d_k не отрицательно, то выполнены н. и д. условия оптимальности, можем выдавать оптимальное значение
	if (NonNegativityOfVector(d_k))
		return SimplexState::OPTIMAL;

	/* если нет, то начинаем строить u_k[N] */
	// найдем j_k
	int j_k = -1;
	for (int i = 0; i < d_k.size(); ++i) {
		if (d_k[i] < 0) {
			j_k = i + 1;
			break;
		}
	}
	column_t u_k(N_number);
	// u_k[Nk] = B[Nk][M] * A[M][j_k]
	for (auto temp : Nk_indexes) {
		for (int i = 0; i < M_number; ++i) 
			u_k[temp - 1] += B_Nk_M[temp - 1][i] * A_M_Nk[i][j_k - 1];
	}
	// u_k[j_k] = -1
	u_k[j_k - 1] = -1;
	// u_k[L_k \ j_k] = 0
	// уже занулены при инициализации

	// если u_k не положительно, то функция не ограничена
	if (NonPositivityOfVector(u_k))
		return SimplexState::UNLIMITED;

	/* если x_k[N] невырожденный, то существует \theta, считаем x_{k+1} и начинаем другую итерацию */
	if (is_cur_X_singular == false) {
		CalculateNextSupportVector(Nk_indexes, u_k, cur_X, N_number);
		return SimplexState::NEXT;
	}
	else { 
		/* если вырожденный, то рассматриваем два случая */
		// если u_k[Nk \ Nk+],то считаем \theta
		if (IsNeedChangeBasis(Nk_indexes, N_plus_indexes, u_k) == false) {
			CalculateNextSupportVector(Nk_indexes, u_k, cur_X, N_number);
			return SimplexState::NEXT;
		}
		else {
			// иначе меняем базис ВОЗМОЖНО ЗАЦИКЛИВАНИЕ

			return SimplexState::NEXT;
		}

	}
}

double SimplexMethod(canon_problem_t& problem, column_t& cur_X) {
	// распаковка кортежа
	matrix_t matrix_A0;
	column_t column_B;
	column_t column_C;
	tie(matrix_A0, column_B, column_C) = problem;

	// транспонирование матрицы для удобства
	int N = matrix_A0.size();
	int M = matrix_A0[0].size();
	matrix_t matrix_A(M);
	for (int i = 0; i < M; ++i) {
		matrix_A[i].resize(N);
		for (int j = 0; j < N; ++j) {
			matrix_A[i][j] = matrix_A0[j][i];
		}
	}
	for (int i = 0; i < M; ++i)
		matrix_A0[i].clear();
	matrix_A0.clear();

	// проверка: является ли матрица полноранговой
	if (RankMatrix(matrix_A) != M) {
		cout << "Матрица ограничений не является полноранговой." << endl << "Досрочное завершение работы." << endl;
		exit(-1);
	}

	column_t cur_Y(M);

	SimplexState state = SimplexState::NEXT;
	while (true) {
		// переходим к следующему опорному вектору
		if (state == SimplexState::NEXT)
			state = IterSimplex(matrix_A, M, column_C, cur_X, cur_Y);

		// если функция не ограничена, то выводим -"бесконечность"
		else if (state == SimplexState::UNLIMITED)
			return -DBL_MAX;

		// если нашли оптимальное решение, то выводим значение функции цели
		else if (state == SimplexState::OPTIMAL) {
			double result = 0;
			for (int i = 0; i < N; ++i)
				result += column_C[i] * cur_X[i];
			return result;
		} 
	}
}

column_t GetInitialApprox(canon_problem_t& problem) {
	// распаковка кортежа
	matrix_t matrix_A;
	column_t column_B;
	column_t column_C;
	tie(matrix_A, column_B, column_C) = problem;

	/* формулируем критерий качества */
	int b_size = column_B.size();
	int c_size = column_C.size();
	column_C.resize(b_size + c_size);
	for (int i = 0; i < c_size; ++i)
		column_C[i] = 0;
	for (int i = c_size; i < b_size + c_size; ++i)
		column_C[i] = 1;

	/* расширяем матрицу A */
	matrix_A.resize(matrix_A.size() + b_size);
	for (int i = matrix_A.size() - b_size; i < matrix_A.size(); ++i)
		matrix_A[i][i - (matrix_A.size() - b_size)] = 1;

	/* после преобразовании задачи у нас b[M] >= 0 */

	/* строим начальное приближение */
	column_t first_approx(column_C.size());
	for (int i = 0; i < c_size; ++i)
		first_approx[i] = 0;
	for (int i = c_size; i < b_size + c_size; ++i)
		first_approx[i] = column_B[i - c_size];

	canon_problem_t help_problem = make_tuple(matrix_A, column_B, column_C);

	cout << "Начало поиска начального приближения методом искусственного базиса:" << endl;

	SimplexMethod(help_problem, first_approx);

	cout << "Конец поиска начального приближения методом искусственного базиса." << endl;

	/* если хотя бы одна компонента y[m] положительна, то исходная задача не имеет ни одной допустимой точки */
	for (int i = c_size; i < c_size + b_size; ++i) {
		if (first_approx[i] > 0) {
			cout << "Исходная задача не имеет ни одной допустимой точки." << endl << "Досрочное завершение работы." << endl;
			exit(-1);
		}
	}
	/* вычленяем из вектора x[N] - опорный вектор для исходной задачи */
	for (int i = 0; i < c_size; ++i) {
		first_approx[i] = first_approx[i + c_size];
	}
	first_approx.resize(c_size);

	return first_approx;
}
