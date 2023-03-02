#include "simplex.h"
using namespace std;

// решение СЛАУ
column_t RotationMethod(matrix_t A, column_t b)
{
	long double c, s;
	int sz = size(b);
	int k = 0;
	column_t X;
	X.resize(sz);
	for (int i = 0; i < sz - 1; i++)
	{
		for (int j = i + 1; j < sz; j++)
		{
			c = A[i][k] / sqrt(A[i][k] * A[i][k] + A[j][k] * A[j][k]);
			s = A[j][k] / sqrt(A[i][k] * A[i][k] + A[j][k] * A[j][k]);
			for (int m = k; m < sz; m++)
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

bool iterSimplexMethod(matrix_t& simplex_table, column_t& cur_solution) {
	int columns_simplex = simplex_table.size();
	int rows_simplex = simplex_table[0].size();

	vector<int> index_solution(columns_simplex - 1);
	int t = 1;
	for (int i = 0; i < columns_simplex - 1; ++i) {
		if (cur_solution[i] != 0) {
			index_solution[i] = t;
			++t;
		}
	}

	// ищем индекс максимально отрицательной компоненты в z столбце
	int leading_column_index = -1; double max_negative = 0;
	for (int i = 0; i < columns_simplex - 1; ++i) {
		if (simplex_table[i][0] < 0) {
			if (simplex_table[i][0] < max_negative) {
				max_negative = simplex_table[i][0];
				leading_column_index = i;
			}
		}
	}
	if (leading_column_index == -1)
		return true;

	// ищем минимум полож. соотношений чисел из столбца "решение" к коэффициентам при найденной переменной
	int leading_row_index = -1; double min_positive = INT_MAX;
	for (int i = 1; i < rows_simplex; ++i) {
		if (simplex_table[leading_column_index][i] > 0) {
			if (simplex_table[columns_simplex - 1][i] / simplex_table[leading_column_index][i] < min_positive) {
				min_positive = simplex_table[columns_simplex - 1][i] / simplex_table[leading_column_index][i];
				leading_row_index = i;
			}
		}
	}

	// преобразовываем симплекс таблицу
	double leading_elem = simplex_table[leading_column_index][leading_row_index];
	// новая ведущая строка
	for (int i = 0; i < columns_simplex; ++i)
		simplex_table[i][leading_row_index] /= leading_elem;
	// остальные строки без ведущей
	for (int i = 0; i < rows_simplex; ++i) {
		if (i == leading_row_index)
			continue;
		leading_elem = simplex_table[leading_column_index][i];
		for (int j = 0; j < columns_simplex; ++j)
			simplex_table[j][i] = simplex_table[j][i] - leading_elem * simplex_table[j][leading_row_index];
	}

	// добавляем новую базисную переменную, удалив одну старую
	int delete_index = -1;
	for (int i = 0; i < columns_simplex - 1; ++i) {
		if (index_solution[i] == leading_row_index) {
			delete_index = i;
			break;
		}
	}
	cur_solution[delete_index] = 0;
	cur_solution[leading_column_index] = 1;

	// находим решение СЛАУ для поиска очередного решения
	int p = 0;
	matrix_t A(rows_simplex - 1);
	for (int i = 0; i < rows_simplex - 1; ++i) {
		A[i].resize(rows_simplex - 1);
		while (cur_solution[p] == 0)
			++p;
		for (int j = 0; j < rows_simplex - 1; ++j) {
			A[i][j] = simplex_table[p][j + 1];
		}
		++p;
	}
	column_t b(rows_simplex - 1);
	for (int i = 0; i < rows_simplex - 1; ++i)
		b[i] = simplex_table[columns_simplex - 1][i + 1];

	column_t not_full_solution = RotationMethod(A, b);

	int j = 0;
	for (int i = 0; i < columns_simplex - 1; ++i) {
		if (cur_solution[i] != 0) {
			cur_solution[i] = not_full_solution[j];
			++j;
		}
	}

	// проверяем функцию цели
	bool flag = true;
	for (int i = 0; i < columns_simplex - 1; ++i) {
		if (simplex_table[i][0] < 0) {
			flag = false;
			break;
		}
	}
	if (flag == true)
		return true;

	return false;
}

column_t simplexMethod(canon_problem_t& problem, column_t& cur_solution) {
	// распаковка кортежа
	matrix_t matrix_A;
	column_t column_B;
	column_t column_C;
	tie(matrix_A, column_B, column_C) = problem;

	// инициализируем симплекс таблицы
	int c_size = column_C.size();
	int columns_simplex = c_size + 1;
	int rows_simplex = matrix_A[0].size() + 1;
	matrix_t simplex_table(columns_simplex);
	for (int i = 0; i < columns_simplex; ++i)
		simplex_table[i].resize(rows_simplex);

	// заполняем симплекс таблицу
	for (int j = 0; j < columns_simplex - 1; ++j)
		simplex_table[j][0] = column_C[j];
	simplex_table[columns_simplex - 1][0] = 0;

	for (int i = 1; i < rows_simplex; ++i) {
		for (int j = 0; j < columns_simplex - 1; ++j)
			simplex_table[j][i] = matrix_A[j][i - 1];
		simplex_table[columns_simplex - 1][i] = column_B[i - 1];
	}

	while (!iterSimplexMethod(simplex_table, cur_solution));

	return cur_solution;
}

double getOptimalValue(canon_problem_t& problem, column_t solution) {
	matrix_t A;
	column_t B, C;
	tie(A, B, C) = problem;
	double sum = 0;
	for (int i = 0; i < C.size(); ++i)
		sum += C[i] * solution[i];
	return sum;
}

column_t getInitialApprox(canon_problem_t& problem) {
	// распаковка кортежа
	matrix_t matrix_A;
	column_t column_B;
	column_t column_C;
	tie(matrix_A, column_B, column_C) = problem;

	/* формулируем критерий качества */
	int b_size = column_B.size();
	for (int i = 0; i < b_size; ++i)
		column_C[i] = 1;

	/* расширяем матрицу A */
	matrix_A.resize(matrix_A.size() + b_size);
	for (int i = matrix_A.size() - b_size; i < matrix_A.size(); ++i)
		matrix_A[i][i - (matrix_A.size() - b_size)] = 1;

	column_t first_approx(column_C.size());
	for (int i = 0; i < matrix_A.size() - b_size; ++i)
		first_approx[i] = 0;
	for (int i = matrix_A.size() - b_size; i < matrix_A.size(); ++i)
		first_approx[i] = column_B[i - (matrix_A.size() - b_size)];

	canon_problem_t help_problem = make_tuple(matrix_A, column_B, column_C);
	
	simplexMethod(help_problem, first_approx);

	return first_approx;
}
