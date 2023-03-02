#include "transformations.h"
using namespace std;

general_problem_t fromFileConvertToGeneral(string filename) {
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

canon_problem_t convertGeneralToCanon(general_problem_t& problem) {
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
