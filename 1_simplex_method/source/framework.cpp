#include "framework.h"
using namespace std;

#pragma warning(disable:4996)
#pragma warning(disable:4267)


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
	temp.resize(n); //������������� ������
	comb.resize(m); //������ ����� ����� - ���� �� ��������� ��������� ����������� m 

	int num = 0; //���������� - �������

	/*������������� ������� ���������*/
	combinations_t matrixCombinations;

	matrixCombinations.resize(1);

	for (int i = 0; i < n; i++)
		temp[i] = i + 1;

	/*��������� ������� ���������*/
	for (int j = 0; j < m; j++)
		comb[j] = temp[j];

	matrixCombinations[num] = comb;

	/*��������� ����������� ���������*/
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

comb_t GetDifferenceSets(comb_t& A, comb_t& B) {
	comb_t res;
	for (auto temp : A)
		if (!IsNumberInCombination(temp, B))
			res.push_back(temp);
	return res;
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

	column_t X(n); ///��������� ������������

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++)
			X[i] += M[i][j] * c[j];
	}

	return X;
}

int RankMatrix(matrix_t& A) {
	return 0;
}

matrix_t InverseMatrix(matrix_t A) {
	int n = size(A);

	column_t Ecolumn; // ������ ��������� �������
	Ecolumn.resize(n);

	matrix_t E; // ������������� ��������� �������
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

	// ��������������� ��������� ������� �������� �������� �������� � ������� � 
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
	// ������ ���
	for (int i = 0; i < n - 1; i++) {
		temp = A[i][i];
		for (int j = i + 1; i < n; i++) {

			for (int k = 0; k < n; k++)
				A[k][j] = A[k][j] - A[i][j] / temp;

			b[j] = b[j] - b[i] / temp;
		}
	}

	// �������� �����������
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
	column_t temp; // ��������� ������ ��� �������� �������
	temp.resize(n + 1);

	// ���� ��� ������ ������������ ���������
	for (int i = 0; i < n; i++) {
		index = i; // ������������� �������

		// ����� �������, ������� ����� ��������� ��������
		while (index < n && matrix[index][i] == 0)
			index++;

		if (index == n) // ���� ���������� ��������� �������
			continue; // ������������ ������� ����� ����

		if (index != i) {
			// ���� ��� ������ ������ ������������� �������� � �������� ��������� ������
			for (int j = 0; j < n; j++)
				swap(matrix[index][j], matrix[i][j]);

			// ���� ������������ ��������, ����� �� �������� ������
			det = det * pow(-1, index - i);
		}

		// ���������� �������� ��������� ������������ ������
		for (int j = 0; j < n; j++)
			temp[j] = matrix[i][j];

		// ����� ������ ������ ��� ������������ ���������
		for (int j = i + 1; j < n; j++) {
			num1 = temp[i]; // �������� ������������� ��������
			num2 = matrix[j][i]; // �������� �������� ��������� ������

			// ����� ������� �������� ������ � ��������� �� ������ ������
			for (int k = 0; k < n; k++)
				matrix[j][k] = (num1 * matrix[j][k]) - (num2 * temp[k]);

			total = total * num1; // Det(kA)=kDet(A);
		}
	}

	// ��������� ������������ ��������� ��� ��������� ������������
	for (int i = 0; i < n; i++)
		det = det * matrix[i][i];

	// Det(kA)/k=Det(A);
	return (det / total);
}

general_problem_t FromFileConvertToGeneral(string filename) {
	// ��������� ����� ������ �� �����
	ifstream input_file_stream(filename);

	string cur_string;
	int N = 0; // ���-�� ����������
	int N1 = 0; // ���-�� �����. ����������
	int M1 = 0; // ������� ����������� � ���� ���������� <= ��� >=
	double cur_number = 0; // �������� �����
	column_t column_C; // ������ ������� ����
	column_t column_B; // ������ ������ ����� ����������� � ���� �������� � ����������
	matrix_t matrix_A; // ������� �����������
	vector<int> N1_vec; // ������� �����. ���������� 

	getline(input_file_stream, cur_string);
	// ���� ���� �� ���������� � min_func: �� ��������� ���������
	if (cur_string != "min_func:") {
		cout << "Error reading file!" << endl;
		exit(-1);
	}
	else {
		getline(input_file_stream, cur_string);
		// ��������� ����� ������ �� ����������� ������
		stringstream str_stream(cur_string);

		// ��������� ������� ����
		while (str_stream >> cur_number)
			column_C.push_back(cur_number);

		// ���������� ���-�� ����������
		N = column_C.size();

		//������� ����� �������� � ������� A
		matrix_A.resize(N);
	}

	getline(input_file_stream, cur_string);

	if (cur_string == ">=") {
		int cur_n = 0;
		getline(input_file_stream, cur_string);

		// ��������� ����������� >= � ������� A
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

		// ��������� ����������� <= � ������� A, ��������� �� � >=
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

		// ��������� ��������� � ������� A
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

		// ���������� ������� �����. ���������� 
		int cur_index = 0;
		while (str_stream >> cur_index)
			N1_vec.push_back(cur_index);
		N1 = N1_vec.size();
	}
	// ��������� ����� �����
	input_file_stream.close();

	/*
	����� �� ������������ � ������� A � ������� C �������/���������� ���,
	����� ������ N1 �������/��������� ��������������� �� ������� x_i, i in N1_vec
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
	
	// ������� ��������������� �������/�������
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
	��������� ������� A, ����� � ��� ��������� �����������,
	��������� � �������������� �����. �����������;
	��������� � (��. ����� ������� � �������)
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
	��������������� ������� ����, �������� ������������� ����������
	(��. ����� ������� � �������)
	*/
	column_C.resize(N + N2 + M1);
	for (int i = N; i < N + N2; ++i)
		column_C[i] = -column_C[i - N2];

	/*
	������������� ������� �����������: b >= 0
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

void GetCurMatrix(matrix_t& A, matrix_t& cur_matrix, comb_t& cur_columns) {
	int j = 0;
	int M = cur_matrix.size();
	for (auto cur_index : cur_columns) {
		for (int i = 0; i < M; ++i)
			cur_matrix[i][j] = A[i][cur_index - 1];

		++j;
	}
}
