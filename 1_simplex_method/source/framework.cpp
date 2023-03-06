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

/*int RankMatrix(matrix_t& m) {
	double tmp;
	unsigned i;
	unsigned x, y;

	int width = m[0].size();
	int height = m.size();
	// �������� �������
	unsigned endi = (width > height) ? height : width;
	// ���� �� ���� ������� ���������
	for (i = 0; i < endi; i++) {
		// ���� ������� �� ��������� ����� 0, �� ���� �� ������� ������� � �������
		if (m[i][i] == 0) {
			// ���� ��� �������� �������, ��������� ����
			if (!m.search(0, false, y, x, i, i)) break;
			// ������ i-�� ������ � y-��
			if (i != y) m.swaprows(i, y);
			// ������ i-�� ������� � x-��
			if (i != x) m.swapcolumns(i, x);
			// ����� �������, � m[i][i], ������ ��������� �������.
		}

		// ������� ������� m[i][i]
		tmp = m[i][i];
		for (x = i; x < width; x++) {
			m[i][x] = m[i][x] / tmp;
		}
		// ����� ������� m[i][i] ������ ����� 1
		// �������� ��� �������� ������� ��� (i, i)-�� � ������ �� ����,
		// ��� ������ ��������� � ���. �������������
		for (y = i + 1; y < height; y++) {
			tmp = m[y][i];
			for (x = i; x < width; x++)
				m[y][x] -= (m[i][x] * tmp);
		}
		for (x = i + 1; x < width; x++) {
			tmp = m[i][x];
			for (y = i; y < height; y++)
				m[y][x] -= (m[y][i] * tmp);
		}
	}

	// ������� ������� �������� �� ������� ���������
	unsigned cnt = 0;
	for (i = 0; i < endi; i++)
		if (m[i][i] == 0) break;
		else cnt++;
	if (!cnt) cnt++;
	return cnt;
}; // rang()*/

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
	// ������������ ������������ ������� � ������������ �������� 
	matrix_t dualA;
	column_t dualb, dualc;

	// ��������� ������� A ������������ ������ 
	dualA = TransformationMatrix(A);

	//������� b ������������ ������  - ��� ������� c ������ ������, ��� ��� �������� ������� �� ������ - 
	/* ��� �������� ������� ����� �������� �� - 1, ����� ������ ���������� ������������ ������ � ����� ������ ��
	(�������� ������ ���� �����������(>= )) */
	dualb.resize(size(c));
	for (int i = 0; i < size(c); i++)
		dualb[i] = -c[i];

	// ������� c ������������ ������ - ��� ������� b ������ ������
	dualc.resize(size(b));
	for (int j = 0; j < size(b); j++)
		dualc[j] = b[j];

	/* ����� N1 � M1 - ��� ���������� ���������� � ������������� �� ����
	� ���������� ���������� � ������� ����������� �������������� */
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

	/* �������� ������������ ������ */
	if (n <= m) {
		cout << "Incorred problem";
		exit(-1);
	}

	/* ��������� ���� ���������� �� n (���������� ����������) �� m (���������� �����������) */
	combinations_t combinations = CombinationsWithoutRepetitions(m, n);

	column_t maxPoint; // ������, � ������� ����������� ����������� �������
	maxPoint.resize(n);
	double max = -1e20; // ������������ �������� ������� ����
	comb_t optComb; // �������� ����������, ��� ������� ����������� ����������� �������
	optComb.resize(m);


	/* ������� ���� ��������� ���������� ��������� */
	for (auto& comb : combinations) {

		double valueInPoint; // �������� ������� ���� ��� ������ ����������
		column_t X; // ������� ������ ��� ������ ���������� 
		X.resize(n);
		column_t columnA1; //������� ��� ������� ������� ���������, �� ������� ��������� ������� , ������ ������� �� ������ � ������ ���������� 

		/*����������� ������� A1 m �� m*/
		columnA1.resize(m);
		matrix_t A1;
		A1.resize(m);

		for (int k = 0; k < m; k++)
			A1[k] = columnA1;

		column_t solveSystem; //������� ���� � �������� A1 � �������� b, ����� ����������� 
		solveSystem.resize(m);


		int numb = 0; // �������� ��� ������� ��������� ������� A1

		/* ������� ��������� ������� A1, ������� �� �������, ������ ������� ���� � ������ ���������� */
		for (int i = 0; i < n; i++) {
			if (IsNumberInCombination(i + 1, comb)) {
				for (int j = 0; j < m; j++)
					A1[j][numb] = A[i][j];

				numb++;
			}
		}

		/* �� ���������� �������� �������, ������ ������� �� ������ � ������ ���������� ��������� ������� �������� */
		for (int i = 0; i < n; i++) {
			if (!IsNumberInCombination(i + 1, comb))
				X[i] = 0;
		}

		// �������� ��������������� ������� 
		if (Determinant(A1) != 0) { 
			/*������� ���� ������� ��������*/
			solveSystem = RotationMethod(A1, b);

			// �������� �� ����������������� ��������� ������� ���� (����� ������ �� ����� ����������)
			if (NonNegativityOfVector(solveSystem)) {
				int numbInSolv = 0; //������� 

			/*�� ���������� �������� �������, ������ ������� ������ � ������ ���������� ��������� ��������������� �������� ������� ����*/
				for (int k = 0; k < n; k++) {
					if (IsNumberInCombination(k + 1, comb)) {
						X[k] = solveSystem[numbInSolv];
						numbInSolv++;
					}
				}

				/*���������� �������� ������� ������� ��� ������ ������� �������*/
				valueInPoint = MultipliedVectors(c, X);

				/*���� ���������� �������� ������� ���� ������ ��������, �� ��� ���������� ���������,
				�������, � ������� ����������� ������� ������������� �������� X*/
				if (valueInPoint > max) {
					maxPoint = X;
					max = valueInPoint;
					optComb = comb;
				}
			}
		}
	}
	/*�������� ������� �� �������, � ������� ����������� ����������� ������� � ������������ �������*/
	solution_t solving = make_tuple(maxPoint, max, optComb);
	return solving;
}

column_t SolvingDualProblem(canon_problem_t problem, column_t X, comb_t optBasis) {
	matrix_t A;
	column_t b, c;
	tie(A, b, c) = problem;

	int n = size(c);
	int m = size(b);

	column_t Y; // ������� ������������ ������ 

	column_t columnBasis; // ������� ������
	columnBasis.resize(m);

	matrix_t basis; //����� ������������ ������� 
	basis.resize(m);

	for (int i = 0; i < m; i++)
		basis[i] = columnBasis;

	column_t cB; // ������ ��������� ������� �������, ��������������� ������ ������������ �������
	cB.resize(m);

	int numb = 0; //�������

	for (int i = 0; i < n; i++) {
		if (IsNumberInCombination(i + 1, optBasis)) {
			for (int j = 0; j < m; j++)
				basis[numb][j] = A[i][j];

			cB[numb] = c[i];
			numb++;
		}
	}

	matrix_t invBasis = InverseMatrix(basis); // ����������� ������� ������

	Y = MultipliedMatrixAndColumn(invBasis, cB);

	double result = 0;
	for (int i = 0; i < size(Y); i++)
		result += Y[i] * b[i];

	cout << endl << " " << result << " " << endl;

	return Y;
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

void CalculateNextSupportVector(comb_t& Nk_indexes, comb_t& Nk_plus_indexes, column_t& u_k, column_t& cur_X, int N_number, int j_k) {
	double theta = DBL_MAX;
	double potential_theta = 0;
	int i_k = -1;
	for (auto temp : Nk_plus_indexes) {
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
	// ������� x_{k+1} = x_k - \theta * u_k
	for (int i = 0; i < N_number; ++i) {
		cur_X[i] = cur_X[i] - theta * u_k[i];
	}
	cur_X[i_k - 1] = 0;
	for (int i = 0; i < Nk_indexes.size(); ++i) {
		if (Nk_indexes[i] == i_k) {
			Nk_indexes[i] = j_k;
			break;
		}
	}
}

comb_t GetIndexesPositiveUk(comb_t& Nk_plus_indexes, column_t& u_k) {
	comb_t u_k_positive_indexes;
	for (auto temp : Nk_plus_indexes) {
		if (u_k[temp - 1] > 0)
			u_k_positive_indexes.push_back(temp);
	}
	return u_k_positive_indexes;
}

void GetLk(comb_t& Lk_indexes, comb_t& Nk_plus_indexes, int N) {
	Lk_indexes.clear();
	for (int i = 1; i <= N; ++i)
		if (!IsNumberInCombination(i, Nk_plus_indexes))
			Lk_indexes.push_back(i);
}

comb_t GetDifferenceSets(comb_t& A, comb_t& B) {
	comb_t res;
	for (auto temp : A)
		if (!IsNumberInCombination(temp, B))
			res.push_back(temp);
	return res;
}

// Nk_indexes - ������� �������� ��������
SimplexState IterSimplex(matrix_t& main_A, int M_number, comb_t& Nk_indexes, column_t& main_C, column_t& cur_X, column_t& cur_Y) {
	// ������������� �������� ������� ��� ���������� �������� �������
	matrix_t A_M_Nk(M_number);
	for (int i = 0; i < M_number; ++i)
		A_M_Nk[i].resize(M_number);

	int N_number = main_C.size();

	// ��������� ����� �������� N+ 
	vector<int> Nk_plus_indexes;
	for (int i = 0; i < N_number; ++i)
		if (cur_X[i] > 0)
			Nk_plus_indexes.push_back(i + 1);

	bool is_cur_X_singular = true;

	// �� ���� ���� ������� ������, det(A[M][Nk]) != 0

	// ���� ������� ������ ���������� (|N+| = |M|), �� ��������� ���������� ������� �� �������� �������� �� N+
	if (M_number == Nk_plus_indexes.size()) {
		GetCurMatrix(main_A, A_M_Nk, Nk_plus_indexes);
		double det = Determinant(A_M_Nk);
		cout << det << endl;
		if (det != 0)
			is_cur_X_singular = false;
	}

	if (is_cur_X_singular) {
		/* ���� �� ��������, �� ����� � ������� �������� N+ ��������� ������� �� N0, ����� ������� ������� ���� ����������� */

		// �������������� ������ ���������
		/*for (int i = 0; i < M_number; ++i)
			Nk_indexes[i] = i + 1;*/

		bool trigger = false;

		while (true) {
			// ���������: ��� �� ������� �� A[][N+] �� �����
			for (int i = 0; i < Nk_plus_indexes.size(); ++i) {
				if (!IsNumberInCombination(Nk_plus_indexes[i], Nk_indexes)) {
					// �� ����� �����-�� ������� - ������� ��������� ���������, ������� ������� � ������ �� �������� ��������
					trigger = NextSet(Nk_indexes, N_number, M_number);
					break;
				}
			}
			// ������� �������� - ��������� ��� � ���� ������ ��������� ������ ��� ���������� ���������
			if (trigger == true) {
				trigger = false;
				continue;
			}
			// ���� ��� �� ��� ������� �� �����, �� ��������� ���������� ������� �� ��������
			GetCurMatrix(main_A, A_M_Nk, Nk_indexes);
			// � ��������� �� ������������
			if (Determinant(A_M_Nk) != 0) {
				break;
			}
			else { // ���� ����������� - ������� ��������� ��������� � ���� ������ ��������� ������ ��� ���������� ���������
				NextSet(Nk_indexes, N_number, M_number);
				continue;
			}
		}
	}

	// ������� ������� ������ C[N_k]
	column_t C_Nk(M_number);
	int i = 0;
	for (auto temp : Nk_indexes) {
		C_Nk[i] = main_C[temp - 1];
		++i;
	}

	/* ������ y^T [M] * A[M][Nk] = c^T[M] */
	/* ��� ����� ������������� ������� 
	   A[Nk][M] * y[M] = c[M]	*/
	TransposeQuadMatrix(A_M_Nk);
	// � ����� ���� ������� ��������
	cur_Y = RotationMethod(A_M_Nk, C_Nk);

	// ������ �������� � A[M][Nk] ������� B[Nk][M]
	matrix_t B_Nk_M = InverseMatrix(A_M_Nk);

	// ������ ������� � ��������� ���������
	TransposeQuadMatrix(A_M_Nk);

	// ���������� ������ �������� L_k 
	comb_t L_k;
	GetLk(L_k, Nk_indexes, N_number);

	// ������� cur_X[L_k] = 0 (������������� ����)

	// ������� d_k[L_k]
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

	// ���� d_k �� ������������, �� ��������� �. � �. ������� �������������, ����� �������� ����������� ��������
	if (NonNegativityOfVector(d_k))
		return SimplexState::OPTIMAL;

	/* ���� ���, �� �������� ������� u_k[N] */
	// ������ j_k
	int j_k = -1;
	int iter = 0;
	for (auto temp : L_k) {
		if (d_k[iter] < 0) {
			j_k = temp;
			break;
		}
		++iter;
	}

	column_t u_k(N_number);
	// u_k[Nk] = B[Nk][M] * A[M][j_k]
	int j = 0;
	for (auto temp : Nk_indexes) {
		for (int i = 0; i < M_number; ++i) {
			u_k[temp - 1] += B_Nk_M[j][i] * main_A[i][j_k - 1];
		}
		++j;
	}
	j = 0;
	u_k[j_k - 1] = -1; // u_k[j_k] = -1
	// u_k[L_k \ j_k] = 0 (��� �������� ��� �������������)

	// ���� u_k �� ������������, �� ������� �� ����������
	if (NonPositivityOfVector(u_k))
		return SimplexState::UNLIMITED;

	bool is_cur_X_supported = true;

	/* ���� x_k[N] �������������, �� ���������� \theta, ������� x_{k+1} � �������� ������ �������� */
	if (is_cur_X_singular == false) {
		CalculateNextSupportVector(Nk_indexes, Nk_plus_indexes, u_k, cur_X, N_number, j_k);
		GetCurMatrix(main_A, A_M_Nk, Nk_indexes);
		if (Determinant(A_M_Nk) == 0)
			is_cur_X_singular = true;
		else
			return SimplexState::NEXT;
	}

	if (is_cur_X_singular == true) {
		/* ���� �����������, �� ������������� ��� ������ */
		// ���� u_k[Nk \ Nk+] <= 0,�� ������� \theta
		comb_t u_k_positive_indexes = GetIndexesPositiveUk(Nk_plus_indexes, u_k);

		if (u_k_positive_indexes.size() == 0) {
			CalculateNextSupportVector(Nk_indexes, Nk_plus_indexes, u_k, cur_X, N_number, j_k);
			GetCurMatrix(main_A, A_M_Nk, Nk_indexes);
			if (Determinant(A_M_Nk) == 0)
				is_cur_X_supported = false;
			else 
				return SimplexState::NEXT;
		}

		if ((u_k_positive_indexes.size() != 0) || is_cur_X_supported == false) {
			// ����� ������ �����
			matrix_t current_main_matrix;
			matrix_t current_sub_matrix;
			comb_t Nk_minus_Nk_plus = GetDifferenceSets(Nk_indexes, Nk_plus_indexes);
			for (auto temp_i : Nk_minus_Nk_plus) {
				for (auto temp_j : L_k) {
					// ������ �������
					current_main_matrix = main_A;
					current_sub_matrix = A_M_Nk;
					for (int i = 0; i < M_number; ++i)
						current_main_matrix[i][temp_i - 1] = main_A[i][temp_j - 1];
					GetCurMatrix(current_main_matrix, current_sub_matrix, Nk_indexes);
					// �������, ����������� ��� ��� �������
					if (Determinant(current_sub_matrix) != 0) {
						// ���� ��, �� ������ �����
						int index;
						for (index = 0; index < Nk_indexes.size(); ++index)
							if (Nk_indexes[index] == temp_i)
								break;

						Nk_indexes[index] = temp_j;
						return SimplexState::NEXT;
					}
					else // ���� ���, �� �������� ��� ������ (�������� ������������)
						continue;
				}
			}
			return SimplexState::NEXT;
		}
		else
			return SimplexState::NEXT;
	}
}

tuple<double, column_t, column_t, comb_t> SimplexMethod(canon_problem_t& problem, column_t& cur_X, comb_t& basis) {
	// ���������� �������
	matrix_t matrix_A0;
	column_t column_B;
	column_t column_C;
	tie(matrix_A0, column_B, column_C) = problem;

	// ���������������� ������� ��� ��������
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

	// ����� ������� ���������� ����� �������
	// ��������: �������� �� ������� �������������
	/*if (RankMatrix(matrix_A) != M) {
		cout << "������� ����������� �� �������� �������������." << endl << "��������� ���������� ������." << endl;
		exit(-1);
	}*/

	column_t cur_Y(M);

	SimplexState state = SimplexState::NEXT;
	while (true) {
		// ��������� � ��������� �������� (����� �������� ������� ��� �������� ������)
		if (state == SimplexState::NEXT)
			state = IterSimplex(matrix_A, M, basis, column_C, cur_X, cur_Y);

		// ���� ������� �� ����������, �� ������� -"�������������"
		else if (state == SimplexState::UNLIMITED)
			return make_tuple(-DBL_MAX, cur_X, cur_Y, basis);

		// ���� ����� ����������� �������, �� ������� �������� ������� ����
		else if (state == SimplexState::OPTIMAL) {
			double result = 0;
			for (int i = 0; i < N; ++i)
				result += column_C[i] * cur_X[i];
			return make_tuple(result, cur_X, cur_Y, basis);
		}
	}
}

tuple<column_t, comb_t> GetInitialApprox(canon_problem_t& problem) {
	// ���������� �������
	matrix_t matrix_A;
	column_t column_B;
	column_t column_C;
	tie(matrix_A, column_B, column_C) = problem;

	/* ����������� �������� �������� */
	int b_size = column_B.size();
	int c_size = column_C.size();
	column_C.resize(b_size + c_size);
	for (int i = 0; i < c_size; ++i)
		column_C[i] = 0;
	for (int i = c_size; i < b_size + c_size; ++i)
		column_C[i] = 1;

	/* ��������� ������� A */
	matrix_A.resize(matrix_A.size() + b_size);
	for (int i = matrix_A.size() - b_size; i < matrix_A.size(); ++i) {
		matrix_A[i].resize(b_size);
		matrix_A[i][i - (matrix_A.size() - b_size)] = 1;
	}
		
	/* ����� �������������� ������ � ��� b[M] >= 0 (��. ConvertGeneralToCanon) */

	/* ������ ��������� ����������� */
	column_t first_approx(column_C.size());
	for (int i = 0; i < c_size; ++i)
		first_approx[i] = 0;
	for (int i = c_size; i < b_size + c_size; ++i)
		first_approx[i] = column_B[i - c_size];

	/* ������ ������������� ����� (��������� �������� �������� ��������� �������) */
	comb_t basis;
	for (int i = c_size; i < c_size + b_size; ++i)
		basis.push_back(i + 1);

	canon_problem_t help_problem = make_tuple(matrix_A, column_B, column_C);

	cout << "������ ������ ���������� ����������� ������� �������������� ������:" << endl;

	SimplexMethod(help_problem, first_approx, basis);

	cout << "����� ������ ���������� ����������� ������� �������������� ������." << endl;

	/* ���� ���� �� ���� ���������� y[m] ������������, �� �������� ������ �� ����� �� ����� ���������� ����� */
	for (int i = c_size; i < c_size + b_size; ++i) {
		if (first_approx[i] > 0) {
			cout << "�������� ������ �� ����� �� ����� ���������� �����." << endl << "��������� ���������� ������." << endl;
			exit(-1);
		}
	}
	/* ��������� �� ������� x[N] - ������� ������ ��� �������� ������ */
	/*for (int i = 0; i < b_size; ++i) {
		first_approx[i] = first_approx[i + c_size];
	}*/
	first_approx.resize(c_size);

	cout << "��������� ����������� ������� �������������� ������ ���� �������." << endl;

	return make_tuple(first_approx, basis);
}

tuple<double, column_t, column_t, matrix_t> SolveProblemWithSimplexMethod(canon_problem_t& problem) {
	column_t X;
	comb_t basis;
	tie(X, basis) = GetInitialApprox(problem);

	double opt_value; 
	column_t Y;
	tie(opt_value, X, Y, basis) = SimplexMethod(problem, X, basis);

	/*matrix_t A;
	column_t B, C;
	tie(A, B, C) = problem;
	matrix_t res(B.size());
	for (int i = 0; i < B.size(); ++i)
		res[i].resize(basis.size());
	int i = 0;
	for (auto temp : basis) {
		for (int j = 0; j < B.size(); ++j) {
			res[j][i] = A[j][temp - 1];
		}
		++i;
	}

	basis.clear();
	B.clear();
	C.clear();
	for (auto& temp : A)
		temp.clear();
	A.clear();*/

	matrix_t res;

	return make_tuple(opt_value, X, Y, res);
}
