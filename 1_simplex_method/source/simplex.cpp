#include "simplex.h"
using namespace std;

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
	sort(Nk_indexes.begin(), Nk_indexes.end());
	// ��������� Nk_plus_indexes
	Nk_plus_indexes.clear();
	for (int i = 0; i < cur_X.size(); ++i) {
		if (cur_X[i] > 0)
			Nk_plus_indexes.push_back(i + 1);
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

void GetNkPlus(comb_t& Nk_plus_indexes, column_t& cur_X) {
	Nk_plus_indexes.clear();
	for (int i = 0; i < cur_X.size(); ++i)
		if (cur_X[i] > 0)
			Nk_plus_indexes.push_back(i + 1);
}

void NormalizeVector(column_t& X) {
	for (int i = 0; i < X.size(); ++i) {
		if (abs(X[i]) < 1e-10)
			X[i] = 0;
	}
}

// Nk_indexes - ������� �������� ��������
SimplexState IterSimplex(matrix_t& main_A, comb_t& Nk_indexes, column_t& main_C, column_t& cur_X, column_t& cur_Y) {
	
	int M_number = main_A.size();
	int N_number = main_C.size();

	// ������������� �������� ������� ��� ���������� �������� �������
	matrix_t A_M_Nk(M_number);
	for (int i = 0; i < M_number; ++i)
		A_M_Nk[i].resize(M_number);

	// ��������� ����� �������� N+ 
	vector<int> Nk_plus_indexes;
	GetNkPlus(Nk_plus_indexes, cur_X);

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
			NextSet(Nk_indexes, N_number, M_number);
			// ���������: ��� �� ������� �� A[][N+] �� �����
			for (int i = 0; i < Nk_plus_indexes.size(); ++i) {
				if (!IsNumberInCombination(Nk_plus_indexes[i], Nk_indexes)) {
					// �� ����� �����-�� ������� - ������� ��������� ���������, ������� ������� � ������ �� �������� ��������
					trigger = NextSet(Nk_indexes, N_number, M_number);
					if (trigger == false) {
						Nk_indexes.clear();
						for (int i = 1; i < M_number; ++i)
							Nk_indexes.push_back(i);
						trigger = true;
					}
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
	NormalizeVector(cur_Y);

	// ������ ������� � ��������� ���������
	TransposeQuadMatrix(A_M_Nk);

	// ������ �������� � A[M][Nk] ������� B[Nk][M]
	// matrix_t B_Nk_M = InverseMatrix(A_M_Nk);

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
	NormalizeVector(d_k);

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

	// u_k[Nk] = B[Nk][M] * A[M][j_k]
	column_t cur_u_k(M_number);
	column_t A_M_jk;// (M_number);
	for (int i = 0; i < M_number; ++i)
		A_M_jk.push_back(main_A[i][j_k - 1]);
	cur_u_k = RotationMethod(A_M_Nk, A_M_jk);
	NormalizeVector(cur_u_k);
	column_t u_k(N_number);
	i = 0;
	for (auto temp : Nk_indexes) {
		u_k[temp - 1] = cur_u_k[i];
		++i;
	}
	u_k[j_k - 1] = -1; // u_k[j_k] = -1
	// u_k[L_k \ j_k] = 0 (��� �������� ��� �������������)

	// ���� u_k �� ������������, �� ������� �� ����������
	if (NonPositivityOfVector(u_k))
		return SimplexState::UNLIMITED;

	bool is_cur_X_supported = true;

	/* ���� x_k[N] �������������, �� ���������� \theta, ������� x_{k+1} � �������� ������ �������� */
	if (is_cur_X_singular == false) {
		CalculateNextSupportVector(Nk_indexes, Nk_plus_indexes, u_k, cur_X, N_number, j_k);
		NormalizeVector(cur_X);
		GetCurMatrix(main_A, A_M_Nk, Nk_indexes);
		if (Determinant(A_M_Nk) == 0)
			is_cur_X_singular = true;
		else
			return SimplexState::NEXT;
	}

	if (is_cur_X_singular == true) {
		/* ���� �����������, �� ������������� ��� ������ */
		// ���� u_k[Nk \ Nk+] <= 0,�� ������� \theta
		comb_t diff_Nk_Nk_plus = GetDifferenceSets(Nk_indexes, Nk_plus_indexes);
		comb_t u_k_positive_indexes = GetIndexesPositiveUk(diff_Nk_Nk_plus, u_k);

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
			GetNkPlus(Nk_plus_indexes, cur_X);
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
						sort(Nk_indexes.begin(), Nk_indexes.end());
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
			state = IterSimplex(matrix_A, basis, column_C, cur_X, cur_Y);

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

tuple<double, column_t, column_t, comb_t> SolveProblemWithSimplexMethod(canon_problem_t& problem) {
	column_t X;
	comb_t basis;
	tie(X, basis) = GetInitialApprox(problem);

	double opt_value;
	column_t Y;
	tie(opt_value, X, Y, basis) = SimplexMethod(problem, X, basis);

	matrix_t A;
	column_t B, C;
	tie(A, B, C) = problem;
	/*
	matrix_t res(B.size());
	for (int i = 0; i < B.size(); ++i)
		res[i].resize(basis.size());
	int i = 0;
	for (auto temp : basis) {
		for (int j = 0; j < B.size(); ++j) {
			res[j][i] = A[temp - 1][j];
		}
		++i;
	}
	*/

	B.clear();
	C.clear();
	for (auto& temp : A)
		temp.clear();
	A.clear();

	return make_tuple(opt_value, X, Y, basis);
}
