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
	
	cout << endl << "��������:" << endl;
	cout << "������� ������: ";
	for (auto temp : Nk_indexes)
		cout << temp << " ";
	cout << endl << "������� x: ";
	for (auto temp : cur_X)
		cout << temp << " ";
	cout << endl << "������� y: ";
	for (auto temp : cur_Y)
		cout << temp << " ";
	cout << endl;

	int M_number = main_A.size();
	int N_number = main_C.size();

	// ������������� �������� ������� ��� ���������� �������� �������
	matrix_t A_M_Nk(M_number);
	for (int i = 0; i < M_number; ++i)
		A_M_Nk[i].resize(M_number);

	// ��������� ����� �������� N+ 
	vector<int> Nk_plus_indexes;
	GetNkPlus(Nk_plus_indexes, cur_X);

	GetCurMatrix(main_A, A_M_Nk, Nk_indexes);

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

	// ���������� ������ �������� L_k 
	comb_t L_k;
	GetLk(L_k, Nk_indexes, N_number);

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
	if (NonNegativityOfVector(d_k)) {
		FreeMatrix(A_M_Nk);
		return SimplexState::OPTIMAL;
	}	

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
	column_t A_M_jk;
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
	if (Nk_indexes.size() == Nk_plus_indexes.size()) {
		CalculateNextSupportVector(Nk_indexes, Nk_plus_indexes, u_k, cur_X, N_number, j_k);
		NormalizeVector(cur_X);
		GetCurMatrix(main_A, A_M_Nk, Nk_indexes);
		if (Determinant(A_M_Nk) != 0)
			return SimplexState::NEXT;
	}
	else {
		/* ���� �����������, �� ������������� ��� ������ */
		// ���� u_k[Nk \ Nk+] <= 0,�� ������� \theta
		comb_t Nk_minus_Nk_plus = GetDifferenceSets(Nk_indexes, Nk_plus_indexes);
		comb_t u_k_positive_indexes = GetIndexesPositiveUk(Nk_minus_Nk_plus, u_k);

		if (u_k_positive_indexes.size() == 0) {
			CalculateNextSupportVector(Nk_indexes, Nk_plus_indexes, u_k, cur_X, N_number, j_k);
			GetCurMatrix(main_A, A_M_Nk, Nk_indexes);
			if (Determinant(A_M_Nk) != 0)
				return SimplexState::NEXT;		
		}
		else {
			// ����� ������ �����
			int cur_m = Nk_minus_Nk_plus.size();
			comb_t cur_diff = Nk_minus_Nk_plus;
			comb_t cur_Nk_indexes = Nk_plus_indexes;

			bool flag = false;
			NextSetCycle(cur_diff, N_number, cur_m);
			while (true) {
				flag = false;
				// ���������: ������������ �� cur_diff � Nk+
				for (auto temp : cur_diff) {
					if (IsNumberInCombination(temp, Nk_plus_indexes)) {
						NextSetCycle(cur_diff, N_number, cur_m);
						flag = true;
						break;
					}
				}
				if (flag == false) {
					for (auto temp : cur_diff)
						cur_Nk_indexes.push_back(temp);
					sort(cur_Nk_indexes.begin(), cur_Nk_indexes.end());
					GetCurMatrix(main_A, A_M_Nk, cur_Nk_indexes);
					if (Determinant(A_M_Nk) != 0) {
						Nk_indexes = cur_Nk_indexes;
						return SimplexState::NEXT;
					}
					else {
						cur_Nk_indexes.clear();
						cur_Nk_indexes = Nk_plus_indexes;
						NextSetCycle(cur_diff, N_number, cur_m);
					}
				}
			}
		}
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

	// ��������: �������� �� ������� �������������
	if (!IsFullRankMatrix(matrix_A)) {
		cout << "������� ����������� �� �������� �������������." << endl << "��������� ���������� ������." << endl;
		exit(-1);
	}

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

	cout << "������ ������ ���������� ����������� ������� �������������� ������." << endl;

	SimplexMethod(help_problem, first_approx, basis);

	cout << endl << "����� ������ ���������� ����������� ������� �������������� ������." << endl;

	/* ���� ���� �� ���� ���������� y[m] ������������, �� �������� ������ �� ����� �� ����� ���������� ����� */
	for (int i = c_size; i < c_size + b_size; ++i) {
		if (first_approx[i] > 0) {
			cout << "�������� ������ �� ����� �� ����� ���������� �����." << endl << "��������� ���������� ������." << endl;
			exit(-1);
		}
	}
	/* ��������� �� ������� x[N] - ������� ������ ��� �������� ������ */
	first_approx.resize(c_size);

	// ��������� ����������� ������� �������������� ������ ���� �������

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

	B.clear();
	C.clear();
	for (auto& temp : A)
		temp.clear();
	A.clear();

	return make_tuple(opt_value, X, Y, basis);
}
