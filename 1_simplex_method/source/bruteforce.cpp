#include "bruteforce.h"
using namespace std;

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
