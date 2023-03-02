#include "dual_problem.h"

using namespace std;

/*���������������� � ��������� �� -1 ������� ��� �������� ������������ ������ � ������������*/
matrix_t TransformationMatrix(matrix_t matrix)
{
	column_t columnTransform;
	matrix_t transform;
	int n = size(matrix);
	int m = size(matrix[1]);
	columnTransform.resize(n);
	transform.resize(m);
	for (int k = 0; k < m;k++)
	{
		transform[k] = columnTransform;
	}
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			transform[j][i] = -matrix[i][j];
		}
	}
	return transform;
}

/*��������� ������������ ������ */
general_problem_t getDualLinearProblem(general_problem_t &problem) 
{
	matrix_t A; 
	column_t b, c;
	int M1, N1;
	tie(A, b, c, M1, N1) = problem;
	// ������������ ������������ ������� � ������������ �������� 
	matrix_t dualA; 
	column_t dualb, dualc;
	
	//��������� ������� A ������������ ������ 
	dualA = TransformationMatrix(A);

	//������� b ������������ ������  - ��� ������� c ������ ������, ��� ��� �������� ������� �� ������ - 
	/*��� �������� ������� ����� �������� �� - 1, ����� ������ ���������� ������������ ������ � ����� ������ ��
	(�������� ������ ���� �����������(>= ))*/
	dualb.resize(size(c));
	for (int i = 0; i < size(c); i++)
	{
		dualb[i] = -c[i];
	}

	//������� c ������������ ������ - ��� ������� b ������ ������
	dualc.resize(size(b));
	for (int j = 0; j < size(b); j++)
	{
		dualc[j] = b[j];
	}

	/*����� N1 � M1 - ��� ���������� ���������� � ������������� �� ����
	� ���������� ���������� � ������� ����������� ��������������*/
	int dualN1 = M1; 
	int dualM1 = N1;

	general_problem_t dualProblem = make_tuple(dualA, dualb, dualc, dualM1, dualN1);

	return dualProblem;
}

