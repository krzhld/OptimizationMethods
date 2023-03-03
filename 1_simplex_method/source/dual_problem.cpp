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

/*�������� ���������� ����� � ������ �� �������� � ���������*/
bool IsNumberInBasis(int number, comb_t comb, int combSize)
{
	for (int i = 0; i < combSize; i++)
	{
		if (number == comb[i])
			return true;
	}
	return false;
}

/*����� �������� ������� ������� ������*/
matrix_t InverseMatrix(matrix_t A)
{
	int n = size(A);

	column_t Ecolumn; // ������ ��������� �������
	Ecolumn.resize(n);

	matrix_t E; // ������������� ��������� �������
	E.resize(n);
	for (int i = 0; i < n; i++)
	{
		E[i] = Ecolumn;
	}

	double temp;
	for (int i = 0; i < n; i++)
	{
		E[i][i] = 1;
	}

	for (int k = 0; k < n; k++)
	{
		temp = A[k][k];

		for (int j = 0; j < n; j++)
		{
			A[k][j] /= temp;
			E[k][j] /= temp;
		}

		for (int i = k + 1; i < n; i++)
		{
			temp = A[i][k];

			for (int j = 0; j < n; j++)
			{
				A[i][j] -= A[k][j] * temp;
				E[i][j] -= E[k][j] * temp;
			}
		}
	}

	for (int k = n - 1; k > 0; k--)
	{
		for (int i = k - 1; i >= 0; i--)
		{
			temp = A[i][k];

			for (int j = 0; j < n; j++)
			{
				A[i][j] -= A[k][j] * temp;
				E[i][j] -= E[k][j] * temp;
			}
		}
	}

	//��������������� ��������� ������� �������� �������� �������� � ������� � 
	return E;
}

column_t MultipliedQuadMatrixAndColumn(matrix_t M, column_t c)
{
	int n = size(c);

	column_t X; ///��������� ������������
	X.resize(n);

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			X[i] += M[i][j] * c[j];
		}
	}

	return X;
}

/*�������������� ������� ������������ ������ �� ������� ������ ������*/
column_t SolvingDualProblem(canon_problem_t problem, column_t X, comb_t optBasis)
{
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
	{
		basis[i] = columnBasis;
	}

	column_t cB; // ������ ��������� ������� �������, ��������������� ������ ������������ �������
	cB.resize(m);

	int numb = 0; //�������

	for (int i = 0; i < n; i++)
	{
		if (IsNumberInBasis(i + 1, optBasis, m))
		{
			for (int j = 0; j < m; j++)
			{
				basis[numb][j] = A[i][j];
			}
			cB[numb] = c[i];
			numb++;
		}
	}

	matrix_t invBasis = InverseMatrix(basis); // ����������� ������� ������

	Y = MultipliedQuadMatrixAndColumn(invBasis, cB); 

	double result = 0;
	for (int i = 0; i < size(Y); i++)
	{
		result += Y[i] * b[i];
	}

	printf("\n %lf \n", result);

	return Y;

}
