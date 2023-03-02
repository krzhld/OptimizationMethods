#include "dual_problem.h"

using namespace std;

/*“ранспонирвоание и умножение на -1 матрицы дл€ сведени€ канонической задачи к двойственной*/
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

/*ѕолучение двойственной задачи */
general_problem_t getDualLinearProblem(general_problem_t &problem) 
{
	matrix_t A; 
	column_t b, c;
	int M1, N1;
	tie(A, b, c, M1, N1) = problem;
	// инциализаци€ двойственной матрицы и двойственных столбцов 
	matrix_t dualA; 
	column_t dualb, dualc;
	
	//получение матрицы A двойственной задачи 
	dualA = TransformationMatrix(A);

	//столбец b двойственной задачи  - это столбец c пр€мой задачи, где все значени€ берутс€ со знаком - 
	/*все элементы столбца нужно умножить на - 1, чтобы свести полученную двойственную задаче к общей задаче Ћѕ
	(получить нужный знак неравенства(>= ))*/
	dualb.resize(size(c));
	for (int i = 0; i < size(c); i++)
	{
		dualb[i] = -c[i];
	}

	//столбец c двойственной задачи - это столбец b пр€мой задачи
	dualc.resize(size(b));
	for (int j = 0; j < size(b); j++)
	{
		dualc[j] = b[j];
	}

	/*числа N1 и M1 - это количество переменных с ограничени€ми на знак
	и количество неравенств в системе ограничений соответственно*/
	int dualN1 = M1; 
	int dualM1 = N1;

	general_problem_t dualProblem = make_tuple(dualA, dualb, dualc, dualM1, dualN1);

	return dualProblem;
}

