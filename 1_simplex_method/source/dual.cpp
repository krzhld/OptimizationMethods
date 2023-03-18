#include "dual.h"

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
	// инциализация двойственной матрицы и двойственных столбцов 
	matrix_t dualA;
	column_t dualb, dualc;

	// получение матрицы A двойственной задачи 
	dualA = TransformationMatrix(A);

	//столбец b двойственной задачи  - это столбец c прямой задачи, где все значения берутся со знаком - 
	/* все элементы столбца нужно умножить на - 1, чтобы свести полученную двойственную задаче к общей задаче ЛП
	(получить нужный знак неравенства(>= )) */
	dualb.resize(size(c));
	for (int i = 0; i < size(c); i++)
		dualb[i] = -c[i];

	// столбец c двойственной задачи - это столбец b прямой задачи
	dualc.resize(size(b));
	for (int j = 0; j < size(b); j++)
		dualc[j] = b[j];

	/* числа N1 и M1 - это количество переменных с ограничениями на знак
	и количество неравенств в системе ограничений соответственно */
	int dualN1 = M1;
	int dualM1 = N1;

	general_problem_t dualProblem = make_tuple(dualA, dualb, dualc, dualM1, dualN1);

	return dualProblem;
}

column_t SolvingDualProblem(canon_problem_t problem, column_t X, comb_t optBasis) {
	matrix_t A;
	column_t b, c;
	tie(A, b, c) = problem;

	int n = size(c);
	int m = size(b);

	column_t Y; // решение двойственной задачи 

	column_t columnBasis; // столбец базиса
	columnBasis.resize(m);

	matrix_t basis; //базис оптимального решения 
	basis.resize(m);

	for (int i = 0; i < m; i++)
		basis[i] = columnBasis;

	column_t cB; // вектор компонент целевой функцию, соответствующий базису оптимального решения
	cB.resize(m);

	int numb = 0; //счетчик

	for (int i = 0; i < n; i++) {
		if (IsNumberInCombination(i + 1, optBasis)) {
			for (int j = 0; j < m; j++)
				basis[numb][j] = A[i][j];

			cB[numb] = c[i];
			numb++;
		}
	}

	matrix_t invBasis = InverseMatrix(basis); // инвертируем матрицу базиса

	Y = MultipliedMatrixAndColumn(invBasis, cB);

	double result = 0;
	for (int i = 0; i < size(Y); i++)
		result += Y[i] * b[i];

	return Y;
}
