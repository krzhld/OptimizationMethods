#include "simplex.h"
#include "cutting_hyperplane_method.h"


solving_linear_problem_t SolvingLinearProblemInFirstIter(Task& t) {
	polyhedron_t polyhedron = t.GetS0();
	matrix_t A0;
	column_t b;
	tie(A0, b) = polyhedron;
	
	column_t c = { 0, 0, 0, 0, -1, 1 };
	int l = b.size();

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < l; j++) {
			A0[j].resize(3 + 3);
			A0[j][3 + i] = -A0[j][i];
		}
	}

	for (int i = 0; i < l; i++)
		c.push_back(0);


	for (int i = 0; i < l; i++) {
		A0[i].resize(3 + 3 + l);
		A0[i][3 + 3 + i] = 1;
	}
	
	column_t x, y;
	int temp;
	comb_t basis;

	matrix_t A(3 + 3 + l);
	for (int i = 0; i < 3 + 3 + l; ++i) {
		A[i].resize(l);
		for (int j = 0; j < l; ++j) {
			A[i][j] = A0[j][i];
		}
	}
	A0.clear();

	canon_problem_t canon_problem = make_tuple(A, b, c);

	tie(temp, x, y, basis) = SolveProblemWithSimplexMethod(canon_problem);
	
	for (int i = 0; i < 3; i++) {
		x[i] = x[i] - x[i + 3];
	}
	x.resize(3);

	return make_tuple(x, y, basis);
}

solving_linear_problem_t SolvingLinearProblem(Task& t, polyhedron_t& Sk, column_t& yk, comb_t basis) {
	canon_problem_t canon_problem;
	matrix_t A;
	column_t b;

	tie(A, b) = Sk;
	column_t c = { 0, 0, -1 };

	canon_problem = make_tuple(A, c, b);

	column_t x, y;
	double res;
	tie(res, y, x, basis) = SimplexMethod(canon_problem, yk, basis);

	return make_tuple(y, x, basis);
}

hyperplane_t GetCuttingHyperplane(column_t& xk, Task& t) {
	column_t subgrk = t.SubgradientLim(xk); // ищем субградиент функции, задающей ограничение 

	double bk = -t.Limit(xk) + MultiplipliedVectors(subgrk, xk); 

	hyperplane_t result = make_tuple(subgrk, bk);

	return result;
}

polyhedron_t AddCuttingHiperplaneInPolyhedron(polyhedron_t& Sk, hyperplane_t& h) {
	matrix_t Ak;
	column_t bk;
	tie(Ak, bk) = Sk;

	column_t ak;
	double b;
	tie(ak, b) = h;

	Ak.push_back(ak);
	bk.push_back(b);

	polyhedron_t NewSk = make_tuple(Ak, bk);
	return NewSk;
}


column_t CuttingHyperplaneMethod(Task& t, double eps) {
	column_t xprev, xk, yk;
	comb_t basis, _;
	polyhedron_t Sk = t.GetS0(); //первоначальное множество, содержащее исходное множество точек
	tie(xk, yk, _) = SolvingLinearProblemInFirstIter(t); //ищем решение прямой и двойственной задачи на первой итерации без использования решения двойственной задачи с прошлой итерации
	_.clear();

	// Формируем базис для двойственной задачи
	basis.push_back(1);
	basis.push_back(3);
	basis.push_back(6);

	do {
		xprev = xk;
		hyperplane_t cuttHyperplane = GetCuttingHyperplane(xk, t); //ищем отсекающую гиперплоскость 

		Sk = AddCuttingHiperplaneInPolyhedron(Sk, cuttHyperplane); //добавляем отсекающую гиперплоскость в множество ограничений задачи ЛП

		yk.push_back(0);
		solving_linear_problem_t solvingLinearProblem = SolvingLinearProblem(t, Sk, yk, basis); //ищем решение прямой и двойственной задачи ЛП

		tie(yk, xk, basis) = solvingLinearProblem; 

	} while (Norm(DiffVector(xk, xprev)) > eps);

	return xk;
}
