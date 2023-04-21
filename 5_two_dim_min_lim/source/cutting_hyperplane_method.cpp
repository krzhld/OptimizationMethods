#include "cutting_hyperplane_method.h"


solving_linear_problem_t SolvingLinearProblemInFirstIter(Task t, polyhedron_t Sk)
{
	//кириешка - ваш выход
}

solving_linear_problem_t SolvingLinearProblem(Task t, polyhedron_t Sk, column_t yk)
{
	//кириешка - ваш выход
}

hyperplane_t GetCuttingHyperplane(column_t xk, Task t)
{
	column_t subgrk = t.SubgradientLim(xk); // ищем субградиент функции, задающей ограничение 

	double bk = -t.Limit(xk) + MultiplipliedVectors(subgrk, xk); 

	hyperplane_t result = make_tuple(subgrk, bk);

	return result;
}

polyhedron_t AddCuttingHiperplaneInPolyhedron(polyhedron_t Sk, hyperplane_t h)
{
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


solving_t CuttingHyperplaneMethod(column_t x0, polyhedron_t S0, Task t, double eps)
{
	column_t xk, yk;
	polyhedron_t Sk = S0;
	tie(xk, yk) = SolvingLinearProblemInFirstIter(t, Sk); //ищем решение прямой и двойственной задачи на первой итерации без использования решения двойственной задачи с прошлой итерации
	
	while (Norm(DiffVector(xk, x0)) > eps)
	{
		x0 = xk;
		hyperplane_t cuttHyperplane = GetCuttingHyperplane(xk, t); //ищем отсекающую гиперплоскость 

		Sk = AddCuttingHiperplaneInPolyhedron(Sk, cuttHyperplane); //добавляем отсекающую гиперплоскость в множество ограничений задачи ЛП

		solving_linear_problem_t solvingLinearProblem = SolvingLinearProblem(t, Sk, yk); //ищем решение прямой и двойственной задачи ЛП

		column_t xk, yk;
		tie(xk, yk) = solvingLinearProblem; 

	}

	column_t result_column = xk;
	double result = t.MinFunc(result_column);

	return make_tuple(result_column, result);
}
