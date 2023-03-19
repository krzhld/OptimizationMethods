#define _CRT_SECURE_NO_WARNINGS
#include "transport_problem.h"
#include "simplex.h"
using namespace std;

int main()
{

	setlocale(LC_ALL, "Russian");

	transport_problem_t problem = ReadFromFileTransportProblem("task.txt");

	matrix_t c;
	column_t b, a;
	tie(a, b, c) = problem;


	printf("Матрица тарифов:\n");
	for (int i = 0; i < size(a); i++)
	{
		for (int j = 0; j < size(b); j++)
		{
			printf("%lf ", c[i][j]);
		}
		printf("\n");
	}
	printf("\nСтолбец поставщиков:\n");
	for (int j = 0; j < size(a); j++)
	{
		printf("%lf\n", a[j]);

	}
	printf("\nСтолбец покупателей:\n");
	for (int j = 0; j < size(b); j++)
	{
		printf("%lf\n", b[j]);

	}
	printf("\n");

	printf("\nРешение задачи методом потенциалов\n");

	solving_t solving = SolveTransportProblem(problem);
	matrix_t X;
	double result;
	tie(X, result) = solving;

	printf("\nНайденное оптимальное решение методом потенциалов:\n");
	for (int i = 0; i < size(a); i++)
	{
		for (int j = 0; j < size(b); j++)
		{
			printf("%lf ", X[i][j]);
		}
		printf("\n");
	}
	printf("\n");
	printf("Минимальные затраты на перевозку, найденные методом потенциалов:\n");
	printf("%lf\n", result);

	// -----------------------------------------------------------------------------------------------------------------------------------------
	
	printf("\nРешение задачи симплекс методом\n");
	canon_problem_t canon_problem = GetCanonProblemFromTransportProblem(problem);
	column_t X_col, Y;
	double opt_value;
	comb_t basis;
	tie(opt_value, X_col, Y, basis) = SolveProblemWithSimplexMethod(canon_problem);

	printf("\nНайденное оптимальное решение симплекс методом:\n");
	for (int i = 0; i < size(a); i++)
	{
		for (int j = 0; j < size(b); j++)
		{
			printf("%lf ", X_col[i * size(b) + j]);
		}
		printf("\n");
	}
	printf("\n");
	printf("Минимальные затраты на перевозку, найденные симплекс методом:\n");
	printf("%lf\n", opt_value);

	// -----------------------------------------------------------------------------------------------------------------------------------------
	
	printf("\n\nРешение задачи в случае недопоставки с учетом штрафов\n");
	transport_problem_t problem_with_imbalance = ReadFromFileTransportProblem("task_with_imbalance.txt");

	matrix_t c_ib;
	column_t b_ib, a_ib;
	tie(a_ib, b_ib, c_ib) = problem_with_imbalance;

	printf("Матрица тарифов:\n");
	for (int i = 0; i < size(a); i++)
	{
		for (int j = 0; j < size(b); j++)
		{
			printf("%lf ", c_ib[i][j]);
		}
		printf("\n");
	}
	printf("\nСтолбец поставщиков:\n");
	for (int j = 0; j < size(a); j++)
	{
		printf("%lf\n", a_ib[j]);

	}
	printf("\nСтолбец покупателей:\n");
	for (int j = 0; j < size(b); j++)
	{
		printf("%lf\n", b_ib[j]);

	}
	printf("\n");

	solving_t solving_ib = SolveTransportProblem(problem_with_imbalance);
	matrix_t X_ib;
	double result_ib;
	tie(X_ib, result_ib) = solving_ib;

	printf("Найденное оптимальное решение с учетом дизбаланса:\n");
	for (int i = 0; i < size(a_ib); i++)
	{
		for (int j = 0; j < size(b_ib); j++)
		{
			printf("%lf ", X_ib[i][j]);
		}
		printf("\n");
	}
	printf("\n");
	printf("Минимальные затраты на перевозку с учетом штрафов:\n");
	printf("%lf\n", result_ib);
	return 0;

}
