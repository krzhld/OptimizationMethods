#define _CRT_SECURE_NO_WARNINGS
#include "transport_problem.h"
#include "simplex.h"
using namespace std;

int main()
{

	setlocale(LC_ALL, "Russian");

	transport_problem_t problem = ReadFromFileTransportProblem("task.txt");

	matrix_t c;
	column_t b, a, f;
	tie(a, b, c, f) = problem;


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

	printf("\nРешение задачи в случае недопоставки с учетом штрафов\n");
	printf("\nПервый случай\n");
	transport_problem_t problem_with_imbalance_1 = ReadFromFileTransportProblem("task_with_imbalance_1.txt");

	matrix_t c_ib1;
	column_t b_ib1, a_ib1, f_ib1;
	tie(a_ib1, b_ib1, c_ib1, f_ib1) = problem_with_imbalance_1;

	double imbalance1 = SumColumn(b_ib1) - SumColumn(a_ib1);

	
	printf("Матрица тарифов:\n");
	for (int i = 0; i < size(a_ib1); i++)
	{
		for (int j = 0; j < size(b_ib1); j++)
		{
			printf("%lf ", c_ib1[i][j]);
		}
		printf("\n");
	}
	printf("\nСтолбец поставщиков:\n");
	for (int j = 0; j < size(a_ib1); j++)
	{
		printf("%lf\n", a_ib1[j]);

	}
	printf("\nСтолбец покупателей:\n");
	for (int j = 0; j < size(b_ib1); j++)
	{
		printf("%lf\n", b_ib1[j]);

	}
	printf("\n");
	printf("\nМера штрафов потребителей:\n");
	for (int j = 0; j < size(f_ib1); j++)
	{
		printf("%lf\n", f_ib1[j]);

	}
	printf("\n");

	printf("Дизбаланс:\n");
	printf("%lf\n", imbalance1);

	solving_t solving_ib1 = SolveTransportProblem(problem_with_imbalance_1);
	matrix_t X_ib1;
	double result_ib1;
	tie(X_ib1, result_ib1) = solving_ib1;

	double fine1 = result_ib1 - FindCost(c_ib1, X_ib1, size(a_ib1), size(b_ib1));

	printf("Найденное оптимальное решение с учетом дизбаланса:\n");
	for (int i = 0; i < size(a_ib1); i++)
	{
		for (int j = 0; j < size(b_ib1); j++)
		{
			printf("%lf ", X_ib1[i][j]);
		}
		printf("\n");
	}
	printf("\n");
	printf("Минимальные затраты на перевозку с учетом штрафов:\n");
	printf("%lf\n", result_ib1);
	printf("Штраф:\n");
	printf("%lf\n", fine1);


	printf("\nВторой случай\n");
	transport_problem_t problem_with_imbalance_2 = ReadFromFileTransportProblem("task_with_imbalance_2.txt");

	matrix_t c_ib2;
	column_t b_ib2, a_ib2, f_ib2;
	tie(a_ib2, b_ib2, c_ib2, f_ib2) = problem_with_imbalance_2;

	double imbalance2 = SumColumn(b_ib2) - SumColumn(a_ib2);


	printf("Матрица тарифов:\n");
	for (int i = 0; i < size(a_ib2); i++)
	{
		for (int j = 0; j < size(b_ib2); j++)
		{
			printf("%lf ", c_ib2[i][j]);
		}
		printf("\n");
	}
	printf("\nСтолбец поставщиков:\n");
	for (int j = 0; j < size(a_ib2); j++)
	{
		printf("%lf\n", a_ib2[j]);

	}
	printf("\nСтолбец покупателей:\n");
	for (int j = 0; j < size(b_ib2); j++)
	{
		printf("%lf\n", b_ib2[j]);

	}
	printf("\n");
	printf("\nМера штрафов потребителей:\n");
	for (int j = 0; j < size(f_ib2); j++)
	{
		printf("%lf\n", f_ib2[j]);

	}
	printf("\n");

	printf("Дизбаланс:\n");
	printf("%lf\n", imbalance2);

	solving_t solving_ib2 = SolveTransportProblem(problem_with_imbalance_2);
	matrix_t X_ib2;
	double result_ib2;
	tie(X_ib2, result_ib2) = solving_ib2;

	double fine2 = result_ib2 - FindCost(c_ib2, X_ib2, size(a_ib2), size(b_ib2));

	printf("Найденное оптимальное решение с учетом дизбаланса:\n");
	for (int i = 0; i < size(a_ib2); i++)
	{
		for (int j = 0; j < size(b_ib2); j++)
		{
			printf("%lf ", X_ib2[i][j]);
		}
		printf("\n");
	}
	printf("\n");
	printf("Минимальные затраты на перевозку с учетом штрафов:\n");
	printf("%lf\n", result_ib2);
	printf("Штраф:\n");
	printf("%lf\n", fine2);
	return 0;

}
