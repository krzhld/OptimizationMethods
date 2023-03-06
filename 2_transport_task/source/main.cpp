#define _CRT_SECURE_NO_WARNINGS
#include "transport_problem.h"

int main()
{
	transport_problem_t problem = ReadFromFileTransportProblem("task.txt");

	matrix_t c;
	column_t b, a;
	tie(a, b, c) = problem;

	for (int i = 0; i < size(a); i++)
	{
		for (int j = 0; j < size(b); j++)
		{
			printf("%lf ", c[i][j]);
		}
		printf("\n");
	}
	printf("\n");
	for (int j = 0; j < size(a); j++)
	{
		printf("%lf\n", a[j]);

	}
	printf("\n");
	for (int j = 0; j < size(b); j++)
	{
		printf("%lf\n", b[j]);

	}
	printf("\n");

	solving_t solving = MethodOfPotentials(problem);
	matrix_t X;
	double result;
	tie(X, result) = solving;

	for (int i = 0; i < size(a); i++)
	{
		for (int j = 0; j < size(b); j++)
		{
			printf("%lf ", X[i][j]);
		}
		printf("\n");
	}
	printf("\n");

	printf("\n%lf\n", result);
	return 0;

}