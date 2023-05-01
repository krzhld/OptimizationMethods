#define _CRT_SECURE_NO_WARNINGS
#include "simplex.h"
#include "dual.h"
#include <iomanip>

using namespace std;

string fix(float x, int p) {
	ostringstream strout;
	strout << fixed << setprecision(p) << x;
	string str = strout.str();
	size_t end = str.find_last_not_of('0') + 1;
	return str.erase(end);
}

int main(void) {
	setlocale(LC_ALL, "Russian");

	general_problem_t problem = FromFileConvertToGeneral("general_task_main.txt");

	/*general_problem_t dual_problem = GetDualLinearProblem(problem);
	// поскольку симплекс метод ищет минимум, а дв задача на максимум, то умножим на -1 критерий качества
	matrix_t dual_A; column_t dual_b, dual_c; int M1, N1;
	tie(dual_A, dual_b, dual_c, M1, N1) = dual_problem;
	for (int i = 0; i < dual_c.size(); ++i)
		dual_c[i] *= -1;
	dual_problem = make_tuple(dual_A, dual_b, dual_c, M1, N1);
	canon_problem_t canon_problem = ConvertGeneralToCanon(dual_problem);

	column_t X, Y;
	double opt_value;
	comb_t basis;
	tie(opt_value, X, Y, basis) = SolveProblemWithSimplexMethod(canon_problem);

	cout << endl << "Оптимальное значение: " << -opt_value << endl << endl;

	cout << "Вектор, сообщающий оптимальное решение прямой задаче: " << endl;
	for (auto temp : X)
		cout << temp << "\t";
	cout << endl << endl;

	cout << "Вектор, сообщающий оптимальное решение двойственной задаче: " << endl;
	for (auto temp : Y)
		cout << temp << "\t";
	cout << endl << endl;

	cout << "Индексы базиса: " << endl;
	for (auto temp : basis)
		cout << temp << "\t";
	cout << endl << endl;*/

	canon_problem_t canon_problem = ConvertGeneralToCanon(problem);

	column_t X, Y;
	double opt_value;
	comb_t basis;
	tie(opt_value, X, Y, basis) = SolveProblemWithSimplexMethod(canon_problem);

	cout << endl << "Оптимальное значение: " << opt_value << endl << endl;

	cout << "Вектор, сообщающий оптимальное решение прямой задаче: " << endl;
	for (auto temp : X)
		cout << temp << "\t";
	cout << endl << endl;

	cout << "Вектор, сообщающий оптимальное решение двойственной задаче: " << endl;
	for (auto temp : Y)
		cout << temp << "\t";
	cout << endl << endl;

	cout << "Индексы базиса: " << endl;
	for (auto temp : basis)
		cout << temp << "\t";
	cout << endl << endl;

	/*
	


	// ЗАПИСЬ ЗАДАЧИ В ОБЩЕМ ВИДЕ
	general_problem_t problem = FromFileConvertToGeneral("general_task_main.txt");
	matrix_t A; column_t b, c; int M1, N1;
	tie(A, b, c, M1, N1) = problem;
	cout << "Общая задача ЛП:" << endl << "Матрица коэффициентов:" << endl;
	for (int i = 0; i < M1; ++i) {
		for (int j = 0; j < A.size(); ++j) {
			cout << A[j][i] << "\t";
		}
		cout << ">=\t" << b[i];
		cout << endl;
	}
	for (int i = M1; i < A[0].size(); ++i) {
		for (int j = 0; j < A.size(); ++j) {
			cout << A[j][i] << "\t";
		}
		cout << "=\t" << b[i];
		cout << endl;
	}
	cout << endl << "Критерий качества:" << endl;
	for (auto temp : c) {
		cout << temp << "\t";
	}
	cout << endl << "Неотрицательны первые " << N1 << " перем." << endl << endl << endl << endl << endl;



	// ПРЕОБРАЗОВАТЬ ОБЩУЮ ЗАДАЧУ К КАНОНИЧЕСКОЙ 
	canon_problem_t canon_problem = ConvertGeneralToCanon(problem);
	tie(A, b, c) = canon_problem;
	cout << "Каноническая задача ЛП:" << endl << "Матрица коэффициентов:" << endl;
	for (int i = 0; i < A[0].size(); ++i) {
		for (int j = 0; j < A.size(); ++j) {
			cout << A[j][i] << "\t";
		}
		cout << "=\t" << b[i];
		cout << endl;
	}
	cout << endl << "Критерий качества:" << endl;
	for (auto temp : c) {
		cout << temp << "\t";
	}
	cout << endl << endl << endl << endl << endl;



	// ПОЛУЧИЛИ РЕШЕНИЕ ПРЯМОЙ ЗАДАЧИ МЕТОДОМ ПЕРЕБОРА КРАЙНИХ ТОЧЕК
	column_t Xtreme;
	double res;
	comb_t comb;
	general_problem_t problem1 = FromFileConvertToGeneral("general_task_main1.txt");
	canon_problem_t canon_problem1 = ConvertGeneralToCanon(problem1);
	tie(Xtreme, res, comb) = IteratingThroughExtremePoints(canon_problem1);
	cout << "Решение прямой задачи методом перебора крайних точек:" << endl;
	cout << "Оптимальное значение: " << -res << endl;
	cout << "Оптимальная точка: " << endl;
	for (auto temp : Xtreme)
		cout << temp << "\t";
	cout << endl << endl << endl << endl << endl;



	// НАЙДЕМ РЕШЕНИЕ ДВОЙСТВЕННОЙ ЗАДАЧИ ЧЕРЕЗ РЕШЕНИЕ ПРЯМОЙ
	column_t Y = SolvingDualProblem(canon_problem1, Xtreme, comb);
	cout << "Решение двойственной задачи, найденное через решение прямой задачи: " << endl;
	for (auto temp : Y)
		cout << fix(temp, 6) << "\t";
	cout << endl << endl;
	double temp = 0;
	for (int i = 0; i < b.size(); ++i)
		temp += b[i] * Y[i];
	cout << "Значение функции цели двойственной задачи: " << -temp << endl << endl << endl << endl;


	
	// ПОЛУЧИТЬ ДВОЙСТВЕННУЮ ЗАДАЧУ В ОБЩЕМ ВИДЕ
	general_problem_t dual_problem = GetDualLinearProblem(problem);
	tie(A, b, c, M1, N1) = dual_problem;
	cout << "Двойственная задача к общей задаче ЛП:" << endl << "Матрица коэффициентов:" << endl;
	for (int i = 0; i < M1; ++i) {
		for (int j = 0; j < A.size(); ++j) {
			cout << A[j][i] << "\t";
		}
		cout << ">=\t" << b[i];
		cout << endl;
	}
	for (int i = M1; i < A[0].size(); ++i) {
		for (int j = 0; j < A.size(); ++j) {
			cout << A[j][i] << "\t";
		}
		cout << "=\t" << b[i];
		cout << endl;
	}
	cout << endl << "Критерий качества:" << endl;
	for (auto temp : c) {
		cout << temp << "\t";
	}
	cout << endl << "Неотрицательны первые " << N1 << " перем." << endl << endl << endl << endl << endl;



	// ПОЛУЧИЛИ РЕШЕНИЕ ДВОЙСТВЕННОЙ ЗАДАЧИ МЕТОДОМ ПЕРЕБОРА КРАЙНИХ ТОЧЕК
	canon_problem_t canon_dual_problem = ConvertGeneralToCanon(dual_problem);
	tie(Xtreme, res, comb) = IteratingThroughExtremePoints(canon_dual_problem);
	cout << "Решение двойственной задачи методом перебора крайних точек:" << endl;
	cout << "Оптимальное значение: " << res << endl;
	cout << "Оптимальная точка: " << endl;
	for (auto temp : Xtreme)
		cout << temp << "\t";
	cout << endl;



	column_t X, Y;
	double opt_value;
	matrix_t basis;
	tie(opt_value, X, Y, basis) = SolveProblemWithSimplexMethod(canon_problem);


	cout << "Оптимальное значение функции цели: " << opt_value << endl;
	cout << "Оптимальный вектор: " << endl;
	for (auto temp : X)
		cout << temp << "\t";
	cout << endl;
	cout << "Базисные вектора: " << endl;
	for (int j = 0; j < basis[0].size(); ++j) {
		for (int i = 0; i < basis.size(); ++j)
			cout << basis[i][j] << "\t";
		cout << endl;
	}
	cout << endl;
	cout << "Двойственный вектор: " << endl;
	for (auto temp : Y)
		cout << temp << "\t";
	cout << endl;
	*/
	
	return 0;
}
