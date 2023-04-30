#pragma once


#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <tuple>
#include <math.h>
#include "task.h"
#include "utils.h"
#include "framework.h"

using namespace std;

/*решение задачи*/
typedef tuple<column_t, double> solving_t;

/*кортеж из решений пр€мой и двойственной задачи соответственно*/
typedef tuple<column_t, column_t, comb_t> solving_linear_problem_t;

/*–ешение задачи методом отсекающей гиперплоскости
* вход: вектор - начальное приближение, кортеж - первоначальное многогранное множество, данные исходной задачи, вещественное число - точность решени€
* выход: кортеж - решение задачи
*/
solving_t CuttingHyperplaneMethod(Task& t, double eps);

/*ѕостроение отсекающей гиперплоскости на данной итерации
* ¬ход: текущее приближение, данные задачи
* ¬ыход: параметры отсекающей гиперплоскости 
*/
hyperplane_t GetCuttingHyperplane(column_t& xk, Task& t);

/*ќтсечение гиперплоскости от многогранного множества
* ¬ход: многогранное множество, отсекаема€ гиперплоскость
* ¬ыход: многогранное множество
*/
polyhedron_t AddCuttingHiperplaneInPolyhedron(polyhedron_t& Sk, hyperplane_t& h);

/*–ешение задачи линейного программировани€ на  данной итерации
* ¬ход: данные задачи, множество ограничений, вектор - решение двойственной задачи с предыдущей итерации
* ¬ыход: кортеж из решений пр€мой и двойственной задачи
*/
solving_linear_problem_t SolvingLinearProblem(Task& t, polyhedron_t& Sk, column_t& yk, comb_t basis);

/*–ешение задачи линейного программировани€ на первой итерации ћќ  (когда неизвестно решение двойственной задачи с предыдущей итерации)
* ¬ход: данные задачи
* ¬ыход: кортеж из решений пр€мой и двойственной задачи
*/
solving_linear_problem_t SolvingLinearProblemInFirstIter(Task& t);
