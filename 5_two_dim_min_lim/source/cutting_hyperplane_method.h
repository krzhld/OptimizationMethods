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

/*кортеж из решений прямой и двойственной задачи соответственно*/
typedef tuple<column_t, column_t, comb_t> solving_linear_problem_t;

/*Решение задачи методом отсекающей гиперплоскости
* вход: вектор - начальное приближение, кортеж - первоначальное многогранное множество, данные исходной задачи, вещественное число - точность решения
* выход: кортеж - решение задачи
*/
column_t CuttingHyperplaneMethod(Task& t, double eps);

/*Построение отсекающей гиперплоскости на данной итерации
* Вход: текущее приближение, данные задачи
* Выход: параметры отсекающей гиперплоскости 
*/
hyperplane_t GetCuttingHyperplane(column_t& xk, Task& t);

/*Отсечение гиперплоскости от многогранного множества
* Вход: многогранное множество, отсекаемая гиперплоскость
* Выход: многогранное множество
*/
polyhedron_t AddCuttingHiperplaneInPolyhedron(polyhedron_t& Sk, hyperplane_t& h);

/*Решение задачи линейного программирования на  данной итерации
* Вход: данные задачи, множество ограничений, вектор - решение двойственной задачи с предыдущей итерации
* Выход: кортеж из решений прямой и двойственной задачи
*/
solving_linear_problem_t SolvingLinearProblem(Task& t, polyhedron_t& Sk, column_t& yk, comb_t basis);

/*Решение задачи линейного программирования на первой итерации МОК (когда неизвестно решение двойственной задачи с предыдущей итерации)
* Вход: данные задачи
* Выход: кортеж из решений прямой и двойственной задачи
*/
solving_linear_problem_t SolvingLinearProblemInFirstIter(Task& t);
