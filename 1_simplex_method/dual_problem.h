#pragma once

#include <vector>
#include <cstdio>
#include <iostream>
#include <tuple>
#include <cmath>
#include "types.h"
#include "transformations.h"

typedef std::vector<int> comb_t;
typedef std::vector<comb_t> combinations_t;

/*
Получение двойственной задачи 
вход: ссылка на кортеж из матрицы и векторов общей задачи ЛП
выход: кортеж из матрицы и векторов общей задачи ЛП
*/
general_problem_t getDualLinearProblem(general_problem_t &problem);