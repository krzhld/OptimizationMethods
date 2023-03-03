#pragma once
#include <vector>
#include <cstdio>
#include <iostream>
#include <tuple>
#include <cmath>
#include "types.h"


typedef std::vector<int> comb_t;
typedef std::vector<comb_t> combinations_t;
typedef std::tuple<column_t, double, comb_t> solving_t;

/*
Решение задачи методом перебора крайних точек
вход: ссылка на кортеж из матрицы и векторов канонической задачи ЛП
выход: кортеж из вектора (вектора, в котором достигается оптимальное решение) 
и число (значение функции цели для вектора, в котором достигается оптимальное решение)
*/
solving_t IteratingThroughExtremePoints(canon_problem_t &problem);

/*
* Нахождение ранга матрицы
* вход: матрица и ее размеры
* выход: целое число
*/
int RankMatrix(matrix_t matrix, int N, int M);