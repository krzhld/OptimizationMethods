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
ѕолучение двойственной задачи 
вход: ссылка на кортеж из матрицы и векторов общей задачи Ћѕ
выход: кортеж из матрицы и векторов общей задачи Ћѕ
*/
general_problem_t getDualLinearProblem(general_problem_t &problem);

/*
* ¬осстановление решени€ двойственной задачи на основе решени€ пр€мой задачи 
* вход: кортеж из матрицы и веткоров канонической задачи Ћѕ, вектор - решение пр€мой задачи, вектор целых чисел - номера базисных компонент решени€ пр€мой задачи 
* выход: вектор - решение двойственной задачи
*/
column_t SolvingDualProblem(canon_problem_t problem, column_t X, comb_t optBasis);


/*
* ѕолучение обратной матрицы 
* вход: квадратна€ матрица 
* выход: квадратна€ матрица
*/
matrix_t InverseMatrix(matrix_t A);
