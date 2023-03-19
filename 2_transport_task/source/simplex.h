#pragma once
#include "framework.h"

enum class SimplexState { OPTIMAL, UNLIMITED, NEXT };

/**/
void CalculateNextSupportVector(comb_t& Nk_indexes, comb_t& Nk_plus_indexes, column_t& u_k, column_t& cur_X, int N_number, int j_k);

/**/
comb_t GetIndexesPositiveUk(comb_t& Nk_plus_indexes, column_t& u_k);

/**/
void GetLk(comb_t& Lk_indexes, comb_t& Nk_indexes, int N);

/**/
void GetNkPlus(comb_t& Nk_plus_indexes, column_t& cur_X);

/**/
void NormalizeVector(column_t& X);

/**/
SimplexState IterSimplex(matrix_t& main_A, comb_t& Nk_indexes, column_t& main_C, column_t& cur_X, column_t& cur_Y);

/*
* Реализация симплекс метода
* вход: кортеж матриц и столбцов канонической ЗЛП
* выход: вектор, сообщающий МИНИМУМ функции цели
*/
std::tuple<double, column_t, column_t, comb_t> SimplexMethod(canon_problem_t& problem, column_t& cur_X, comb_t& basis);

/*
* Нахождение начального приближения для решения исходной канонической ЗЛП симплекс-методом
* вход: кортеж матриц и столбцов канонической ЗЛП
* выход: допустимый вектор для старта симплекс-метода для исходной задачи
*/
std::tuple<column_t, comb_t> GetInitialApprox(canon_problem_t& problem);

/**/
std::tuple<double, column_t, column_t, comb_t> SolveProblemWithSimplexMethod(canon_problem_t& problem);
