#pragma once
#include "transformations.h"

typedef std::tuple<matrix_t, std::vector<int>> simplex_table_t;

void iterSimplexMethod(canon_problem_t& problem, column_t& first_approx);

/*
реализация симплекс метода

вход: кортеж матриц и столбцов канонической ЗЛП

выход: вектор, сообщающий МИНИМУМ функции цели
*/
column_t simplexMethod(canon_problem_t& problem, column_t& first_approx);


/*
нахождение начального приближения для решения исходной 
канонической ЗЛП симплекс методом

вход: кортеж матриц и столбцов канонической ЗЛП

выход: допустимый вектор для старта симплекс-метода для исходной задачи
*/
column_t getInitialApprox(canon_problem_t& problem);

double getOptimalValue(canon_problem_t& problem, column_t solution);
