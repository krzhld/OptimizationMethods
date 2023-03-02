#pragma once
#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <tuple>

// вектор вещественных чисел
typedef std::vector<double> column_t;

// матрица вещественных чисел
typedef std::vector<column_t> matrix_t;

// кортеж из матрицы A[M,N], векторов b[M], c[N], числа M1 и N1 для общей ЗЛП 
typedef std::tuple<matrix_t, column_t, column_t, int, int> general_problem_t;

// кортеж из матрицы A[M,N], векторов b[M], c[N] канонической ЗЛП
typedef std::tuple<matrix_t, column_t, column_t> canon_problem_t;
