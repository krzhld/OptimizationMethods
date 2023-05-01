#pragma once
#define _CRT_SECURE_NO_WARNINGS

#include <algorithm>
#include <vector>
#include <cstdio>
#include <iostream>
#include <sstream>
#include <fstream>
#include <tuple>
#include <cmath>

/* Вектор вещественных чисел */
typedef std::vector<double> column_t;

/* Матрица вещественных чисел */
typedef std::vector<column_t> matrix_t;

/* Вектор некоторого сочетания */
typedef std::vector<int> comb_t;

/* Совокупность векторов сочетаний */
typedef std::vector<comb_t> combinations_t;


/* Кортеж из матрицы A[M,N], векторов b[M], c[N], числа M1 и N1 для общей ЗЛП */
typedef std::tuple<matrix_t, column_t, column_t, int, int> general_problem_t;

/* Кортеж из матрицы A[M,N], векторов b[M], c[N] для канонической ЗЛП */
typedef std::tuple<matrix_t, column_t, column_t> canon_problem_t;

/* Кортеж из вектора решения и оптимального значения функции цели */
typedef std::tuple<column_t, double, comb_t> solution_t;


/* 
* Следующая комбинация в поиске сочетаний без повторений
* вход:
* выход:
*/
bool NextSet(comb_t& a, int n, int m);

/**/
void NextSetCycle(comb_t& a, int n, int m);

/* 
* Поиск всех сочетаний из n по m без повторений 
* вход: 
* выход:
*/
combinations_t CombinationsWithoutRepetitions(int m, int n);

/* 
* Проверка совпадения числа и одного из значений в сочетании
* вход:
* выход:
*/
bool IsNumberInCombination(int number, comb_t comb);

/**/
comb_t GetDifferenceSets(comb_t& A, comb_t& B);

/* 
* Проверка неотрицательности компонент вектора 
* вход: вектор
* выход: булево значение
*/
bool NonNegativityOfVector(column_t v);

/*
* Проверка неположительности компонент вектора
* вход: вектор
* выход: булево значение
*/
bool NonPositivityOfVector(column_t v);

/* 
* Перемножение матрицы на вектор
* вход: квадратная матрица и подходящий по размеру вектор
* выход: вектор
*/
column_t MultipliedMatrixAndColumn(matrix_t M, column_t c);

/*
* Является ли матрица полноранговой (считаем, что столбцов больше чем строк)
* вход: матрица
* выход: булевое значение
*/
bool IsFullRankMatrix(matrix_t& A);

/*
* Получение обратной матрицы методом Гаусса
* вход: квадратная матрица
* выход: квадратная матрица
*/
matrix_t InverseMatrix(matrix_t A);

/*
* Метод вращений решения СЛАУ
* вход: матрица и правый столбец СЛАУ
* выход: столбец решения СЛАУ
*/
column_t RotationMethod(matrix_t A, column_t b);

/*
* Метод Гаусса решения СЛАУ
* вход: матрица и правый столбец СЛАУ
* выход: столбец решения СЛАУ
*/
column_t GaussMethod(matrix_t A, column_t b);

/* 
* Умножение транспонированного вектора на вектор такой же размерности
* вход: два вектора
* выход: их скалярное произведение
*/
double MultipliedVectors(column_t v1, column_t v2);

/* 
* Определитель матрицы
* вход: квадратная матрица
* выход: её определитель
*/
double Determinant(matrix_t matrix);

/*
читать из файла задачу ЛП и преобразовать её к общей задаче ЛП 
пример файла:
----------------------------------------------------------------
min_func:
-5	-4
>=
<=
6	4	24
1	2	6
-1	1	1
0	1	2
==
x_i >= 0
1	1
----------------------------------------------------------------
вход: путь файла
выход: кортеж из матрицы A[M,N], векторов b[M], c[N], числа M1 и N1:
первые M1 строк матрицы A - неравенства,
первые N1 координат C[N] соответствуют x_i >= 0
*/
general_problem_t FromFileConvertToGeneral(std::string filename);

/*
* Преобразовать общую задачу ЛП в каноническую задачу ЛП
* вход: ссылка на кортеж матриц, векторов и чисел из общей задачи ЛП
* выход: кортеж из матрицы A[M,N], векторов b[M] >= 0, c[N] канонической ЗЛП
*/
canon_problem_t ConvertGeneralToCanon(general_problem_t& problem);

void TransposeQuadMatrix(matrix_t& matrix);

/*
* Получить квадратную матрицу из прямоугольной по индексам столбцов
* вход:
* выход:
*/
void GetCurMatrix(matrix_t& A, matrix_t& cur_matrix, comb_t& cur_columns);

void FreeMatrix(matrix_t& A);
