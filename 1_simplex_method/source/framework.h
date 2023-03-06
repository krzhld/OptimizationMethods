#pragma once
#define _CRT_SECURE_NO_WARNINGS


#include <vector>
#include <cstdio>
#include <iostream>
#include <sstream>
#include <fstream>
#include <tuple>
#include <cmath>

enum class SimplexState {OPTIMAL, UNLIMITED, NEXT};

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
* Нахождение ранга матрицы
* вход: матрица
* выход: ранг матрицы - целое число
*/
//int RankMatrix(matrix_t A);

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
* Транспонирование и умножение на -1 матрицы для сведения канонической задачи к двойственной
* вход: матрица
* выход: транспонированная матрица, умноженная на -1
*/
matrix_t TransformationMatrix(matrix_t matrix);

/*
* Получение двойственной задачи к общей ЗЛП
* вход: ссылка на кортеж из матрицы и векторов общей задачи ЛП
* выход: кортеж из матрицы и векторов общей задачи ЛП
*/
general_problem_t GetDualLinearProblem(general_problem_t& problem);

/*
* Решение задачи методом перебора крайних точек
* вход: ссылка на кортеж из матрицы и векторов канонической задачи ЛП
* выход: кортеж из вектора (вектора, в котором достигается оптимальное решение)
* и число (значение функции цели для вектора, в котором достигается оптимальное решение)
*/
solution_t IteratingThroughExtremePoints(canon_problem_t& problem);

/*
* Восстановление решения двойственной задачи на основе решения прямой задачи
* вход: кортеж из матрицы и веткоров канонической задачи ЛП, вектор - решение прямой задачи, вектор целых чисел - номера базисных компонент решения прямой задачи
* выход: вектор - решение двойственной задачи
*/
column_t SolvingDualProblem(canon_problem_t problem, column_t X, comb_t optBasis);

/*
* Получить квадратную матрицу из прямоугольной по индексам столбцов
* вход:
* выход:
*/
void GetCurMatrix(matrix_t& A, matrix_t& cur_matrix, comb_t& cur_columns);

/**/
void CalculateNextSupportVector(comb_t& Nk_indexes, comb_t& Nk_plus_indexes, column_t& u_k, column_t& cur_X, int N_number, int j_k);

/**/
comb_t GetIndexesPositiveUk(comb_t& Nk_plus_indexes, column_t& u_k);

/**/
void GetLk(comb_t& Lk_indexes, comb_t& Nk_indexes, int N);

/**/
comb_t GetDifferenceSets(comb_t& A, comb_t& B);

/**/
SimplexState IterSimplex(matrix_t& main_A, int M_number, comb_t& Nk_indexes, column_t& main_C, column_t& cur_X, column_t& cur_Y);

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

std::tuple<double, column_t, column_t, matrix_t> SolveProblemWithSimplexMethod(canon_problem_t& problem);
