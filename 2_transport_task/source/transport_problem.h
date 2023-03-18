#pragma once

#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <tuple>

using namespace std;

// вектор вещественных чисел
typedef vector<double> column_t;

// матрица вещественных чисел
typedef vector<column_t> matrix_t;

// кортеж из столбцов данных о поставщиках, о покупателях и матрицы стоимостей перевозок
typedef tuple<column_t, column_t, matrix_t> transport_problem_t;

// кортеж из матрицы оптимальных объемов перевозок и числа - оптимальной стоимости перевозок
typedef tuple<matrix_t, double> solving_t;

// кортеж из чисел - индексов минимального числа дельта в методе потенциалов и самого минимального дельта 
typedef tuple<int, int, double> min_delta_t;

// кортеж из матрицы A[M,N], векторов b[M], c[N] канонической ЗЛП
typedef std::tuple<matrix_t, column_t, column_t> canon_problem_t;

/*
* Чтение из файла транспортной задачи 
* Пример файла: 
* --------------------------------------
* a:
11 15 10 10 
b: 
11 8 10 9 8
c: 
8 5 7 5 4 
9 6 8 2 3
7 4 7 1 4 
5 4 7 1 2
----------------------------------------
Вход: имя файла
Выход: кортеж из данных транспортной задачи
*/
transport_problem_t ReadFromFileTransportProblem(string filename);

/*
* Решение транспортной задачи методом потенциалов 
* Вход: кортеж из данных транспортной задачи
* Выход: кортеж из данных решения транспортной задачи
*/
solving_t MethodOfPotentials(transport_problem_t problem);

/*Поиск ячеек в цикле перестроения по строке*/
double FindCeelsLine(int inext, int jnext, int imin, int jmin, int N, int M, matrix_t &X, matrix_t& temp, int numb, double teta);

/*Поиск ячеек в цикле перестроения по столбцу*/
double FindCeelsColumn(int inext, int jnext, int imin, int jmin, int N, int M, matrix_t &X, matrix_t& temp, int numb, double teta);

/*Получение из транспортной задачи канонической задачи в классической постановке*/
canon_problem_t GetCanonProblemFromTransportProblem(transport_problem_t& problem);

/*Решение транспортной задачи с учетом возможности дизбаланса данных*/
solving_t SolveTransportProblem(transport_problem_t problem);