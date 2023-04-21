#pragma once

#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <tuple>
#include <math.h>

using namespace std;

/*вектор вещественных чисел*/
typedef vector<double> column_t;

/*матрица вещественных чисел*/
typedef vector<column_t> matrix_t;

/*Перемножение транспонированного вектора на вектор
* Вход: два вектора
* Выход: число
*/
double MultiplipliedVectors(column_t v1, column_t v2);

/*Сумма векторов
* Вход: два вектора
* Выход: вектор
*/
column_t SumVector(column_t v1, column_t v2);

/*Разность векторов
* Вход: два вектора
* Выход: вектор
*/
column_t DiffVector(column_t v1, column_t v2);

/*Умножение вектора на число
* Вход: вектор и число
* Выход: вектор
*/
column_t MultVectorAndNumber(column_t vector, double number);


/*Норма вектора
* Вход: вектор вещественных чисел
* Выход: вещественное число
*/
double Norm(column_t vector);