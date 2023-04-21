#pragma once

#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <tuple>
#include <math.h>
#include "utils.h"

using namespace std;


/*Исходная задача, приведенная к постановке для решения методом отсекающей гиперплоскости и вспомогательные функции к ней*/
class Task
{
public:

	/*минимизуруемая функция*/
	double MinFunc(column_t x);

	/*функция - ограничение*/
	double Limit(column_t x);

	/*субградиент функции - ограничения*/
	column_t SubgradientLim(column_t x);

private:

	/*функция цели в исходной постановке до преобразования*/
	double Phi0(column_t x);

	/*индекс наибольшего ограничения в данной точке*/
	int Index(column_t x);


};