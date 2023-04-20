#pragma once

#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <tuple>
#include <math.h>

using namespace std;


/*Исходная задача, приведенная к постановке для решения методом отсекающей гиперплоскости и вспомогательные функции к ней*/
class Task
{
public:

	/*минимизуруемая функция*/
	double MinFunc(double x, double y, double z);

	/*функция - ограничение*/
	double Limit(double x, double y, double z);

	/*субградиент функции - ограничения*/
	vector<double> SubgradientLim(double x, double y, double z);

private:

	/*функция цели в исходной постановке до преобразования*/
	double Phi0(double x, double y);

	/*индекс наибольшего ограничения в данной точке*/
	int Index(double x, double y, double z);


};