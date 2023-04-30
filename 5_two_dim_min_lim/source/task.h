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


/*�������� ������, ����������� � ���������� ��� ������� ������� ���������� �������������� � ��������������� ������� � ���*/
class Task {
public:

	/*�������������� �������*/
	double MinFunc(column_t x);

	/*������� - �����������*/
	double Limit(column_t x);

	/*����������� ������� - �����������*/
	column_t SubgradientLim(column_t x);

	/*������������ ���������, ���������� �������� ��������� �����*/
	polyhedron_t GetS0();

private:

	/*������� ���� � �������� ���������� �� ��������������*/
	double Phi0(column_t x);

	/*������ ����������� ����������� � ������ �����*/
	int Index(column_t x);


};
