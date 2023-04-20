#pragma once

#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <tuple>
#include <math.h>

using namespace std;


/*�������� ������, ����������� � ���������� ��� ������� ������� ���������� �������������� � ��������������� ������� � ���*/
class Task
{
public:

	/*�������������� �������*/
	double MinFunc(double x, double y, double z);

	/*������� - �����������*/
	double Limit(double x, double y, double z);

	/*����������� ������� - �����������*/
	vector<double> SubgradientLim(double x, double y, double z);

private:

	/*������� ���� � �������� ���������� �� ��������������*/
	double Phi0(double x, double y);

	/*������ ����������� ����������� � ������ �����*/
	int Index(double x, double y, double z);


};