#pragma once

#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <tuple>
#include <math.h>

using namespace std;

/*������ ������������ �����*/
typedef vector<double> column_t;

/*������� ������������ �����*/
typedef vector<column_t> matrix_t;

/*������ �� ������� � �������, ������������ ������������ ���������*/
typedef tuple<matrix_t, column_t> polyhedron_t;

/*������ �� ������ ��������������*/
typedef tuple<column_t, double> hyperplane_t;

/*������������ ������������������ ������� �� ������
* ����: ��� �������
* �����: �����
*/
double MultiplipliedVectors(column_t v1, column_t v2);

/*����� ��������
* ����: ��� �������
* �����: ������
*/
column_t SumVector(column_t v1, column_t v2);

/*�������� ��������
* ����: ��� �������
* �����: ������
*/
column_t DiffVector(column_t v1, column_t v2);

/*��������� ������� �� �����
* ����: ������ � �����
* �����: ������
*/
column_t MultVectorAndNumber(column_t vector, double number);


/*����� �������
* ����: ������ ������������ �����
* �����: ������������ �����
*/
double Norm(column_t vector);
