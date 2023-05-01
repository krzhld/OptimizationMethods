#pragma once


#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <tuple>
#include <math.h>
#include "task.h"
#include "utils.h"
#include "framework.h"

using namespace std;

/*������ �� ������� ������ � ������������ ������ ��������������*/
typedef tuple<column_t, column_t, comb_t> solving_linear_problem_t;

/*������� ������ ������� ���������� ��������������
* ����: ������ - ��������� �����������, ������ - �������������� ������������ ���������, ������ �������� ������, ������������ ����� - �������� �������
* �����: ������ - ������� ������
*/
column_t CuttingHyperplaneMethod(Task& t, double eps);

/*���������� ���������� �������������� �� ������ ��������
* ����: ������� �����������, ������ ������
* �����: ��������� ���������� �������������� 
*/
hyperplane_t GetCuttingHyperplane(column_t& xk, Task& t);

/*��������� �������������� �� ������������� ���������
* ����: ������������ ���������, ���������� ��������������
* �����: ������������ ���������
*/
polyhedron_t AddCuttingHiperplaneInPolyhedron(polyhedron_t& Sk, hyperplane_t& h);

/*������� ������ ��������� ���������������� ��  ������ ��������
* ����: ������ ������, ��������� �����������, ������ - ������� ������������ ������ � ���������� ��������
* �����: ������ �� ������� ������ � ������������ ������
*/
solving_linear_problem_t SolvingLinearProblem(Task& t, polyhedron_t& Sk, column_t& yk, comb_t basis);

/*������� ������ ��������� ���������������� �� ������ �������� ��� (����� ���������� ������� ������������ ������ � ���������� ��������)
* ����: ������ ������
* �����: ������ �� ������� ������ � ������������ ������
*/
solving_linear_problem_t SolvingLinearProblemInFirstIter(Task& t);
