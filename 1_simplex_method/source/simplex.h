#pragma once
#include "transformations.h"

typedef std::tuple<matrix_t, std::vector<int>> simplex_table_t;

void iterSimplexMethod(canon_problem_t& problem, column_t& first_approx);

/*
���������� �������� ������

����: ������ ������ � �������� ������������ ���

�����: ������, ���������� ������� ������� ����
*/
column_t simplexMethod(canon_problem_t& problem, column_t& first_approx);


/*
���������� ���������� ����������� ��� ������� �������� 
������������ ��� �������� �������

����: ������ ������ � �������� ������������ ���

�����: ���������� ������ ��� ������ ��������-������ ��� �������� ������
*/
column_t getInitialApprox(canon_problem_t& problem);

double getOptimalValue(canon_problem_t& problem, column_t solution);
