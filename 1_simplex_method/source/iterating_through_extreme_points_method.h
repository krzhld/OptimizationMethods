#pragma once
#include <vector>
#include <cstdio>
#include <iostream>
#include <tuple>
#include <cmath>
#include "types.h"


typedef std::vector<int> comb_t;
typedef std::vector<comb_t> combinations_t;
typedef std::tuple<column_t, double, comb_t> solving_t;

/*
������� ������ ������� �������� ������� �����
����: ������ �� ������ �� ������� � �������� ������������ ������ ��
�����: ������ �� ������� (�������, � ������� ����������� ����������� �������) 
� ����� (�������� ������� ���� ��� �������, � ������� ����������� ����������� �������)
*/
solving_t IteratingThroughExtremePoints(canon_problem_t &problem);

/*
* ���������� ����� �������
* ����: ������� � �� �������
* �����: ����� �����
*/
int RankMatrix(matrix_t matrix, int N, int M);