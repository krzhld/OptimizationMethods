#pragma once
#include <vector>
#include <cstdio>
#include <iostream>
#include <tuple>
#include <cmath>
#include "types.h"


typedef std::pair<column_t, double> solving_t;
typedef std::vector<int> comb_t;
typedef std::vector<comb_t> combinations_t;

/*
������� ������ ������� �������� ������� �����
����: ������ �� ������ �� ������� � �������� ������������ ������ ��
�����: ������ �� ������� (�������, � ������� ����������� ����������� �������) 
� ����� (�������� ������� ���� ��� �������, � ������� ����������� ����������� �������)
*/
solving_t IteratingThroughExtremePoints(canon_problem_t &problem);
