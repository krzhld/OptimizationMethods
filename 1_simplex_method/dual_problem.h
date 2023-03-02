#pragma once

#include <vector>
#include <cstdio>
#include <iostream>
#include <tuple>
#include <cmath>
#include "types.h"
#include "transformations.h"

typedef std::vector<int> comb_t;
typedef std::vector<comb_t> combinations_t;

/*
��������� ������������ ������ 
����: ������ �� ������ �� ������� � �������� ����� ������ ��
�����: ������ �� ������� � �������� ����� ������ ��
*/
general_problem_t getDualLinearProblem(general_problem_t &problem);