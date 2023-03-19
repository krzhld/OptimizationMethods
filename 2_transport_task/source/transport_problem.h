#pragma once

#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <tuple>

using namespace std;

// ������ ������������ �����
typedef vector<double> column_t;

// ������� ������������ �����
typedef vector<column_t> matrix_t;

// ������ �� �������� ������ � �����������, � ����������� � ������� ���������� ���������, � ����� ������� ������ � ������� �� ������������ 
typedef tuple<column_t, column_t, matrix_t, column_t> transport_problem_t;

// ������ �� ������� ����������� ������� ��������� � ����� - ����������� ��������� ���������
typedef tuple<matrix_t, double> solving_t;

// ������ �� ����� - �������� ������������ ����� ������ � ������ ����������� � ������ ������������ ������ 
typedef tuple<int, int, double> min_delta_t;

// ������ �� ������� A[M,N], �������� b[M], c[N] ������������ ���
typedef std::tuple<matrix_t, column_t, column_t> canon_problem_t;

/*
* ������ �� ����� ������������ ������ 
* ������ �����: 
* --------------------------------------
* a:
11 15 10 10 
b: 
11 8 10 9 8
c: 
8 5 7 5 4 
9 6 8 2 3
7 4 7 1 4 
5 4 7 1 2
----------------------------------------
����: ��� �����
�����: ������ �� ������ ������������ ������
*/
transport_problem_t ReadFromFileTransportProblem(string filename);

/*
* ������� ������������ ������ ������� ����������� 
* ����: ������ �� ������ ������������ ������
* �����: ������ �� ������ ������� ������������ ������
*/
solving_t MethodOfPotentials(transport_problem_t problem);

/*����� ����� � ����� ������������ �� ������*/
double FindCeelsLine(int inext, int jnext, int imin, int jmin, int N, int M, matrix_t &X, matrix_t& temp, int numb, double teta);

/*����� ����� � ����� ������������ �� �������*/
double FindCeelsColumn(int inext, int jnext, int imin, int jmin, int N, int M, matrix_t &X, matrix_t& temp, int numb, double teta);

/*��������� �� ������������ ������ ������������ ������ � ������������ ����������*/
canon_problem_t GetCanonProblemFromTransportProblem(transport_problem_t& problem);

/*������� ������������ ������ � ������ ����������� ���������� ������*/
solving_t SolveTransportProblem(transport_problem_t problem);

/*����� ��������� ������� */
double SumColumn(column_t col);

/*����� ��������� ���������*/
double FindCost(matrix_t& c, matrix_t& X, int N, int M);