#pragma once
#define _CRT_SECURE_NO_WARNINGS


#include <vector>
#include <cstdio>
#include <iostream>
#include <sstream>
#include <fstream>
#include <tuple>
#include <cmath>

enum class SimplexState {OPTIMAL, UNLIMITED, NEXT};

/* ������ ������������ ����� */
typedef std::vector<double> column_t;

/* ������� ������������ ����� */
typedef std::vector<column_t> matrix_t;

/* ������ ���������� ��������� */
typedef std::vector<int> comb_t;

/* ������������ �������� ��������� */
typedef std::vector<comb_t> combinations_t;


/* ������ �� ������� A[M,N], �������� b[M], c[N], ����� M1 � N1 ��� ����� ��� */
typedef std::tuple<matrix_t, column_t, column_t, int, int> general_problem_t;

/* ������ �� ������� A[M,N], �������� b[M], c[N] ��� ������������ ��� */
typedef std::tuple<matrix_t, column_t, column_t> canon_problem_t;

/* ������ �� ������� ������� � ������������ �������� ������� ���� */
typedef std::tuple<column_t, double, comb_t> solution_t;


/* 
* ��������� ���������� � ������ ��������� ��� ����������
* ����:
* �����:
*/
bool NextSet(comb_t& a, int n, int m);

/* 
* ����� ���� ��������� �� n �� m ��� ���������� 
* ����: 
* �����:
*/
combinations_t CombinationsWithoutRepetitions(int m, int n);

/* 
* �������� ���������� ����� � ������ �� �������� � ���������
* ����:
* �����:
*/
bool IsNumberInCombination(int number, comb_t comb);

/* 
* �������� ����������������� ��������� ������� 
* ����: ������
* �����: ������ ��������
*/
bool NonNegativityOfVector(column_t v);

/*
* �������� ����������������� ��������� �������
* ����: ������
* �����: ������ ��������
*/
bool NonPositivityOfVector(column_t v);

/* 
* ������������ ������� �� ������
* ����: ���������� ������� � ���������� �� ������� ������
* �����: ������
*/
column_t MultipliedMatrixAndColumn(matrix_t M, column_t c);

/*
* ���������� ����� �������
* ����: �������
* �����: ���� ������� - ����� �����
*/
//int RankMatrix(matrix_t A);

/*
* ��������� �������� ������� ������� ������
* ����: ���������� �������
* �����: ���������� �������
*/
matrix_t InverseMatrix(matrix_t A);

/*
* ����� �������� ������� ����
* ����: ������� � ������ ������� ����
* �����: ������� ������� ����
*/
column_t RotationMethod(matrix_t A, column_t b);

/*
* ����� ������ ������� ����
* ����: ������� � ������ ������� ����
* �����: ������� ������� ����
*/
column_t GaussMethod(matrix_t A, column_t b);

/* 
* ��������� ������������������ ������� �� ������ ����� �� �����������
* ����: ��� �������
* �����: �� ��������� ������������
*/
double MultipliedVectors(column_t v1, column_t v2);

/* 
* ������������ �������
* ����: ���������� �������
* �����: � ������������
*/
double Determinant(matrix_t matrix);

/*
������ �� ����� ������ �� � ������������� � � ����� ������ �� 
������ �����:
----------------------------------------------------------------
min_func:
-5	-4
>=
<=
6	4	24
1	2	6
-1	1	1
0	1	2
==
x_i >= 0
1	1
----------------------------------------------------------------
����: ���� �����
�����: ������ �� ������� A[M,N], �������� b[M], c[N], ����� M1 � N1:
������ M1 ����� ������� A - �����������,
������ N1 ��������� C[N] ������������� x_i >= 0
*/
general_problem_t FromFileConvertToGeneral(std::string filename);

/*
* ������������� ����� ������ �� � ������������ ������ ��
* ����: ������ �� ������ ������, �������� � ����� �� ����� ������ ��
* �����: ������ �� ������� A[M,N], �������� b[M] >= 0, c[N] ������������ ���
*/
canon_problem_t ConvertGeneralToCanon(general_problem_t& problem);

void TransposeQuadMatrix(matrix_t& matrix);

/* 
* ���������������� � ��������� �� -1 ������� ��� �������� ������������ ������ � ������������
* ����: �������
* �����: ����������������� �������, ���������� �� -1
*/
matrix_t TransformationMatrix(matrix_t matrix);

/*
* ��������� ������������ ������ � ����� ���
* ����: ������ �� ������ �� ������� � �������� ����� ������ ��
* �����: ������ �� ������� � �������� ����� ������ ��
*/
general_problem_t GetDualLinearProblem(general_problem_t& problem);

/*
* ������� ������ ������� �������� ������� �����
* ����: ������ �� ������ �� ������� � �������� ������������ ������ ��
* �����: ������ �� ������� (�������, � ������� ����������� ����������� �������)
* � ����� (�������� ������� ���� ��� �������, � ������� ����������� ����������� �������)
*/
solution_t IteratingThroughExtremePoints(canon_problem_t& problem);

/*
* �������������� ������� ������������ ������ �� ������ ������� ������ ������
* ����: ������ �� ������� � �������� ������������ ������ ��, ������ - ������� ������ ������, ������ ����� ����� - ������ �������� ��������� ������� ������ ������
* �����: ������ - ������� ������������ ������
*/
column_t SolvingDualProblem(canon_problem_t problem, column_t X, comb_t optBasis);

/*
* �������� ���������� ������� �� ������������� �� �������� ��������
* ����:
* �����:
*/
void GetCurMatrix(matrix_t& A, matrix_t& cur_matrix, comb_t& cur_columns);

/**/
void CalculateNextSupportVector(comb_t& Nk_indexes, comb_t& Nk_plus_indexes, column_t& u_k, column_t& cur_X, int N_number, int j_k);

/**/
comb_t GetIndexesPositiveUk(comb_t& Nk_plus_indexes, column_t& u_k);

/**/
void GetLk(comb_t& Lk_indexes, comb_t& Nk_indexes, int N);

/**/
comb_t GetDifferenceSets(comb_t& A, comb_t& B);

/**/
SimplexState IterSimplex(matrix_t& main_A, int M_number, comb_t& Nk_indexes, column_t& main_C, column_t& cur_X, column_t& cur_Y);

/*
* ���������� �������� ������
* ����: ������ ������ � �������� ������������ ���
* �����: ������, ���������� ������� ������� ����
*/
std::tuple<double, column_t, column_t, comb_t> SimplexMethod(canon_problem_t& problem, column_t& cur_X, comb_t& basis);

/*
* ���������� ���������� ����������� ��� ������� �������� ������������ ��� ��������-�������
* ����: ������ ������ � �������� ������������ ���
* �����: ���������� ������ ��� ������ ��������-������ ��� �������� ������
*/
std::tuple<column_t, comb_t> GetInitialApprox(canon_problem_t& problem);

std::tuple<double, column_t, column_t, matrix_t> SolveProblemWithSimplexMethod(canon_problem_t& problem);
