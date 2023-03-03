#pragma once
#include "types.h"


/*
читать из файла задачу ЛП и преобразовать её к общей задаче ЛП 

пример файла:
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


вход: путь файла

выход: кортеж из матрицы A[M,N], векторов b[M], c[N], числа M1 и N1:
первые M1 строк матрицы A - неравенства, 
первые N1 координат C[N] соответствуют x_i >= 0
*/
general_problem_t fromFileConvertToGeneral(std::string filename);


/*
преобразовать общую задачу ЛП в каноническую задачу ЛП

вход: ссылка на кортеж матриц, векторов и чисел из общей задачи ЛП

выход: кортеж из матрицы A[M,N], векторов b[M] >= 0, c[N] канонической ЗЛП
*/
canon_problem_t convertGeneralToCanon(general_problem_t& problem);