#include "bruteforce.h"
using namespace std;

solution_t IteratingThroughExtremePoints(canon_problem_t& problem) {
	matrix_t A;
	column_t b, c;
	tie(A, b, c) = problem;
	int n = size(c), m = size(b);

	/* проверка корректности задачи */
	if (n <= m) {
		cout << "Incorred problem";
		exit(-1);
	}

	/* получение всех комбинаций из n (количество переменных) по m (количество ограничений) */
	combinations_t combinations = CombinationsWithoutRepetitions(m, n);

	column_t maxPoint; // вектор, в котором достигается оптимальное решение
	maxPoint.resize(n);
	double max = -1e20; // максимальное значение функции цели
	comb_t optComb; // базисная комбинация, при которой достигается оптимальное решение
	optComb.resize(m);


	/* перебор всех возможных комбинаций сочетаний */
	for (auto& comb : combinations) {

		double valueInPoint; // значение функции цели при данной комбинации
		column_t X; // опорный вектор при данной комбинации 
		X.resize(n);
		column_t columnA1; //матрица для решения системы уравнений, из которой отброшены столбцы , номера которых не входят в данную комбинацию 

		/*размерность матрицы A1 m на m*/
		columnA1.resize(m);
		matrix_t A1;
		A1.resize(m);

		for (int k = 0; k < m; k++)
			A1[k] = columnA1;

		column_t solveSystem; //решение СЛАУ с матрицей A1 и вектором b, имеет размерность 
		solveSystem.resize(m);


		int numb = 0; // счетсчик для задания компонент матрицы A1

		/* задание компонент матрицы A1, берутся те столбцы, номера которых есть в данной комбинации */
		for (int i = 0; i < n; i++) {
			if (IsNumberInCombination(i + 1, comb)) {
				for (int j = 0; j < m; j++)
					A1[j][numb] = A[i][j];

				numb++;
			}
		}

		/* те компоненты опорного вектора, номера которых не входят в данную комбинацию принимают нулевые значения */
		for (int i = 0; i < n; i++) {
			if (!IsNumberInCombination(i + 1, comb))
				X[i] = 0;
		}

		// проверка невырожденности матрицы 
		if (Determinant(A1) != 0) {
			/*решение СЛАУ методом вращений*/
			solveSystem = RotationMethod(A1, b);

			// проверка на неотрицательность компонент решения СЛАУ (иначе вектор не будет допустимым)
			if (NonNegativityOfVector(solveSystem)) {
				int numbInSolv = 0; //счетчик 

			/*те компоненты опорного вектора, номера которых входят в данную комбинацию принимают соответствующие значения решения СЛАУ*/
				for (int k = 0; k < n; k++) {
					if (IsNumberInCombination(k + 1, comb)) {
						X[k] = solveSystem[numbInSolv];
						numbInSolv++;
					}
				}

				/*вычисление значения целевой функции при данном опорном векторе*/
				valueInPoint = MultipliedVectors(c, X);

				/*если полученное значение функции цели меньше минимума, то оно становится минимумом,
				вектору, в котором достигается минимум присваивается значение X*/
				if (valueInPoint > max) {
					maxPoint = X;
					max = valueInPoint;
					optComb = comb;
				}
			}
		}
	}
	/*создание кортежа из вектора, в котором достигается оптимальное решение и оптимального решения*/
	solution_t solving = make_tuple(maxPoint, max, optComb);

	return solving;
}
