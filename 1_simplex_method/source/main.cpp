#define _CRT_SECURE_NO_WARNINGS
#include "framework.h"
#include<windows.h>

using namespace std;

int main(void) {
	setlocale(LC_ALL, "Russian");
	SetConsoleOutputCP(866);

	general_problem_t problem = FromFileConvertToGeneral("general_task_inet.txt");
	canon_problem_t canon_problem = ConvertGeneralToCanon(problem);


	/*column_t Xtreme;
	double res;
	comb_t comb;
	tie(Xtreme, res, comb) = IteratingThroughExtremePoints(canon_problem);*/

	column_t X, Y;
	double opt_value;
	matrix_t basis;
	tie(opt_value, X, Y, basis) = SolveProblemWithSimplexMethod(canon_problem);

	/*cout << "Оптимальное значение функции цели: " << opt_value << endl;
	cout << "Оптимальный вектор: " << endl;
	for (auto temp : X)
		cout << temp << "\t";
	cout << endl;
	cout << "Базисные вектора: " << endl;
	for (int j = 0; j < basis[0].size(); ++j) {
		for (int i = 0; i < basis.size(); ++j)
			cout << basis[i][j] << "\t";
		cout << endl;
	}
	cout << endl;
	cout << "Двойственный вектор: " << endl;
	for (auto temp : Y)
		cout << temp << "\t";
	cout << endl;*/
	
	return 0;
}
