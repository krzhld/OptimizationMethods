#define _CRT_SECURE_NO_WARNINGS
#include "framework.h"
#include <iomanip>
#include <windows.h>

using namespace std;

string fix(float x, int p) {
	ostringstream strout;
	strout << fixed << setprecision(p) << x;
	string str = strout.str();
	size_t end = str.find_last_not_of('0') + 1;
	return str.erase(end);
}

int main(void) {
	setlocale(LC_ALL, "Russian");
	SetConsoleOutputCP(866);

	general_problem_t problem = FromFileConvertToGeneral("general_task_studfile.txt");
	canon_problem_t canon_problem = ConvertGeneralToCanon(problem);

	column_t X, Y;
	double opt_value;
	matrix_t basis;
	tie(opt_value, X, Y, basis) = SolveProblemWithSimplexMethod(canon_problem);

	// ������ ������ � ����� ����
	/*general_problem_t problem = FromFileConvertToGeneral("general_task_main.txt");
	matrix_t A; column_t b, c; int M1, N1;
	tie(A, b, c, M1, N1) = problem;
	cout << "����� ������ ��:" << endl << "������� �������������:" << endl;
	for (int i = 0; i < M1; ++i) {
		for (int j = 0; j < A.size(); ++j) {
			cout << A[j][i] << "\t";
		}
		cout << ">=\t" << b[i];
		cout << endl;
	}
	for (int i = M1; i < A[0].size(); ++i) {
		for (int j = 0; j < A.size(); ++j) {
			cout << A[j][i] << "\t";
		}
		cout << "=\t" << b[i];
		cout << endl;
	}
	cout << endl << "�������� ��������:" << endl;
	for (auto temp : c) {
		cout << temp << "\t";
	}
	cout << endl << "�������������� ������ " << N1 << " �����." << endl << endl << endl << endl << endl;*/



	// ������������� ����� ������ � ������������ 
	/*canon_problem_t canon_problem = ConvertGeneralToCanon(problem);
	tie(A, b, c) = canon_problem;
	cout << "������������ ������ ��:" << endl << "������� �������������:" << endl;
	for (int i = 0; i < A[0].size(); ++i) {
		for (int j = 0; j < A.size(); ++j) {
			cout << A[j][i] << "\t";
		}
		cout << "=\t" << b[i];
		cout << endl;
	}
	cout << endl << "�������� ��������:" << endl;
	for (auto temp : c) {
		cout << temp << "\t";
	}
	cout << endl << endl << endl << endl << endl;*/



	
	// �������� ������� ������ ������ ������� �������� ������� �����
	/*column_t Xtreme;
	double res;
	comb_t comb;
	general_problem_t problem1 = FromFileConvertToGeneral("general_task_main1.txt");
	canon_problem_t canon_problem1 = ConvertGeneralToCanon(problem1);
	tie(Xtreme, res, comb) = IteratingThroughExtremePoints(canon_problem1);
	cout << "������� ������ ������ ������� �������� ������� �����:" << endl;
	cout << "����������� ��������: " << -res << endl;
	cout << "����������� �����: " << endl;
	for (auto temp : Xtreme)
		cout << temp << "\t";
	cout << endl << endl << endl << endl << endl;*/



	// ������ ������� ������������ ������ ����� ������� ������
	/*column_t Y = SolvingDualProblem(canon_problem1, Xtreme, comb);
	cout << "������� ������������ ������, ��������� ����� ������� ������ ������: " << endl;
	for (auto temp : Y)
		cout << fix(temp, 6) << "\t";
	cout << endl << endl;
	double temp = 0;
	for (int i = 0; i < b.size(); ++i)
		temp += b[i] * Y[i];
	cout << "�������� ������� ���� ������������ ������: " << -temp << endl << endl << endl << endl;*/


	// �������� ������������ ������ � ����� ����
	/*general_problem_t dual_problem = GetDualLinearProblem(problem);
	tie(A, b, c, M1, N1) = dual_problem;
	cout << "������������ ������ � ����� ������ ��:" << endl << "������� �������������:" << endl;
	for (int i = 0; i < M1; ++i) {
		for (int j = 0; j < A.size(); ++j) {
			cout << A[j][i] << "\t";
		}
		cout << ">=\t" << b[i];
		cout << endl;
	}
	for (int i = M1; i < A[0].size(); ++i) {
		for (int j = 0; j < A.size(); ++j) {
			cout << A[j][i] << "\t";
		}
		cout << "=\t" << b[i];
		cout << endl;
	}
	cout << endl << "�������� ��������:" << endl;
	for (auto temp : c) {
		cout << temp << "\t";
	}
	cout << endl << "�������������� ������ " << N1 << " �����." << endl << endl << endl << endl << endl;*/



	// �������� ������� ������������ ������ ������� �������� ������� �����
	/*canon_problem_t canon_dual_problem = ConvertGeneralToCanon(dual_problem);
	tie(Xtreme, res, comb) = IteratingThroughExtremePoints(canon_dual_problem);
	cout << "������� ������������ ������ ������� �������� ������� �����:" << endl;
	cout << "����������� ��������: " << res << endl;
	cout << "����������� �����: " << endl;
	for (auto temp : Xtreme)
		cout << temp << "\t";
	cout << endl;*/



	/*column_t X, Y;
	double opt_value;
	matrix_t basis;
	tie(opt_value, X, Y, basis) = SolveProblemWithSimplexMethod(canon_problem);


	cout << "����������� �������� ������� ����: " << opt_value << endl;
	cout << "����������� ������: " << endl;
	for (auto temp : X)
		cout << temp << "\t";
	cout << endl;
	cout << "�������� �������: " << endl;
	for (int j = 0; j < basis[0].size(); ++j) {
		for (int i = 0; i < basis.size(); ++j)
			cout << basis[i][j] << "\t";
		cout << endl;
	}
	cout << endl;
	cout << "������������ ������: " << endl;
	for (auto temp : Y)
		cout << temp << "\t";
	cout << endl;*/
	
	return 0;
}
