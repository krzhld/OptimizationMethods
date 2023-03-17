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

	// ÇÀÏÈÑÜ ÇÀÄÀ×È Â ÎÁÙÅÌ ÂÈÄÅ
	/*general_problem_t problem = FromFileConvertToGeneral("general_task_main.txt");
	matrix_t A; column_t b, c; int M1, N1;
	tie(A, b, c, M1, N1) = problem;
	cout << "Îáùàÿ çàäà÷à ËÏ:" << endl << "Ìàòðèöà êîýôôèöèåíòîâ:" << endl;
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
	cout << endl << "Êðèòåðèé êà÷åñòâà:" << endl;
	for (auto temp : c) {
		cout << temp << "\t";
	}
	cout << endl << "Íåîòðèöàòåëüíû ïåðâûå " << N1 << " ïåðåì." << endl << endl << endl << endl << endl;*/



	// ÏÐÅÎÁÐÀÇÎÂÀÒÜ ÎÁÙÓÞ ÇÀÄÀ×Ó Ê ÊÀÍÎÍÈ×ÅÑÊÎÉ 
	/*canon_problem_t canon_problem = ConvertGeneralToCanon(problem);
	tie(A, b, c) = canon_problem;
	cout << "Êàíîíè÷åñêàÿ çàäà÷à ËÏ:" << endl << "Ìàòðèöà êîýôôèöèåíòîâ:" << endl;
	for (int i = 0; i < A[0].size(); ++i) {
		for (int j = 0; j < A.size(); ++j) {
			cout << A[j][i] << "\t";
		}
		cout << "=\t" << b[i];
		cout << endl;
	}
	cout << endl << "Êðèòåðèé êà÷åñòâà:" << endl;
	for (auto temp : c) {
		cout << temp << "\t";
	}
	cout << endl << endl << endl << endl << endl;*/



	
	// ÏÎËÓ×ÈËÈ ÐÅØÅÍÈÅ ÏÐßÌÎÉ ÇÀÄÀ×È ÌÅÒÎÄÎÌ ÏÅÐÅÁÎÐÀ ÊÐÀÉÍÈÕ ÒÎ×ÅÊ
	/*column_t Xtreme;
	double res;
	comb_t comb;
	general_problem_t problem1 = FromFileConvertToGeneral("general_task_main1.txt");
	canon_problem_t canon_problem1 = ConvertGeneralToCanon(problem1);
	tie(Xtreme, res, comb) = IteratingThroughExtremePoints(canon_problem1);
	cout << "Ðåøåíèå ïðÿìîé çàäà÷è ìåòîäîì ïåðåáîðà êðàéíèõ òî÷åê:" << endl;
	cout << "Îïòèìàëüíîå çíà÷åíèå: " << -res << endl;
	cout << "Îïòèìàëüíàÿ òî÷êà: " << endl;
	for (auto temp : Xtreme)
		cout << temp << "\t";
	cout << endl << endl << endl << endl << endl;*/



	// ÍÀÉÄÅÌ ÐÅØÅÍÈÅ ÄÂÎÉÑÒÂÅÍÍÎÉ ÇÀÄÀ×È ×ÅÐÅÇ ÐÅØÅÍÈÅ ÏÐßÌÎÉ
	/*column_t Y = SolvingDualProblem(canon_problem1, Xtreme, comb);
	cout << "Ðåøåíèå äâîéñòâåííîé çàäà÷è, íàéäåííîå ÷åðåç ðåøåíèå ïðÿìîé çàäà÷è: " << endl;
	for (auto temp : Y)
		cout << fix(temp, 6) << "\t";
	cout << endl << endl;
	double temp = 0;
	for (int i = 0; i < b.size(); ++i)
		temp += b[i] * Y[i];
	cout << "Çíà÷åíèå ôóíêöèè öåëè äâîéñòâåííîé çàäà÷è: " << -temp << endl << endl << endl << endl;*/


	// ÏÎËÓ×ÈÒÜ ÄÂÎÉÑÒÂÅÍÍÓÞ ÇÀÄÀ×Ó Â ÎÁÙÅÌ ÂÈÄÅ
	/*general_problem_t dual_problem = GetDualLinearProblem(problem);
	tie(A, b, c, M1, N1) = dual_problem;
	cout << "Äâîéñòâåííàÿ çàäà÷à ê îáùåé çàäà÷å ËÏ:" << endl << "Ìàòðèöà êîýôôèöèåíòîâ:" << endl;
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
	cout << endl << "Êðèòåðèé êà÷åñòâà:" << endl;
	for (auto temp : c) {
		cout << temp << "\t";
	}
	cout << endl << "Íåîòðèöàòåëüíû ïåðâûå " << N1 << " ïåðåì." << endl << endl << endl << endl << endl;*/



	// ÏÎËÓ×ÈËÈ ÐÅØÅÍÈÅ ÄÂÎÉÑÒÂÅÍÍÎÉ ÇÀÄÀ×È ÌÅÒÎÄÎÌ ÏÅÐÅÁÎÐÀ ÊÐÀÉÍÈÕ ÒÎ×ÅÊ
	/*canon_problem_t canon_dual_problem = ConvertGeneralToCanon(dual_problem);
	tie(Xtreme, res, comb) = IteratingThroughExtremePoints(canon_dual_problem);
	cout << "Ðåøåíèå äâîéñòâåííîé çàäà÷è ìåòîäîì ïåðåáîðà êðàéíèõ òî÷åê:" << endl;
	cout << "Îïòèìàëüíîå çíà÷åíèå: " << res << endl;
	cout << "Îïòèìàëüíàÿ òî÷êà: " << endl;
	for (auto temp : Xtreme)
		cout << temp << "\t";
	cout << endl;*/



	/*column_t X, Y;
	double opt_value;
	matrix_t basis;
	tie(opt_value, X, Y, basis) = SolveProblemWithSimplexMethod(canon_problem);


	cout << "Îïòèìàëüíîå çíà÷åíèå ôóíêöèè öåëè: " << opt_value << endl;
	cout << "Îïòèìàëüíûé âåêòîð: " << endl;
	for (auto temp : X)
		cout << temp << "\t";
	cout << endl;
	cout << "Áàçèñíûå âåêòîðà: " << endl;
	for (int j = 0; j < basis[0].size(); ++j) {
		for (int i = 0; i < basis.size(); ++j)
			cout << basis[i][j] << "\t";
		cout << endl;
	}
	cout << endl;
	cout << "Äâîéñòâåííûé âåêòîð: " << endl;
	for (auto temp : Y)
		cout << temp << "\t";
	cout << endl;*/
	
	return 0;
}
