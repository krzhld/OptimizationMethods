#define _CRT_SECURE_NO_WARNINGS
#include "framework.h"

using namespace std;

int main(void) {
	general_problem_t problem = FromFileConvertToGeneral("general_task.txt");
	canon_problem_t canon_problem = ConvertGeneralToCanon(problem);
	column_t cur_X = GetInitialApprox(canon_problem);
	SimplexMethod(canon_problem, cur_X);
	
	return 0;
}
