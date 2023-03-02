#define _CRT_SECURE_NO_WARNINGS
#include "transformations.h"
#include "simplex.h"
/*#include "iterating_through_extreme_points_method.h"
#include "dual_problem.h"*/

int main(void) {
	general_problem_t problem = fromFileConvertToGeneral("general_task.txt");
	canon_problem_t canon_problem = convertGeneralToCanon(problem);
	column_t first_approx = { 0, 0, 24, 6, 1, 2 };
	column_t solution = simplexMethod(canon_problem, first_approx);
	std::cout << getOptimalValue(canon_problem, solution);
	return 0;
}
