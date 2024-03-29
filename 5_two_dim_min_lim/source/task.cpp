#include "task.h"

#define PI 3.1416

double Task::MinFunc(column_t x) {
	return x[2];
}

double Task::Limit(column_t x) {
	vector<double> limits;

	limits.push_back(4 * x[0] + 5 * x[1] + 7);
	limits.push_back(-4 * x[0] - 5 * x[1] - 9);
	limits.push_back(2 * (x[0] + 2) * (x[0] + 2) - x[1] - 2);
	limits.push_back(Phi0(x) - x[2]);
	
	double result = limits[0];

	for (auto lim : limits) {
		if (result < lim)
			result = lim;
	}

	return result;
}

column_t Task::SubgradientLim(column_t x) {
	int index = Index(x);
	double xs, ys, zs;
	if (index == 0) {
		xs = 4;
		ys = 5;
		zs = 0;
	}
	else if (index == 1) {
		xs = -4;
		ys = -5;
		zs = 0;
	}
	else if (index == 2) {
		xs = 4 * (x[0] + 2);
		ys = -1;
		zs = 0;
	}
	else if (index == 3) {
		xs = 2 * x[0] + 4 * cos(4 * x[0] + 5 * x[1]) + 3;
		ys = 10 * x[1] + 5 * cos(4 * x[0] + 5 * x[1]) + 2;
		zs = -1;
	}
	else
		exit(-1);

	column_t result;
	result.push_back(xs);
	result.push_back(ys);
	result.push_back(zs);

	return result;
}


double Task::Phi0(column_t x) {
	return x[0] * x[0] + 5 * x[1] * x[1] + sin(4 * x[0] + 5 * x[1]) + 3 * x[0] + 2 * x[1];
}

int Task::Index(column_t x) {
	vector<double> limits;
	limits.push_back(4 * x[0] + 5 * x[1] + 7);
	limits.push_back(-4 * x[0] - 5 * x[1] - 9);
	limits.push_back(2 * (x[0] + 2) * (x[0] + 2) - x[1] - 2);
	limits.push_back(Phi0(x) - x[2]);

	int result = 0;
	double max = limits[0];
	for (int i = 1; i < size(limits); i++) {
		if (max < limits[i]) {
			max = limits[i];
			result = i;
		}
			
	}

	return result; 
}

polyhedron_t Task::GetS0() {
	matrix_t A; // матрица A, определ¤юща¤ многогранное множество
	A = { 
		{4, 5, 0},
		{-4, - 5, 0},
		{1, -1, 0},
		{-1, 1, 0},
		{0, 0, 1},
		{0, 0, -1}
	};

	column_t b; //столбец b, определ¤ющий многограное множество
	b = {
		-2 * PI,
		3 * PI,
		0,
		6,
		100,
		100
	};

	return make_tuple(A, b);
}
