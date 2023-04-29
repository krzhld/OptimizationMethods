#include "task.h"

double Task::MinFunc(column_t x)
{
	return x[2];
}

double Task::Limit(column_t x)
{
	vector<double> limits;
	limits.push_back(3 * x[0] + 2 * x[1] - 4);
	limits.push_back(-10 * x[0] + 7 * x[1] - 16);
	limits.push_back((x[0] + 1) * (x[0] + 1) + (x[1] - 1) * (x[1] - 1) - 2.5);
	limits.push_back(Phi0(x) - x[2]);
	
	double result = limits[0];

	for (auto lim : limits)
	{
		if (result < lim)
			result = lim;
	}

	return result;
}

column_t Task::SubgradientLim(column_t x)
{
	int index = Index(x);
	int xs, ys, zs;
	if (index == 0)
	{
		xs = 3;
		ys = 2;
		zs = 0;
	}
	else if (index == 1)
	{
		xs = -10;
		ys = 7;
		zs = 0;
	}
	else if (index == 2)
	{
		xs = 2 * (x[0] + 1);
		ys = 2 * (x[1] - 1);
		zs = 0;
	}
	else if (index == 3)
	{
		xs = 2 * x[0] + 4 * cos(4 * x[0] + 5 * x[1]) + 3;
		ys = 10 * x[1] + 5 * cos(4 * x[0] + 5 * x[1]) + 2;
		zs = -1;
	}
	else
	{
		exit(-1);
	}
	column_t result;
	result.push_back(xs);
	result.push_back(ys);
	result.push_back(zs);

	return result;
}


double Task::Phi0(column_t x)
{
	return x[0] * x[0] + 5 * x[1] * x[1] + sin(4 * x[0] + 5 * x[1]) + 3 * x[0] + 2 * x[1];
}

int Task::Index(column_t x)
{
	vector<double> limits;
	limits.push_back(3 * x[0] + 2 * x[1] - 4);
	limits.push_back(-10 * x[0] + 7 * x[1] - 16);
	limits.push_back((x[0] + 1) * (x[0] + 1) + (x[1] - 1) * (x[1] - 1) - 2.5);
	limits.push_back(Phi0(x) - x[2]);

	int result = 0;
	int max = limits[0];
	for (int i = 1; i < size(limits); i++)
	{
		if (max < limits[i])
		{
			max = limits[i];
			result = i;
		}
			
	}

	return result; 
}

polyhedron_t Task::GetS0()
{
	matrix_t A; // матрица A, определ¤юща¤ многогранное множество
	A = { {3,2,0}, {-10,7,0}, {0,0,-1}, {0,-1,0}, {0,0,1} };

	column_t b; //столбец b, определ¤ющий многограное множество
	b = { -4, 16, 3.5, 0.7, 0 };

	return make_tuple(A, b);
}