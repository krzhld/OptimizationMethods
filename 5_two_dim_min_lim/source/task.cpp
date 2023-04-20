#include "task.h"

double Task::MinFunc(double x, double y, double z)
{
	return z;
}

double Task::Limit(double x, double y, double z)
{
	vector<double> limits;
	limits.push_back(3 * x + 2 * y - 4);
	limits.push_back(-10 * x + 7 * y - 16);
	limits.push_back((x + 1) * (x + 1) + (y - 1) * (y - 1) - 2.5);
	limits.push_back(Phi0(x, y) - z);
	
	double result = limits[0];

	for (auto lim : limits)
	{
		if (result < lim)
			result = lim;
	}

	return result;
}

vector<double> Task::SubgradientLim(double x, double y, double z)
{
	int index = Index(x, y, z);
	int xs, ys, zs;
	if (index == 0)
	{
		// прописать градиент
	}
	else if (index == 1)
	{
		// прописать градиент
	}
	else if (index == 2)
	{
		// прописать градиент
	}
	else if (index == 3)
	{
		// прописать градиент
	}
	else
	{
		exit(-1);
	}
	vector<double> result;
	result.push_back(xs);
	result.push_back(ys);
	result.push_back(zs);

	return result;
}


double Task::Phi0(double x, double y)
{
	return x * x + 5 * y * y + sin(4 * x + 5 * y) + 3 * x + 2 * y;
}

int Task::Index(double x, double y, double z)
{
	vector<double> limits;
	limits.push_back(3 * x + 2 * y - 4);
	limits.push_back(-10 * x + 7 * y - 16);
	limits.push_back((x + 1) * (x + 1) + (y - 1) * (y - 1) - 2.5);
	limits.push_back(Phi0(x, y) - z);

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
