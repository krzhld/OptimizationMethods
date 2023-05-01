#include "utils.h"


double MultiplipliedVectors(column_t v1, column_t v2) {
	if (v1.size() != v2.size())
		exit(-1);
	double result = 0;
	for (int i = 0; i < v1.size(); i++)
		result += v1[i] * v2[i];

	return result;
}

double Norm(column_t vector) {
	double result = 0;
	for (auto v : vector)
		result += v * v;

	return sqrt(result);
}

column_t SumVector(column_t v1, column_t v2) {
	if (v1.size() != v2.size())
		exit(-1);
	column_t result;
	result.resize(v1.size());
	for (int i = 0; i < result.size(); i++)
		result[i] = v1[i] + v2[i];

	return result;
}

column_t DiffVector(column_t v1, column_t v2) {
	if (v1.size() != v2.size())
		exit(-1);
	column_t result;
	result.resize(v1.size());
	for (int i = 0; i < result.size(); i++)
		result[i] = v1[i] - v2[i];

	return result;
}

column_t MultVectorAndNumber(column_t vector, double number) {
	for (auto v : vector)
		v *= number;

	return vector;
}
