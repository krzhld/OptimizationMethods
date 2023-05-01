#include "task.h"
#include "cutting_hyperplane_method.h"

int main() {
	Task t;
	setlocale(LC_ALL, "Russian");
	double eps = 1e-3;
	column_t x = CuttingHyperplaneMethod(t, eps);
	cout << "����������� �����: (" << x[0] << ", " << x[1] << ")" << endl;
	cout << "����������� ��������: " << x[2];
	return 0;
}
