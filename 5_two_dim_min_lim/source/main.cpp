#include "task.h"
#include "cutting_hyperplane_method.h"

int main() {
	Task t;
	setlocale(LC_ALL, "Russian");
	column_t x;
	double res;
	tie(x, res) = CuttingHyperplaneMethod(t, 0.00001);
	cout << "����������� �����: (" << x[0] << ", " << x[1] << ")" << endl;
	cout << "����������� ��������: " << x[2];
	return 0;
}
