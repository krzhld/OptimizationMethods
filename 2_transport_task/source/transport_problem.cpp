#include "transport_problem.h"

using namespace std;

transport_problem_t ReadFromFileTransportProblem(string filename)
{
	// ��������� ����� ������ �� �����
	ifstream input_file_stream(filename);

	string cur_string;

	int N = 0; // ���������� �����������
	int M = 0; // ���������� �����������
	double cur_number = 0; // �������� �����
	column_t a; // ������ ���������� ������ � �����������
	column_t b; // ������ ������� ������ � �����������
	matrix_t c; // ������� ���������� ��������� �� ������� ����������� � ������ �����������

	getline(input_file_stream, cur_string);

	// ���� ���� �� ���������� � a: �� ��������� ���������
	if (cur_string != "a:") {
		cout << "Error reading file!" << endl;
		exit(-1);
	}
	else {
		getline(input_file_stream, cur_string);
		// ��������� ����� ������ �� ����������� ������
		stringstream str_stream(cur_string);

		// ��������� ���������� ������ � �����������
		while (str_stream >> cur_number)
			a.push_back(cur_number);

		//������ ���������� �����������
		N = a.size();
	}

	getline(input_file_stream, cur_string);
	// ���� � ����� ��� b: �� ��������� ���������
	if (cur_string != "b:") {
		cout << "Error reading file!" << endl;
		exit(-1);
	}
	else {
		getline(input_file_stream, cur_string);
		// ��������� ����� ������ �� ����������� ������
		stringstream str_stream(cur_string);

		// ��������� ����������� ������ ������������
		while (str_stream >> cur_number)
			b.push_back(cur_number);

		//������ ���������� �����������
		M = b.size();

		//������ ������ ������� ��������� ���������
		c.resize(N); 
	}
	getline(input_file_stream, cur_string);
	if (cur_string != "c:") {
		cout << "Error reading file!" << endl;
		exit(-1);
	}
	else {
		for (int i = 0; i < N; i++)
		{
			getline(input_file_stream, cur_string);
			stringstream str_stream(cur_string);
			//c[i].resize(N);
			while (str_stream >> cur_number)
			   c[i].push_back(cur_number);
			
		}
		
	}
	// ��������� ����� �����
	input_file_stream.close();
	transport_problem_t problem = make_tuple(a, b, c);
	return problem;

}

/*����� ������-��������� ���� ������ ���������� �����������*/
void NorthWestCornerMethod(column_t a, column_t b, int N, int M, matrix_t &X)
{
	int i = 0, j = 0;
	while ((i < N) || (j < M))
	{
		if (a[i] <= b[j]) //���� ���������� ������ � i-�� ���������� ������ ������� j-�� �����������
		{
			X[i][j] = a[i]; //����� ��������� �� i � j ����� ����� ������, ������� ������� � i-�� ����������
			if ((i == N - 1) && (j == M - 1))
			{
				return;
			}
			b[j] -= a[i];  // ������ j-�� ����������� ����������� �� �����, ��������� � i-�� ����������
			a[i] = 0; //�����, ��������� � i-�� ���������� ��������� ��������
			for (int k = j+1; k < M; k++)
			{
				X[i][k] = -1; // ���������� ������ i-�� ������ ������������ ������� ������ (����������� -1)
			}
			i++; // ������� �� ��������� ������ (� ���������� ����������)
		
		}

		else // ��� �� �� �����, ������ � ������, ���� ������ j-�� ����������� ������ ������, ���������� � i-�� ����������
		{
			X[i][j] = b[j];
			if ((i == N - 1) && (j == M - 1))
			{
				return;
			}
			a[i] -= b[j];
			b[j] = 0;
			for (int k = i+1; k < N; k++)
			{
				X[k][j] = -1;
			}
			j++;
		
		}


	}
}

/*����� �����������*/
void FindPotentials(matrix_t &c, int N, int M, matrix_t &X, column_t& u, column_t& v)
{
	u[0] = 0; // ������ �������� ������ ��������� ������ ���� 

	vector<bool> u1, v1; // ������� �����������: ��������� �� ������ �������� �����������
	u1.resize(N);
	v1.resize(M);

	u1[0] = true;
	for (int i = 1; i < N; i++) 
	{
		u1[i] = false;
	}
	for (int i = 0; i < M; i++)
	{
		v1[i] = false;
	}

	for (int k = 0; k < M + N; k++) //���� ������ ���� �����������
	{
		for (int i = 0; i < N; i++) //���� ������ ����������� �����������
		{
			if (u1[i]) // ���� i-�� ��������� ���������� ������ 
			{
				for (int j = 0; j < M; j++) //����� ����������� ����������, ������� ���������� ���� �� j-�� ������
				{
					if (X[i][j] != -1)
					{
						if (!v1[j])
						{
							v[j] = c[i][j] -u[i];
							v1[j] = true;
						}
					}
				}
			}
		}


		for (int l = 0; l < M; l++) //���� ������ ����������� �����������
		{
			if (v1[l]) // ���� l-�� ��������� ���������� ������ 
			{
				for (int k = 0; k < N; k++) //����� ����������� �����������, ������� ���������� ���� � l-�� �����
				{
					if (X[k][l] != -1)
					{
						if (!u1[k])
						{
							u[k] = c[k][l] - v[l];
							u1[k] = true;
						}
					}
				}
			}
		}

	}

}

/*����� ������������ ����� ������ ��� �������� ������������� �������*/
min_delta_t FindMinDelta(matrix_t &c, column_t &u, column_t &v, int N, int M, matrix_t &X)
{
	double mindelta = 10e20;
	int mini = -1, minj = -1;

	column_t colDelta;
	colDelta.resize(M);
	matrix_t delta; //������ ����� ������
	delta.resize(N);
	for (int i = 0; i < N; i++)
	{
		delta[i] = colDelta;
	}

	//���� ��������� ���� ������
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < M; j++)
		{
			if (X[i][j] == -1)
			{
				delta[i][j] =  c[i][j] - (v[j] + u[i]);
				if (delta[i][j] < mindelta)
				{
					mindelta = delta[i][j]; 
					mini = i;
					minj = j;
				}
					
			}
		}
	}
	
	min_delta_t minDelta = make_tuple(mini, minj, mindelta);
	return minDelta;
}

/*����� ��������� ���������*/
double FindCost(matrix_t &c, matrix_t &X, int N, int M)
{
	double result = 0;
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < M; j++)
		{
			result += c[i][j] * X[i][j];
		}
	}
	return result;
}

/*����� ����� � ����� ������������ �� ������*/
double FindCeelsLine(int inext, int jnext, int imin, int jmin, int N, int M, matrix_t &X, matrix_t& temp, int numb, double teta)
{
	double result = -1;
	//���� ������ ���� ����� ������� ������
	for (int j = 0; j < M; j++)
	{
		//����� ����������� ������ (�������� ������� �������) ��� ��������� ������
		if (((X[inext][j] != -1) && (j != jnext)) || ((j == jmin) && (inext == imin) && (numb != 0)))
		{
			numb++; //����� ������ � �����
			double tetaold = -1; //�������� ���������� ��� ������������ �������� � ������ ������� ������
			//���� ������ ����� �������� ����� � �����
			if (numb % 2)
			{
				tetaold = teta; //�������� ����������� �������� � ����� �� ������ ������� ������� ������
				if (teta < 0)
				{
					teta = X[inext][j]; //���� ��� ������ ��������� ����������� ������, ������ �� �� ������� 
				}
				//���� ��� ������ ��������� � �� �������� ������ �������� �������� - ��������� �� �������
				else if ((X[inext][j] < teta) && (X[inext][j]!=-1)) 
				{
					teta = X[inext][j];
				}
			}
			//���� ������� ���������� �� ��������� ������ � ���������� ������ � ��� �����
			if ((j == jmin) && (inext == imin) && (numb % 2 == 0))
			{
				temp[imin][jmin] = teta; //����������� ��������� ������ �������� ����������� ��������
				return teta;
			}
			//���� ������� �� ����������, ������ �������, ��������� � ������ �� �������
			else
			{
				result = FindCeelsColumn(inext, j, imin, jmin, N, M, X, temp, numb, teta);
			}
			

			//���� ���� ���������, ��������� ��������������� �������, �������� ��� ��������
			if (result != -1)
			{
				//���� ����� ������, �� ����� ����������� �������� �� ������ ����
				if (numb % 2 == 0)
				{
					temp[inext][j] = temp[imin][jmin];
				}
				// ���� ��������, �� �� ������ �����
				else
				{
					temp[inext][j] = -temp[imin][jmin];
				}
				break;
			}
			// ���� ���� �� ��������� - �������, ������������ �����
			else
			{
				numb--;
				if (tetaold >= 0)
				{
					teta = tetaold;
				}
			}
		}

		
	}
	return result;
}

/*����� ����� � ����� ������������ �� �������*/
double FindCeelsColumn(int inext, int jnext, int imin, int jmin, int N, int M, matrix_t &X, matrix_t& temp, int numb, double teta)
{
	double result = -1;
	//���� ������ ���� ����� �������� �������
	for (int i = 0; i < N; i++)
	{
		//����� ����������� ������ (�������� ������� �������) ��� ��������� ������
		if (((X[i][jnext] != -1) && (i != inext)) || ((i == imin) && (jnext == jmin) && (numb != 0)))
		{
			numb++; //����� ������ � �����
			double tetaold = -1; //�������� ���������� ��� ������������ �������� � ������ ������� ������
			//���� ������ ����� �������� ����� � �����
			if (numb % 2)
			{
				tetaold = teta; //�������� ����������� �������� � ����� �� ������ ������� ������� ������
				if (teta < 0)
				{
					teta = X[i][jnext]; //���� ��� ������ ��������� ����������� ������, ������ �� �� ������� 
				}
				//���� ��� ������ ��������� � �� �������� ������ �������� �������� - ��������� �� �������
				else if ((X[i][jnext] < teta) && (X[i][jnext] != -1))
				{
					teta = X[i][jnext];
				}
			}
			//���� ������� ���������� �� ��������� ������ � ���������� ������ � ��� �����
			if ((i == imin) && (jnext == jmin) && (numb % 2 == 0))
			{
				temp[imin][jmin] = teta; //����������� ��������� ������ �������� ����������� ��������
				return teta;
			}
			//���� ������� �� ����������, ������ �������, ��������� � ������ �� ������
			else
			{
				result = FindCeelsLine(i, jnext, imin, jmin, N, M, X, temp, numb, teta);
			}


			//���� ���� ���������, ��������� ��������������� �������, �������� ��� ��������
			if (result != -1)
			{
				//���� ����� ������, �� ����� ����������� �������� �� ������ ����
				if (numb % 2 == 0)
				{
					temp[i][jnext] = temp[imin][jmin];
				}
				// ���� ��������, �� �� ������ �����
				else
				{
					temp[i][jnext] = -temp[imin][jmin];
				}
				break;
			}
			// ���� ���� �� ��������� - �������, ������������ �����
			else
			{
				numb--;
				if (tetaold >= 0)
				{
					teta = tetaold;
				}
			}
		}

	}
	return result;
}

void NewX(matrix_t &X, matrix_t &temp, int N, int M, int imin, int jmin)
{
	int flag = 0; //���������, ������� �� ���������� �� ������
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < M; j++)
		{
			if ((temp[i][j] != -10e10) && ((i != imin)||(j != jmin)))
			{
				X[i][j] += temp[i][j]; //������ �������� �������� ���������� � ������������ �� ��������� ����, ���������� ��� ��������� �����
				if ((temp[i][j] <= 0) && (X[i][j] == 0) && (flag == 0))
				{
					X[i][j] = -1; //������� ���������� �� ������
					flag = 1; //������ �������� ����������
				}
			}
		}
	}
	X[imin][jmin] += 1; //���������� ������� � ����� �������� ����������
	X[imin][jmin] = temp[imin][jmin]; //������ �������� ����� �������� ����������
}

/*����� ����������� ������� ������������ ������ ��������� ����*/
solving_t MethodOfPotentials(transport_problem_t problem)
{
	column_t a, b; 
	matrix_t c;
	tie(a, b, c) = problem;
	int N = size(a);
	int M = size(b);

	column_t columnX;
	columnX.resize(M);
	matrix_t X;
	X.resize(N);
	for (int i = 0; i < N; i++)
	{
		X[i] = columnX;
	}

	NorthWestCornerMethod(a, b, N, M, X); // ���� ������ ����������� ������� ������-��������� ���� 

	column_t v, u;
	v.resize(M); // ���������� �����������
	u.resize(N); // ���������� �����������

	FindPotentials(c, N, M, X, u, v); //���� ���������� 

	int imin, jmin;
	double mindelta;


	tie(imin, jmin, mindelta) = FindMinDelta(c, u, v, N, M, X);


	while (mindelta < 0) //�������� ���� ������, ��������� ����������, ����� ���������� ����������� ������� (��� ������������� ������)
	{
		
		matrix_t tmp; //��������������� ������� ��� ������������ ����� ���������
		column_t tmpColumn;
		tmpColumn.resize(M);
		tmp.resize(N);
		for (int i = 0; i < N; i++)
		{
			tmp[i] = tmpColumn;
		}

		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < M; j++)
			{
				tmp[i][j] = -10e10;
			}
		}


		//������������� ����, ����� ����� �������� ���������� � ��������� imin, jmin
		double teta = FindCeelsLine(imin, jmin, imin, jmin, N, M, X, tmp, 0, -1);
		

		NewX(X, tmp, N, M, imin, jmin);


		v.clear();
		u.clear();

		//column_t v, u;
		v.resize(M); // ���������� �����������
		u.resize(N); // ���������� �����������

		FindPotentials(c, N, M, X, u, v); //���� ���������� 


		tie(imin, jmin, mindelta) = FindMinDelta(c, u, v, N, M, X);

	}
	
	for (int i = 0; i < N; i++) //������ ������ -1 ���������� � ����������� ������� ������� 0
	{
		for (int j = 0; j < M; j++)
		{
			if (X[i][j] == -1)
			{
				X[i][j] = 0;
			}
			
		}
	}

	double resultCost = FindCost(c, X, N, M); //������� ����������� ���������
	solving_t result = make_tuple(X, resultCost);
	return result;
}