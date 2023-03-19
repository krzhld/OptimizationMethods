#include "transport_problem.h"

using namespace std;

transport_problem_t ReadFromFileTransportProblem(string filename)
{
	// открываем поток чтения из файла
	ifstream input_file_stream(filename);

	string cur_string;

	int N = 0; // количество поставщиков
	int M = 0; // количество покупателей
	double cur_number = 0; // буферное число
	column_t a; // вектор количества товара у поставщиков
	column_t b; // вектор запроса товара у покупателей
	matrix_t c; // матрица стоимостей перевозок из пунктов поставщиков в пункты покупателей

	getline(input_file_stream, cur_string);

	// если файл не начинается с a: то завершаем программу
	if (cur_string != "a:") {
		cout << "Error reading file!" << endl;
		exit(-1);
	}
	else {
		getline(input_file_stream, cur_string);
		// открываем поток чтения из прочитанной строки
		stringstream str_stream(cur_string);

		// считываем количество товара у поставщиков
		while (str_stream >> cur_number)
			a.push_back(cur_number);

		//задаем количество поставщиков
		N = a.size();
	}

	getline(input_file_stream, cur_string);
	// если в файле нет b: то завершаем программу
	if (cur_string != "b:") {
		cout << "Error reading file!" << endl;
		exit(-1);
	}
	else {
		getline(input_file_stream, cur_string);
		// открываем поток чтения из прочитанной строки
		stringstream str_stream(cur_string);

		// считываем потребность товара покупателями
		while (str_stream >> cur_number)
			b.push_back(cur_number);

		//задаем количество покупателей
		M = b.size();

		//задаем размер матрицы стоимости перевозок
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
	// закрываем поток файла
	input_file_stream.close();
	transport_problem_t problem = make_tuple(a, b, c);
	return problem;

}

/*Получение из транспортной задачи канонической задачи в классической постановке*/
canon_problem_t GetCanonProblemFromTransportProblem(transport_problem_t& transportProblem)
{
	matrix_t tC;
	column_t ta, tb; // данные транспортной задачи
	tie(ta, tb, tC) = transportProblem;

	int n = size(ta);
	int m = size(tb);

	matrix_t cA; // матрица А канонической ЗЛП
	cA.resize(n * m); //в транспортной задаче n*m переменных
	column_t cAcolumn; 
	cAcolumn.resize(n + m - 1); //в канонической задаче будет n+m-1 уравнений

	for (int i = 0; i < n * m; i++)
	{
		cA[i] = cAcolumn;
	}

	//заполнение первых n строк матрицы 
	int j = 0;
	for (int i = 0; i < n; i++) 
	{
		while (j < n*m)
		{
			cA[j][i] = 1;
			j++;
			if (j % m == 0)
				break;
		}
	}

	//заполнение последних строк матрицы
	int k = 0;
	for (int i = n; i < n + m - 1; i++)
	{
		int l = 0;
		while (l < n * m)
		{
			if ((l-k) % m == 0)
				cA[l][i] = 1;
			l++;
		}
		k++;
	}

	column_t cb; //столбец b канонической ЗЛП
	cb.resize(n + m - 1);

	//запонение первых n элементов столбца b
	for (int i = 0; i < n; i++)
	{
		cb[i] = ta[i];
	}

	// запонление остальных элементов столбца b
	k = 0;
	for (int i = n; i < n + m - 1; i++)
	{
		cb[i] = tb[k];
		k++;
	}

	column_t cc; //целевая функция канонической ЗЛП
	cc.resize(n * m);
	k = 0;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			cc[k] = tC[i][j];
			k++;
		}
	
	}
	canon_problem_t canonProblem = make_tuple(cA, cb, cc);
	return canonProblem;

}

double SumColumn(column_t col)
{
	double res = 0;
	for (auto& c : col)
	{
		res += c;
	}
	return res;
}

double SumMatrix(matrix_t matrix, int N, int M)
{
	double res=0;
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < M; j++)
		{
			res += matrix[i][j];
		}
	}
	return res;
}

solving_t SolveTransportProblem(transport_problem_t problem)
{
	column_t a, b;
	matrix_t c;
	tie(a, b, c) = problem;
	int N = size(a);
	int M = size(b);

	solving_t solving; 

	double sumA = SumColumn(a); //суммарное количество товара у поставщиков
	double sumB = SumColumn(b); //суммарная потребность всех потребителей

	if (sumA == sumB) //если задача сбалансирована
	{
		solving = MethodOfPotentials(problem);
		return solving;
	}
	if (sumA < sumB) //если запасы поставщиков меньше потребностей покупателей
	{
		double difference = sumB - sumA; //разность суммарной потребности и суммарных запасов

		double meanCost = SumMatrix(c, N, M) / (N * M); //средняя стоимость перевозок 

		double fine = double(int(meanCost) * 5); //штрафы за недопоставку 
		
	
		a.push_back(difference); // создаем фиктивного поставщика, запасы которого равны разности потребностей покупателей и суммарных запасов поставщиков
		c.resize(N + 1);
		c[N].resize(M); 
		for (int i = 0; i < M;i++) 
		{
			c[N][i] = fine; //стоимости перевозок от фиктивных поставщиков равны штрафу за недопоставку
		}

		transport_problem_t new_problem = make_tuple(a, b, c);

		matrix_t X;
		double result;
		tie(X, result) = MethodOfPotentials(new_problem); //решение полученной сбалансированной задачи методом потенциалов 
		X.resize(N);
		solving = make_tuple(X, result);
		// суммарная стоимость перевозок выводится с учетом штрафов
		return solving;
	}

}

/*Метод северо-западного угла поиска начального приближения*/
void NorthWestCornerMethod(column_t a, column_t b, int N, int M, matrix_t &X)
{
	int i = 0, j = 0;
	while ((i < N) || (j < M))
	{
		if (a[i] <= b[j]) //если количество товара у i-го поставщика меньше запроса j-го потребителя
		{
			X[i][j] = a[i]; //объем перевозки из i в j пункт равен объему, который имеется у i-го поставщика
			if ((i == N - 1) && (j == M - 1))
			{
				return;
			}
			b[j] -= a[i];  // запрос j-го потребителя уменьшается на объем, имеющийся у i-го поставщика
			a[i] = 0; //объем, имеющийся у i-го поставщика полностью исчерпан
			for (int k = j+1; k < M; k++)
			{
				X[i][k] = -1; // оставщиеся клетки i-ой строки транспортной таблицы пустые (заполняются -1)
			}
			i++; // переход на следующую строку (к следующему поставщику)
		
		}

		else // все то же самое, только в случае, если запрос j-го потребителя меньше объема, имеющегося у i-го поставщика
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


/*Поиск потенциалов*/
void FindPotentials(matrix_t &c, int N, int M, matrix_t &X, column_t& u, column_t& v)
{
	u[0] = 0; // обычно полагают данный потенциал равным нулю 

	vector<bool> u1, v1; // векторы индикаторов: вычислены ли данные значения потенциалов
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

	for (int k = 0; k < M + N; k++) //цикл поиска всех потенциалов
	{
		for (int i = 0; i < N; i++) //цикл поиска потенциалов покупателей
		{
			if (u1[i]) // если i-ый потенциал поставщика найден 
			{
				for (int j = 0; j < M; j++) //поиск потенциалов покупателя, которым поставляют груз из j-го пункта
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


		for (int l = 0; l < M; l++) //цикл поиска потенциалов поставщиков
		{
			if (v1[l]) // если l-ый потенциал покупателя найден 
			{
				for (int k = 0; k < N; k++) //поиск потенциалов поставщиков, которые поставляют груз в l-ый пункт
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

/*Поиск минимального числа дельта для проверки оптимальности решения*/
min_delta_t FindMinDelta(matrix_t &c, column_t &u, column_t &v, int N, int M, matrix_t &X)
{
	double mindelta = 10e20;
	int mini = -1, minj = -1;

	column_t colDelta;
	colDelta.resize(M);
	matrix_t delta; //массив чисел дельта
	delta.resize(N);
	for (int i = 0; i < N; i++)
	{
		delta[i] = colDelta;
	}

	//цикл пересчета всех дельта
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

/*Общая стоимость перевозок*/
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

/*Поиск ячеек в цикле перестроения по строке*/
double FindCeelsLine(int inext, int jnext, int imin, int jmin, int N, int M, matrix_t &X, matrix_t& temp, int numb, double teta)
{
	double result = -1;
	//цикл обхода хода вдоль текущей строки
	for (int j = 0; j < M; j++)
	{
		//поиск заполненной ячейки (исключая текущую позицию) или стартовой ячейки
		if (((X[inext][j] != -1) && (j != jnext)) || ((j == jmin) && (inext == imin) && (numb != 0)))
		{
			numb++; //номер ячейки в цикле
			double tetaold = -1; //запасная переменная для минимального значения в случае неудачи обхода
			//если ячейка имеет нечетный номер в цикле
			if (numb % 2)
			{
				tetaold = teta; //значение минимальной поставки в цикле на случай неудачи данного обхода
				if (teta < 0)
				{
					teta = X[inext][j]; //если это первая встречная заполненная клетка, примем ее за минимум 
				}
				//если эта клетка заполнена и ее значение меньше текущего минимума - принимаем за минимум
				else if ((X[inext][j] < teta) && (X[inext][j]!=-1)) 
				{
					teta = X[inext][j];
				}
			}
			//если ломаная замкнулась на начальной клетке и количество клеток в ней четно
			if ((j == jmin) && (inext == imin) && (numb % 2 == 0))
			{
				temp[imin][jmin] = teta; //присваиваем стартовой клетке значение минимальной поставки
				return teta;
			}
			//если ломаная не замкнулась, делаем поворот, переходим к поиску по столбцу
			else
			{
				result = FindCeelsColumn(inext, j, imin, jmin, N, M, X, temp, numb, teta);
			}
			

			//если цикл замкнулся, заполняем вспомогательную матрицу, обратный ход рекурсии
			if (result != -1)
			{
				//если номер четный, то берем минимальную поставку со знаком плюс
				if (numb % 2 == 0)
				{
					temp[inext][j] = temp[imin][jmin];
				}
				// если нечетный, то со знаком минус
				else
				{
					temp[inext][j] = -temp[imin][jmin];
				}
				break;
			}
			// если цикл не замкнулся - неудача, откатываемся назад
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

/*Поиск ячеек в цикле перестроения по столбцу*/
double FindCeelsColumn(int inext, int jnext, int imin, int jmin, int N, int M, matrix_t &X, matrix_t& temp, int numb, double teta)
{
	double result = -1;
	//цикл обхода хода вдоль текущего столбца
	for (int i = 0; i < N; i++)
	{
		//поиск заполненной ячейки (исключая текущую позицию) или стартовой ячейки
		if (((X[i][jnext] != -1) && (i != inext)) || ((i == imin) && (jnext == jmin) && (numb != 0)))
		{
			numb++; //номер ячейки в цикле
			double tetaold = -1; //запасная переменная для минимального значения в случае неудачи обхода
			//если ячейка имеет нечетный номер в цикле
			if (numb % 2)
			{
				tetaold = teta; //значение минимальной поставки в цикле на случай неудачи данного обхода
				if (teta < 0)
				{
					teta = X[i][jnext]; //если это первая встречная заполненная клетка, примем ее за минимум 
				}
				//если эта клетка заполнена и ее значение меньше текущего минимума - принимаем за минимум
				else if ((X[i][jnext] < teta) && (X[i][jnext] != -1))
				{
					teta = X[i][jnext];
				}
			}
			//если ломаная замкнулась на начальной клетке и количество клеток в ней четно
			if ((i == imin) && (jnext == jmin) && (numb % 2 == 0))
			{
				temp[imin][jmin] = teta; //присваиваем стартовой клетке значение минимальной поставки
				return teta;
			}
			//если ломаная не замкнулась, делаем поворот, переходим к поиску по строке
			else
			{
				result = FindCeelsLine(i, jnext, imin, jmin, N, M, X, temp, numb, teta);
			}


			//если цикл замкнулся, заполняем вспомогательную матрицу, обратный ход рекурсии
			if (result != -1)
			{
				//если номер четный, то берем минимальную поставку со знаком плюс
				if (numb % 2 == 0)
				{
					temp[i][jnext] = temp[imin][jmin];
				}
				// если нечетный, то со знаком минус
				else
				{
					temp[i][jnext] = -temp[imin][jmin];
				}
				break;
			}
			// если цикл не замкнулся - неудача, откатываемся назад
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
	int flag = 0; //индикатор, удалена ли переменная из базиса
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < M; j++)
		{
			if ((temp[i][j] != -10e10) && ((i != imin)||(j != jmin)))
			{
				X[i][j] += temp[i][j]; //меняем значения базисных переменных в соответствии со значением тета, полученным при пересчете цикла
				if ((temp[i][j] <= 0) && (X[i][j] == 0) && (flag == 0))
				{
					X[i][j] = -1; //удаляем переменную из базиса
					flag = 1; //меняем значение индикатора
				}
			}
		}
	}
	X[imin][jmin] += 1; //прибавляем единицу к новой базисной переменной
	X[imin][jmin] = temp[imin][jmin]; //задаем значение новой базисной переменной
}

/*Метод потенциалов решения транспортной задачи закрытого типа*/
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

	NorthWestCornerMethod(a, b, N, M, X); // ищем первое приближение методом северо-западного угла 

	column_t v, u;
	v.resize(M); // потенциалы покупателей
	u.resize(N); // потенциалы поставщиков

	FindPotentials(c, N, M, X, u, v); //ищем потенциалы 

	int imin, jmin;
	double mindelta;


	tie(imin, jmin, mindelta) = FindMinDelta(c, u, v, N, M, X);


	while (mindelta < 0) //основной цикл метода, остановка просиходит, когда достигнуто оптимальное решение (нет отрицательных дельта)
	{
		
		matrix_t tmp; //вспомогательная матрица для перестроения цикла пересчета
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


		//перестраиваем цикл, вводя новую базисную переменную с индексами imin, jmin
		double teta = FindCeelsLine(imin, jmin, imin, jmin, N, M, X, tmp, 0, -1);
		

		NewX(X, tmp, N, M, imin, jmin);


		v.clear();
		u.clear();

		//column_t v, u;
		v.resize(M); // потенциалы покупателей
		u.resize(N); // потенциалы поставщиков

		FindPotentials(c, N, M, X, u, v); //ищем потенциалы 


		tie(imin, jmin, mindelta) = FindMinDelta(c, u, v, N, M, X);

	}
	
	for (int i = 0; i < N; i++) //делаем равные -1 переменные в оптимальном решении равными 0
	{
		for (int j = 0; j < M; j++)
		{
			if (X[i][j] == -1)
			{
				X[i][j] = 0;
			}
			
		}
	}


	double resultCost = FindCost(c, X, N, M); //находим оптимальную стоимость
	solving_t result = make_tuple(X, resultCost);
	return result;
}