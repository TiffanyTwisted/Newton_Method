// Newton.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include "pch.h"

using namespace std;


double** createMatrix(int x)
{
	double **A = new double*[x];
	for (int i = 0; i < x; i++)
		A[i] = new double[x + 1];
	return A;
}

void deleteMatrix(double **X, int x)
{
	for (int i = 0; i < x; i++)
		delete X[i];
	delete[] X;
}

void copyMatrix(double **X, double **copyX, int x)
{
	for (int i = 0; i < x; i++)
	{
		for (int j = 0; j < x + 1; j++)
			copyX[i][j] = X[i][j];
	}
}

void viewMatrix(double **X, int x)
{
	for (int i = 0; i < x; i++)
	{
		for (int j = 0; j < x + 1; j++)
			cout << setw(12) << X[i][j];
		cout << endl;
	}
	cout << endl << endl;
}

void viewAnswer(double *X, int x)
{
	for (int i = 0; i < x; i++)
		cout << setw(12) << X[i];
	cout << endl;
}

bool Gauss(double *An, double **X, int x)
{
	for (int k = 0; k < x; k++)
	{
		double max = fabs(X[k][k]);
		int remeber = k;		//запоминаем строку, чтобы не поменять саму себя
		for (int i = k + 1; i < x; i++)
		{
			if (max < fabs(X[i][k]))		//находим максимальный по модулю элемент в столбце
			{
				max = fabs(X[i][k]);
				remeber = i;
			}
		}

		if (fabs(max - 0) < 1e-6)
		{
			return 0;
		}

		if (k != remeber)				//меняем строки местами
		{
			double *temp = X[k];
			X[k] = X[remeber];
			X[remeber] = temp;
		}

		//viewMatrix(X, x);

		double lead = X[k][k];			//запоминаем ведущий элемент
		for (int r = k; r < x + 1; r++)
		{
			X[k][r] /= lead;
		}
		//начиная со следующей строки приводим исходную матрицу к диагональному виду
		for (int i = k + 1; i < x; i++)
		{
			double temp = X[i][k];
			for (int j = k; j < x + 1; j++)
			{
				X[i][j] -= X[k][j] * temp;
			}
		}
		//viewMatrix(X, x);
	}

	An[x - 1] = X[x - 1][x + 1 - 1];				//обратный ход
	for (int i = x - 2; i >= 0; i--)
	{
		An[i] = X[i][x + 1 - 1];
		for (int j = i + 1; j < x + 1 - 1; j++)
		{
			An[i] -= X[i][j] * An[j];
		}
	}
	return 1;
}
//первое уравнение системы
double f1(double x1, double x2) {
   return 1.5*x1*x1*x1 - x2*x2 - 1;
}
//второе уравнение системы
double f2(double x1, double x2) {
   return x1*x2*x2*x2 - x2 - 4;
	
}



typedef double(*pf)(double, double); //pf указатель на функцию возвращающую double
//подсчет производной через дифференциал (как lim)
//Аргументы - (уравнение системы, переменные, по какой переменной ищем)
double Differential(pf f, double x1, double x2, int x)
{
	if (x == 1)											//по первой неизвестной
		return (f(x1 + M, x2) - f(x1, x2)) / M;
	else                                                //по второй неизвестной
		return (f(x1, x2 + M) - f(x1, x2)) / M;
}
double difff1x1( double x1) {   // производная f1 по x1
	return (4.5*x1*x1);
}
double difff1x2(double x2) {   // производная f1 по x2
	return (-2 * x2);
}
double difff2x1(double x2) {  // производная f2 по x1
	return (x2*x2*x2);
}
double difff2x2(double x1, double x2) {    // производная f2 по x2
	return (3 * x2 * x1 - 1);
}


//метод Ньютона
double Newton(double *function, double *approximate, int n)
{
	double **F = createMatrix(n);
	double *An = new double[n];
	double e = 1e-9, b1 = 1, b2 = 1;

	cout << endl << setw(15)  <<  "k" << setw(15)<< "b1      " << setw(15) << "b2      " << endl << endl;
	int itr = 1;

	do {
		//заполнение матрицы Якоби 
		// 1 способ
		
		  F[0][0] = Differential(f1, approximate[0], approximate[1], 1);
		  F[0][1] = Differential(f1, approximate[0], approximate[1], 2);
		  F[0][2] = -f1(approximate[0], approximate[1]);					
		  F[1][0] = Differential(f2, approximate[0], approximate[1], 1);
	      F[1][1] = Differential(f2, approximate[0], approximate[1], 2);
	 	  F[1][2] = -f2(approximate[0], approximate[1]);					
		  

		// 2 способ 
		/*
		 F[0][0] = difff1x1( approximate[0]) ;
		 F[0][0] = difff1x1( approximate[0]) ;
         F[0][1] = difff1x2( approximate[1]) ;
         F[0][2] = -f1(approximate[0], approximate[1]);	;
         F[1][0] = difff2x1( approximate[1]) ;
         F[1][1] = difff2x2( approximate[0], approximate[1]) ;
         F[1][2] = -f2(approximate[0], approximate[1]);	
		 
		 */


		double **copyF = createMatrix(n);					// копия для сохранения оригинала
		copyMatrix(F, copyF, n);
		//viewMatrix(copyF, n);

		if (!Gauss(An, copyF, n))				//Проверка метода Гаусса
		{
			deleteMatrix(copyF, n);		//очистка памяти
			deleteMatrix(F, n);
			delete An;
			//	cout << itr;
			return 0;
		}

		approximate[0] += An[0];				//уточнение решения
		approximate[1] += An[1];

		if (fabs(f1(approximate[0], approximate[1])) > fabs(f2(approximate[0], approximate[1]))) 
			b1 = fabs(f1(approximate[0], approximate[1]));
		else
			b1 = fabs(f2(approximate[0], approximate[1]));

		for (int i = 0; i < n; i++)									
		{
			if (fabs(An[i]) < 1)
				b2 = fabs(An[i]);
			else if (fabs(An[i]) >= 1)
				b2 = fabs(An[i] / approximate[i]);
		}

		cout << setw(15) << itr << setw(15) << b1 << setw(15) << b2 << endl;
		itr++;
	} while ((b1 > e || b2 > e) && (itr < itr_max));

	function[n - 2] = approximate[0];
	function[n - 1] = approximate[1];

	//deleteMatrix(copyF, n);
	deleteMatrix(F, n);				//очистка памяти
	delete An;

	return *function;	//возвращаем решение СНЛАУ
}

int main()
{
	const int n = 2;
	double *Approximate = new double[n];
	double *Function = new double[n];
	Approximate[0] = 1; Approximate[1] = 1;
	*Function = Newton(Function, Approximate, n);
	if (*Function != 0)
	{
		cout << "\n\n Answer: ";
		viewAnswer(Function, n);
	}
	else
		cout << "\n\n There is no solution!\n\n";

	
}

