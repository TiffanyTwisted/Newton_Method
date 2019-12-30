// Советы по началу работы 
//   1. В окне обозревателя решений можно добавлять файлы и управлять ими.
//   2. В окне Team Explorer можно подключиться к системе управления версиями.
//   3. В окне "Выходные данные" можно просматривать выходные данные сборки и другие сообщения.
//   4. В окне "Список ошибок" можно просматривать ошибки.
//   5. Последовательно выберите пункты меню "Проект" > "Добавить новый элемент", чтобы создать файлы кода, или "Проект" > "Добавить существующий элемент", чтобы добавить в проект существующие файлы кода.
//   6. Чтобы снова открыть этот проект позже, выберите пункты меню "Файл" > "Открыть" > "Проект" и выберите SLN-файл.

#ifndef PCH_H
#define PCH_H

#include <iostream>
#include <iomanip>
const double M = 0.000001; // 0.01, 0.05 , 0.1
const int itr_max = 100;

using namespace std;

double** createMatrix(int);
void deleteMatrix(double **, int);
void viewMatrix(double **, int);
void viewAnswer(double *, int);
void copyMatrix(double **, double **, int);
bool Gauss(double *, double **, int);
double Newton(double , double , int );
double difff1x1(double );
double difff1x2(double );
double difff2x1(double);
double difff2x2(double , double );

double f1(double , double );
double f2(double , double );


#endif //PCH_H
