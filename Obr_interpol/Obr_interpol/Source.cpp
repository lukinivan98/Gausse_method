#include <iostream>
#include <cmath>
#include <iomanip>

using namespace std;

const double Separation_eps = 0.0005;
const double Newton_eps = 0.0000001;
const double der_eps = 0.000001;

double Left[20]; /// предполагаемые промежутки, где могут находиться корни
double Right[20];
int counter = 0; /// счетчик корней

double *X;
double *Y;
int n;
double h = 0.1;

double NewtonRoots[20];
int Ncounter = 0;

double lagrange1(double* x, double* y, double _x)
{
	double result = 0.0;

	for (int i = 0; i < n; i++)
	{
		double P = 1.0;

		for (int j = 0; j < n; j++)
			if (j != i)
				P *= (_x - x[j]) / (x[i] - x[j]);

		result += P * y[i];
	}

	return result;
}

double f(double x)
{
	return lagrange1(X, Y, x) - 0.8;
}

double fder(double x)
{
	return (f(x + der_eps) - f(x)) / der_eps;
}

double find_max_fder(double a, double b)
{
	double eps = 0.005;
	double M = fder(a);
	double f;
	a += eps;
	while (a <= b)
	{
		f = fder(a);
		if (f > M) M = f;
		a += eps;
	}
	return M;
}

void Separation(double a0, double b0)
{
	double a, b = b0;

	double fleft, fright;

	while (b >= a0 + Separation_eps)
	{
		a = b - Separation_eps;
		fleft = f(a);
		fright = f(b);
		if (fleft * fright < 0)
		{
			Left[counter] = a;
			Right[counter] = b;
			counter++;
		}
		b -= Separation_eps;
	}
}

double NewtonMethod(double a, double b)
{
	double x1 = (a + b) / 2, x0 = b;
	double delta;

	while (abs(x1 - x0) > Newton_eps)
	{
		x0 = x1;
		delta = -f(x0) / fder(x0);
		x1 = x0 + delta;
	}
	return x1;
}

void NewtonMethod()
{
	for (int i = 0; i < counter; i++)
		NewtonRoots[Ncounter++] = NewtonMethod(Left[i], Right[i]);
}

void printNewton()
{
	cout << "Корни по Методу Ньютона:\n";
	for (int i = 0; i < Ncounter; i++)
		cout << setprecision(10) << NewtonRoots[i] << "\n";
}

void print_intervals()
{
	cout << "Интервалы, где есть ровно один корень:\n";
	for (int i = 0; i < counter; i++)
		cout << setprecision(10) << Left[i] << " " << Right[i] << "\n";
}

int main()
{

	n = 11;

	X = new double[n];
	Y = new double[n];

	X[0] = 0.0;
	Y[0] = 0.0;
	for (int i = 1; i < n; i++)
	{
		X[i] = X[i - 1] + h;
		Y[i] = asin(X[i]);
	}

	setlocale(LC_CTYPE, "Russian");

	cout << "Обратная интерполяция многочленом Лагранжа в точке 0.8: " << lagrange1(Y, X, 0.8) << "\n";

	Separation(-2, 2);
	print_intervals();

	NewtonMethod();
	printNewton();

	system("pause");
    return 0;
}
