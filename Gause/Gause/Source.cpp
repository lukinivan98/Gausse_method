#include <iostream>
#include <cmath>
#include <iomanip>

using namespace std;

double f(double x) 
{
	return exp(x);
}

void get_moments(double* mom)
{
	mom[0] = 2.0 / 3.0;
	mom[1] = 4.0 / 15.0;
	mom[2] = 16.0 / 105.0;
	mom[3] = 32.0 / 315.0;
}

void GetPoly(double coef[], double &a1, double & a2)
{
	double main_det = coef[1] * coef[1] - coef[2] * coef[0];
	double det1 = coef[0] * coef[3] - coef[1] * coef[2];
	double det2 = coef[2] * coef[2] - coef[1] * coef[3];

	a1 = det1 / main_det;
	a2 = det2 / main_det;
}

void FindRoots(double a1, double a2, double &x1, double &x2)
{
	double a = 1;
	double b = a1;
	double c = a2;
	double d = b * b - 4 * a*c;
	x1 = 0.5 * (-b + sqrt(d)) / a;
	x2 = 0.5 * (-b - sqrt(d)) / a;
}

void CalculateA(double x1, double x2, double* moment, double &A1, double &A2)
{
	A1 = (moment[1] - x2 * moment[0]) / (x1 - x2);
	A2 = (moment[1] - x1 * moment[0]) / (x2 - x1);
}

double CalculateGause(double A1, double A2, double x1, double x2)
{
	return A1 * f(x1) + A2 * f(x2);
}

double move(double t, double a, double b)
{
	return 0.5*(b - a)* t + 0.5*(b + a);
}

double FullGause(double coef[])
{
	double a1, a2, x1, x2, A1, A2;
	GetPoly(coef, a1, a2);
	FindRoots(a1, a2, x1, x2);
	CalculateA(x1, x2, coef, A1, A2);
	return CalculateGause(A1, A2, x1, x2);
}

int main()
{
	double a = 0;
	double b = 1;
	int N = 2;
	double rez = 1.003008;
	double moments[4];
	get_moments(moments);
	double G = FullGause(moments);

	cout << "int (e^x * sqrt(1-x) * dx) from 0 to 1 = " << rez << "\n";
	cout << "Gause formule =  \n" << G << "\n";
	cout << "error = " << abs(rez - G) << "\n" << "\n";

	

	
	double x0 = 0;
	double x1 = 1.0 / 2.0;
	double x2 = 1;
	
	double A0 = 0.171429;
	double A1 = 0.45714;
	double A2 = 0.038;

	double Interp;
	Interp = A0 * f(x0) + A1 * f(x1) + A2 * f(x2);


	cout << "int (e^x * sqrt(1-x) * dx) from 0 to 1 = " << rez << "\n";
	cout << "Interpol formule =  \n" << Interp << "\n";
	cout << "error = " << abs(rez - Interp) << "\n";
	system("pause");
	return 0;
}