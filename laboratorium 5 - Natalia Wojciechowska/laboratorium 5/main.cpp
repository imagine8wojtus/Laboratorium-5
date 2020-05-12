#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "rk4.h"

#define g 9.81  //przypisanie wartosci dla g

void rhs_fun(double t, double* X, double* F);
void veuler(double t, double* X, double h, int n, void (* fun)(double, double*,double*), double* X1);

void main()
{
	FILE* f = fopen("Dane.csv", "w");   //wpisanie wartosci do pliku
	
	double h, tk, t, alfa0, omega0, E, m, l;
	int n = 2;
	h = 0.1;
	l = 3.;
	m = 4.;
	tk = 10.;
	alfa0 = 0.1;
	omega0 = 0.;
	double* X = (double*)malloc(n*sizeof(double));    //tablica wartosci zmiennych zaleznych omega i alfa w kroku t (dla metody Eulera)
	double* X1 = (double*)malloc(n*sizeof(double));   //tablica wartosci zmiennych zaleznych omega i alfa w kroku t+h (dla metody Eulera)
	double* Y = (double*)malloc(n*sizeof(double));    //tablica wartosci zmiennych zaleznych omega i alfa w kroku t (dla metody RK4)
	double* Y1 = (double*)malloc(n*sizeof(double));   //tablica wartosci zmiennych zaleznych omega i alfa w kroku t+h (dla metody RK4)
	X[0] = omega0;
	X[1] = alfa0;
	Y[0] = omega0;
	Y[1] = alfa0;

		for (t = 0; t <= tk; t+=h) {
		vrk4(t, Y, h, n, rhs_fun, Y1);       //funkcja metody rungego-kutty
		veuler(t, X, h, n, rhs_fun, X1);
		E = (m * l * l) / 2 * Y[0] * Y[0] + m * g * l * (1 - cos(Y[1]));     //obliczanie energii calkowitej wahadla
		Y[0] = Y1[0];                                            //przypisanie nowych wartosci 
		Y[1] = Y1[1];
		fprintf(f, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", t, X1[0], X1[1], Y1[0], Y1[1], E); 
																						
	}


	free(X);
	free(X1);
	free(Y);
	free(Y1);
	fclose(f);
}

void rhs_fun(double t, double* X, double* F) {   //funkcja do obliczania wartoœci prawych stron równañ ró¿niczkowych
	double l;
	l = 1.0;
	F[0] = (-1.0) * g / l * sin(X[1]);   //omega
	F[1] = X[0];						 //alfa
}

void veuler(double t, double* X, double h, int n, void (* fun)(double, double*, double*), double* X1) {
	double* F;
	F = (double*)malloc(n * sizeof(double)); 
	rhs_fun(t, X, F);
	X1[0] = X[0] + (h * F[0]);    //wzor dla kolejnych wartosci omega
	X1[1] = X[1] + (h * F[1]);    //wzor dla kolejnych wartosci alfa

	X[0] = X1[0];
	X[1] = X1[1];

	free(F);
}