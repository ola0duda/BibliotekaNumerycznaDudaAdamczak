#include "CalkowanieNumeryczne.h"

#include <iostream>

using namespace std;

//do obliczania normalnie z liczb

//metoda prostok¹tów
double calka_prostokat(double* przedzial, double* ai, int stopien, int n) {
	double xi, fx;
	double h = (przedzial[1] - przedzial[0]) / n;
	double suma = 0.0;
	double wynik;

	for (int i = 0; i < n; i++) {
		xi = przedzial[0] + i * h;
		fx = horner(xi, ai, stopien + 1);
		suma += fx;
	}

	wynik = suma * h;

	return wynik;
}

// wzór trapezów
double calka_trapez(double* przedzial, double* ai, int stopien, int n) {
	double xi, fxi, xi1, fxi1;
	double h = (przedzial[1] - przedzial[0]) / n;
	double iloczyn = 0.0, wynik = 0.0;

	xi1 = przedzial[0];
	fxi1 = horner(xi1, ai, stopien + 1);

	for (int i = 1; i <= n; i++) {
		xi = przedzial[0] + i * h;
		fxi = horner(xi, ai, stopien + 1);
		iloczyn += (xi - xi1) * (fxi + fxi1);			//do pola trapezu (pó¿niej dzielenie przez 2)

		xi1 = xi;			//przesuniecie xi na xi-1
		fxi1 = fxi;			//przesuniecie f(xi) na f(x-1)
	}
	wynik = iloczyn / 2;
	return wynik;
}

//wzór Simpsona
double calka_simpson(double* przedzial, double* ai, int stopien, int n) {
	//n to parzysta liczba przedzia³ów (podwójne kroki, 2n punktów)
	double h = (przedzial[1] - przedzial[0]) / (2 * n);
	double suma1 = 0.0, suma2 = 0.0;		//sumy dla f(x2i-1) i f(x2i)
	double wynik = 0.0;

	double fx0 = horner(przedzial[0], ai, stopien + 1);
	double fx2n = horner(przedzial[1], ai, stopien + 1);

	for (int i = 1; i < 2 * n; i++) {
		double xi = przedzial[0] + i * h;
		double fxi = horner(xi, ai, stopien + 1);
		if (i % 2 == 0) {
			suma2 += fxi;
		}
		else {
			suma1 += fxi;
		}
	}

	wynik = (h / 3.0) * (fx0 + 4 * suma1 + 2 * suma2 + fx2n);
	return wynik;
}
//------------------
//z funkcja jako argumentem
	//metoda prostok¹tów
double calka_prostokat_funkcja(double* przedzial, double(*f)(double), int n) {
	double xi;
	double h = (przedzial[1] - przedzial[0]) / n;
	double suma = 0.0;
	double wynik;

	for (int i = 0; i < n; i++) {
		xi = przedzial[0] + i * h;
		suma += f(xi);
	}

	wynik = suma * h;

	return wynik;
}

//wzór trapezów
double calka_trapez_funkcja(double* przedzial, double(*f)(double), int n) {
	double xi, xi1;
	double h = (przedzial[1] - przedzial[0]) / n;
	double iloczyn = 0.0, wynik = 0.0;

	xi1 = przedzial[0];

	for (int i = 1; i <= n; i++) {
		xi = przedzial[0] + i * h;
		iloczyn += (xi - xi1) * (f(xi) + f(xi1));			//do pola trapezu (pó¿niej dzielenie przez 2)

		xi1 = xi;			//przesuniecie xi na xi-1
	}
	wynik = iloczyn / 2;
	return wynik;
}

//wzór Simpsona
double calka_simpson_funkcja(double* przedzial, double(*f)(double), int n) {
	//n to parzysta liczba przedzia³ów (podwójne kroki, 2n punktów)
	double h = (przedzial[1] - przedzial[0]) / (2 * n);
	double suma1 = 0.0, suma2 = 0.0;		//sumy dla f(x2i-1) i f(x2i)
	double wynik = 0.0;

	double fx0 = f(przedzial[0]);
	double fx2n = f(przedzial[1]);

	for (int i = 1; i < 2 * n; i++) {
		double xi = przedzial[0] + i * h;
		if (i % 2 == 0) {
			suma2 += f(xi);
		}
		else {
			suma1 += f(xi);
		}
	}

	wynik = (h / 3.0) * (fx0 + 4 * suma1 + 2 * suma2 + fx2n);
	return wynik;
}

//wspolrzedne i wagi w kwadraturze G-L (2, 3 i 4 wêze³)
void GL(int n, double x[], double w[]) {
	if (n == 2) {
		x[0] = -1.0 / sqrt(3);
		x[1] = 1.0 / sqrt(3);
		w[0] = 1.0;
		w[1] = 1.0;
	}
	else if (n == 3) {
		x[0] = -sqrt(3.0 / 5.0);
		x[1] = 0.0;
		x[2] = sqrt(3.0 / 5.0);
		w[0] = 5.0 / 9.0;
		w[1] = 8.0 / 9.0;
		w[2] = 5.0 / 9.0;
	}
	else if (n == 4) {
		x[0] = -0.861136;
		x[1] = -0.339981;
		x[2] = 0.339981;
		x[3] = 0.861136;
		w[0] = 0.347855;
		w[1] = 0.652145;
		w[2] = 0.652145;
		w[3] = 0.347855;
	}
}

//funkcja obliczaj¹ca kwadraturê - jednym z jej argumentów jest funkcja podca³kowa
double kwadraturaGL1(double (*f)(double), double a, double b, int n, int m) {                //m - liczba podprzedzialow
	double x[5], w[5];
	GL(n, x, w);

	double wynik = 0.0;
	double h = (b - a) / m;

	for (int j = 0; j < m; ++j) {
		double aj = a + j * h;
		double bj = aj + h;

		for (int i = 0; i < n; ++i) {
			double xi = ((bj - aj) / 2) * x[i] + (aj + bj) / 2;         //dostosowanie wartosci pod przedzial calkowania (wagi podane w innej funkcji s¹ odpowiednie dla [-1, 1])
			wynik += w[i] * f(xi);           //wynik obliczany zgodnie z funkcj¹ wi*f(xi)
		}
	}

	wynik *= (b - a) / (2 * m);

	return wynik;
}

double ai[];		//zmienna globalna ze wspolczynnikami wielomianu, ¿eby wykorzystaæ je w funkcji podca³kowej wielomian()

double wielomian(double x) {
	return ai[0] + ai[1] * x + ai[2] * pow(x, 2) + ai[3] * pow(x, 3) + ai[4] * pow(x, 4);
}
