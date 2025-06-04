#include "Interpolacja.h"

#include <iostream>

//interpolacja lagrange'a--------------------------------------
double lagrange(double punkt, double rozmiar, double* tX, double* tY, int odlegloscWezly) {
	double L = 0;
	for (int i = 0; i < rozmiar; i += odlegloscWezly) {
		double licznik = 1, mianownik = 1;
		for (int j = 0; j < rozmiar; j += odlegloscWezly) {
			if (j != i) {
				licznik *= (punkt - tX[j]);
				mianownik *= (tX[i] - tX[j]);
			}
		}
		L += tY[i] * (licznik / mianownik);			//wielomian interpolacyjny Lagrange'a
	}
	return L;
}

//obliczanie wartoœci wielomianu zadanego w postaci naturalnej dla wartoœci ?? wczytanej z pliku
double naturalna(double wartosc, double* wspolczynniki, int rozmiar) {
	double w = wspolczynniki[0];
	for (int i = 1; i < rozmiar; i++) {
		w += wspolczynniki[i] * pow(wartosc, i);
	}
	return w;
}

//obliczanie wartoœci wielomianu wed³ug schematu Hornera dla wartoœci ?? wczytanej z pliku
double horner(double wartosc, double* wspolczynniki, int rozmiar) {
	double w = wspolczynniki[rozmiar - 1];		//wn=an
	for (int i = rozmiar - 2; i >= 0; i--) {
		w = wspolczynniki[i] + wartosc * w;
	}
	return w;
}

//interpolacja newtona-------------------
//obliczanie ilorazów ró¿nicowych potrzebnych do wyznaczenia wspó³czynników ai - utworzenie tabeli, gdzie dwie pierwsze kolumny to wartoœci wêz³ów interpolacji i znane wartoœci funkcji
void ilorazy(double** f, double* x, double* y, int n) {
	for (int i = 0; i < n; i++) {
		f[i][0] = y[i];
	}
	for (int j = 1; j < n; j++) {
		for (int i = j; i < n; i++) {
			f[i][j] = (f[i][j - 1] - f[i - 1][j - 1]) / (x[i] - x[i - j]);
		}
	}
}

//obliczanie wartosci interpolacji funkcji wielomianem w postaci Newtona
double newton(double wartosc, double* a, double* x, int rozmiar) {
	double w = a[0];
	double mnozenie = 1.0;
	for (int i = 1; i < rozmiar; i++) {             //W(x)=a0+a1(x-x0)+a2(x-x1)(x-x0)+...
		mnozenie *= (wartosc - x[i - 1]);
		w += a[i] * mnozenie;
	}
	return w;
}
