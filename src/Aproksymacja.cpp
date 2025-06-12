#include "Aproksymacja.h"
#include "RozwiazywanieUkladowRownanLiniowych.h" // Potrzebne dla funkcji trojkat
#include <iostream>
#include <iomanip>
#include <functional>

using namespace std;



// Funkcja bazowa (potêga)
double baza(double x, int n) {
	return pow(x, n);
}

// Kwadratura Gaussa-Legendre'a (4 wêz³y) do obliczania ca³ek
double kwadraturaGL2(std::function<double(double)> f, double a, double b) {
	double x[4], w[4];
	// Wspó³rzêdne i wagi w kwadraturze G-L dla 4 wêz³ów
	x[0] = -0.861136;
	x[1] = -0.339981;
	x[2] = 0.339981;
	x[3] = 0.861136;
	w[0] = 0.347855;
	w[1] = 0.652145;
	w[2] = 0.652145;
	w[3] = 0.347855;

	double wynik = 0.0;
	double h = (b - a) / 10000;

	for (int j = 0; j < 10000; ++j) {
		double aj = a + j * h;
		double bj = aj + h;

		for (int i = 0; i < 4; ++i) {
			double xi = ((bj - aj) / 2) * x[i] + (aj + bj) / 2;
			wynik += w[i] * f(xi);
		}
	}

	wynik *= (b - a) / (2 * 10000);

	return wynik;
}

// Eliminacja Gaussa z podstawianiem wstecznym
void gauss2(double** A, double* b, double* x, int n) {
	for (int i = 0; i < n; i++) {
		int max = i;
		for (int j = i + 1; j < n; j++) {
			if (fabs(A[j][i]) > fabs(A[max][i])) {
				max = j;
			}
		}
		swap(A[i], A[max]);
		swap(b[i], b[max]);

		if (fabs(A[i][i]) < 1e-9) {
			cout << "Otrzymano wiersz zerowy - problem z rozwiazaniem ukladu rownan" << endl;
			return;
		}

		double pivot = A[i][i];
		for (int j = i; j < n; j++) {
			A[i][j] /= pivot;
		}
		b[i] /= pivot;

		for (int j = i + 1; j < n; j++) {
			double wspolczynnik = A[j][i] / A[i][i];
			for (int k = i; k < n; k++) {
				A[j][k] -= wspolczynnik * A[i][k];
			}
			b[j] -= wspolczynnik * b[i];
		}
	}
	trojkat(A, b, x, n);
}

// G³ówna funkcja do aproksymacji
void aproksymacja(double a, double b, int stopien, double (*funkcja)(double)) {
	double** A = new double* [stopien];
	double* B = new double[stopien];
	double* x = new double[stopien];
	for (int i = 0; i < stopien; i++) {
		A[i] = new double[stopien];
	}

	for (int i = 0; i < stopien; i++) {
		B[i] = kwadraturaGL2([=](double x) {return funkcja(x) * baza(x, i); }, a, b);
		for (int j = 0; j < stopien; j++) {
			A[i][j] = kwadraturaGL2([=](double x) { return baza(x, i) * baza(x, j); }, a, b);
		}
	}

	gauss2(A, B, x, stopien);

	cout << "Wspolczynniki wielomianu aproksymacyjnego:" << endl;
	for (int i = 0; i < stopien; i++) {
		cout << "a" << i << " = " << x[i] << endl;
	}
	cout << fixed << setprecision(6);
	cout << endl;
	cout << "x\tF(x)\t\tP(x)\t\tBlad: " << endl;
	for (double xx = a; xx <= b; xx += 0.1) {
		double fx = funkcja(xx);
		double px = 0;
		for (int i = 0; i < stopien; i++) {
			px += x[i] * baza(xx, i);
		}
		cout << xx << "\t" << fx << "\t" << px << "\t" << fabs(fx - px) << endl;
	}
	for (int i = 0; i < stopien; i++) {
		delete[] A[i];
	}
	delete[] A;
	delete[] B;
	delete[] x;
}