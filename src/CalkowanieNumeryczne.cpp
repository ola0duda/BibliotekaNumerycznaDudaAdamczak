#include "CalkowanieNumeryczne.h"
#include "Interpolacja.h" // Do³¹czenie dla funkcji horner
#include <iostream>

using namespace std;

// Definicja zmiennej globalnej 'ai', aby biblioteka by³a samowystarczalna.
double ai[5] = { 0 };


// --- Ca³kowanie wielomianu zadanego przez wspó³czynniki ---

// Metoda prostok¹tów
double calka_prostokat(double* przedzial, double* ai_coeffs, int stopien, int n) {
	double h = (przedzial[1] - przedzial[0]) / n;
	double suma = 0.0;
	for (int i = 0; i < n; i++) {
		double xi = przedzial[0] + i * h;
		suma += horner(xi, ai_coeffs, stopien + 1);
	}
	return suma * h;
}

// Metoda trapezów
double calka_trapez(double* przedzial, double* ai_coeffs, int stopien, int n) {
	double h = (przedzial[1] - przedzial[0]) / n;
	double suma = horner(przedzial[0], ai_coeffs, stopien + 1) + horner(przedzial[1], ai_coeffs, stopien + 1);
	for (int i = 1; i < n; i++) {
		suma += 2 * horner(przedzial[0] + i * h, ai_coeffs, stopien + 1);
	}
	return suma * h / 2.0;
}

// Metoda Simpsona
double calka_simpson(double* przedzial, double* ai_coeffs, int stopien, int n) {
	if (n % 2 != 0) {
		// W metodzie Simpsona liczba podprzedzia³ów musi byæ parzysta.
		n++;
	}
	double h = (przedzial[1] - przedzial[0]) / n;
	double suma = horner(przedzial[0], ai_coeffs, stopien + 1) + horner(przedzial[1], ai_coeffs, stopien + 1);

	for (int i = 1; i < n; i += 2) {
		suma += 4 * horner(przedzial[0] + i * h, ai_coeffs, stopien + 1);
	}
	for (int i = 2; i < n - 1; i += 2) {
		suma += 2 * horner(przedzial[0] + i * h, ai_coeffs, stopien + 1);
	}
	return suma * h / 3.0;
}


// --- Ca³kowanie z u¿yciem wskaŸnika na funkcjê ---

// Metoda prostok¹tów
double calka_prostokat_funkcja(double* przedzial, double(*f)(double), int n) {
	double h = (przedzial[1] - przedzial[0]) / n;
	double suma = 0.0;
	for (int i = 0; i < n; i++) {
		suma += f(przedzial[0] + i * h);
	}
	return suma * h;
}

// Metoda trapezów
double calka_trapez_funkcja(double* przedzial, double(*f)(double), int n) {
	double h = (przedzial[1] - przedzial[0]) / n;
	double suma = f(przedzial[0]) + f(przedzial[1]);
	for (int i = 1; i < n; i++) {
		suma += 2 * f(przedzial[0] + i * h);
	}
	return suma * h / 2.0;
}

// Metoda Simpsona
double calka_simpson_funkcja(double* przedzial, double(*f)(double), int n) {
	if (n % 2 != 0) {
		n++;
	}
	double h = (przedzial[1] - przedzial[0]) / n;
	double suma = f(przedzial[0]) + f(przedzial[1]);

	for (int i = 1; i < n; i += 2) {
		suma += 4 * f(przedzial[0] + i * h);
	}
	for (int i = 2; i < n - 1; i += 2) {
		suma += 2 * f(przedzial[0] + i * h);
	}
	return suma * h / 3.0;
}


// --- Kwadratura Gaussa-Legendre'a ---

void GL(int n, double x[], double w[]) {
	if (n == 2) {
		x[0] = -1.0 / sqrt(3); x[1] = 1.0 / sqrt(3);
		w[0] = 1.0; w[1] = 1.0;
	}
	else if (n == 3) {
		x[0] = -sqrt(3.0 / 5.0); x[1] = 0.0; x[2] = sqrt(3.0 / 5.0);
		w[0] = 5.0 / 9.0; w[1] = 8.0 / 9.0; w[2] = 5.0 / 9.0;
	}
	else if (n == 4) {
		x[0] = -0.861136; x[1] = -0.339981; x[2] = 0.339981; x[3] = 0.861136;
		w[0] = 0.347855; w[1] = 0.652145; w[2] = 0.652145; w[3] = 0.347855;
	}
}

double kwadraturaGL1(double (*f)(double), double a, double b, int n, int m) {
	double x[5], w[5];
	GL(n, x, w);

	double calka = 0.0;
	double h = (b - a) / m;
	double c1 = (b - a) / (2.0 * m);

	for (int j = 0; j < m; ++j) {
		double srodek_podprzedzialu = a + (j + 0.5) * h;
		for (int i = 0; i < n; ++i) {
			calka += w[i] * f(srodek_podprzedzialu + c1 * x[i]);
		}
	}
	return calka * c1;
}

// Ta funkcja zale¿y od globalnej tablicy `ai`, która jest teraz zdefiniowana na górze pliku.
double wielomian(double x) {
	// U¿ywamy hornera z zdefiniowanym rozmiarem 5, zgodnym z globaln¹ tablic¹ `ai`
	return horner(x, ai, 5);
}