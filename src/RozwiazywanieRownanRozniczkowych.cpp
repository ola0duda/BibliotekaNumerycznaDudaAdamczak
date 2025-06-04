#include "RozwiazywanieRownanRozniczkowych.h"

#include <iostream>

using namespace std;

double euler(int N, double T0_val, double t_tab[], double T_tab[], int max_rozmiar_tab, double (*rownanie)(double), double CZAS) {			//+zad.2 - obliczenia dla trzech ró¿nych kroków ca³kowania: h=1s, 10s, 100s
	double h = CZAS / N;			//N to liczba krokow
	T_tab[0] = T0_val;
	t_tab[0] = 0.0;

	for (int i = 0; i < N; ++i) {
		t_tab[i + 1] = (static_cast<double>(i) + 1.0) * h;
		T_tab[i + 1] = T_tab[i] + h * rownanie(T_tab[i]);

		if (T_tab[i + 1] < 0) T_tab[i + 1] = 0;
	}

	//upewnienie sie, ze ostatni czas to czas docelowy
	if (N > 0) {
		t_tab[N] = CZAS;
	}
	return N + 1;
}

int heun(int N, double T0, double t[], double T[], double (*rownanie)(double), double CZAS) {
	double h = static_cast<double>(CZAS) / N;
	t[0] = 0.0;
	T[0] = T0;

	for (int i = 0; i < N; ++i) {
		double k1 = rownanie(T[i]);
		double y = T[i] + h * k1;
		double k2 = rownanie(y);
		T[i + 1] = T[i] + (h / 2.0) * (k1 + k2);
		t[i + 1] = t[i] + h;
		if (T[i + 1] < 0) {
			T[i + 1] = 0;
		}
	}
	return N + 1;
}

int midpoint(int N, double T0, double t[], double T[], double (*rownanie)(double), double CZAS) {
	double h = static_cast<double>(CZAS) / N;
	t[0] = 0.0;
	T[0] = T0;

	for (int i = 0; i < N; ++i) {
		double k1 = rownanie(T[i]);
		double y = T[i] + (h / 2.0) * k1;
		double k2 = rownanie(y);
		T[i + 1] = T[i] + h * k2;
		t[i + 1] = t[i] + h;
		if (T[i + 1] < 0) {
			T[i + 1] = 0;
		}
	}
	return N + 1;
}

int runge_kutta(int N, double T0, double t[], double T[], double (*rownanie)(double), double CZAS) {
	double h = static_cast<double>(CZAS) / N;
	t[0] = 0.0;
	T[0] = T0;

	for (int i = 0; i < N; ++i) {
		double k1 = h * rownanie(T[i]);
		double k2 = h * rownanie(T[i] + k1 / 2.0);
		double k3 = h * rownanie(T[i] + k2 / 2.0);
		double k4 = h * rownanie(T[i] + k3);
		T[i + 1] = T[i] + (1.0 / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4);
		t[i + 1] = t[i] + h;
		if (T[i + 1] < 0) {
			T[i + 1] = 0;
		}
	}
	return N + 1;
}

