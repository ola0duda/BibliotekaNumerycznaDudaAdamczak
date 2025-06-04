#include "RozwiazywanieUkladowRownanLiniowych.h"

#include <iostream>
#include <iomanip>

using namespace std;

//wypisywanie macierzy
void pokazMacierz1(double** A, double* b, int n) {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			cout << setw(8) << A[i][j] << " ";
		}
		cout << "| " << setw(8) << b[i] << endl;
	}
	cout << endl;
}

//algorytm eliminacji Gaussa uwzgl dniaj cy przypadki: liniowej zaleznosci wierszy, cz ciowy pivoting (redukcja wierszy)
void gauss(double** A, double* b, int n) {
	cout << "Eliminacja Gaussa" << endl;

	for (int i = 0; i < n; i++) {
		//szukanie wiersza z najwiekszym elementem w kolumnie i zamiana
		int max = i;
		for (int j = i + 1; j < n; j++) {
			if (fabs(A[j][i]) > fabs(A[max][i])) {
				max = j;
			}
		}
		swap(A[i], A[max]);
		swap(b[i], b[max]);

		//sprawdzenie, czy nie otrzymany zostanie wiersz zerowy
		if (fabs(A[i][i]) < 1e-9) {					//<1e-9 a nie =0 poniewa  nie ma mo liwo ci por wnania warto ci do warto ci bezwzgl dnej fabs()
			cout << "Otrzymano wiersz zerowy - problem z rozwiazaniem ukladu rownan" << endl;
			return;
		}

		//normalizacja wiersza - jedynki na przek tnej
		double pivot = A[i][i];
		for (int j = i; j < n; j++) {
			A[i][j] /= pivot;
		}
		b[i] /= pivot;

		//zerowanie elementow pod przekatna
		for (int j = i + 1; j < n; j++) {
			double wspolczynnik = A[j][i] / A[i][i];
			for (int k = i; k < n; k++) {
				A[j][k] -= wspolczynnik * A[i][k];
			}
			b[j] -= wspolczynnik * b[i];
		}

		//wypisanie macierzy po ka dym kroku eliminacji
		cout << "Po " << i + 1 << " kroku algorytmu eliminacji Gaussa" << endl;
		pokazMacierz1(A, b, n);
	}
}

//podstawienie wsteczne
void trojkat(double** A, double* b, double* x, int n) {
	for (int i = n - 1; i >= 0; i--) {
		x[i] = b[i];
		for (int j = i + 1; j < n; j++) {
			x[i] -= A[i][j] * x[j];
		}
		x[i] /= A[i][i];
	}
}

//metoda LU-------------------------------------------
void pokazMacierz2(double** A, int n) {			//wypisywanie macierzy
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			cout << setw(8) << A[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl;
}

void rozkladLU(double** A, double** L, double** U, double* b, int n) {				//nowy argument funkcji - wektor wyrazów wolnych b, poniewa¿ w pivotingu konsekwentnie musz¹ byæ zamieniane tak¿e wyrazy wolne
	for (int i = 0; i < n; i++) {
		//pivoting----------------------------
		int max = i;			//szukanie wiersza z najwiêkszym elementem w kolumnie i
		for (int j = i + 1; j < n; j++) {
			if (fabs(A[j][i]) > fabs(A[max][i])) {			//porównanie wartoœci bezwzglêdnych w kolumnie i
				max = j;
			}
		}

		//zamiana wiersza i z wierszem max
		swap(A[i], A[max]);
		swap(b[i], b[max]);

		//zamiana tak¿e odpowiednich, ju¿ obliczonych elementów w macierzy L (poniewa¿ zmienione zosta³o przypisanie wierszy w A, a L zale¿y od poprzednich kroków opartych na A)
		for (int k = 0; k < i; k++) {
			swap(L[i][k], L[max][k]);
		}
		//koniec pivotingu---------------------

		for (int j = i; j < n; j++) {
			U[i][j] = A[i][j];
			for (int k = 0; k < i; k++) {
				U[i][j] -= L[i][k] * U[k][j];
			}
		}

		//sprawdzenie, czy element na przek¹tnej nie jest zbyt bliski 0
		if (fabs(U[i][i]) <= 1e-12) {
			cout << "Blad: macierz osobliwa, brak rozwiazan." << endl;
			cout << "U[" << i << "][" << i << "] = " << U[i][i] << " po kroku " << i << ", stan macierzy A w momencie skonczenia obliczen:" << endl;
			pokazMacierz2(A, n);
			return;
		}

		for (int j = i; j < n; j++) {
			if (i == j) {
				L[i][i] = 1;
			}
			else {
				L[j][i] = A[j][i];
				for (int k = 0; k < i; k++) {
					L[j][i] -= L[j][k] * U[k][i];
				}
				if (fabs(U[i][i]) <= 1e-12) {
					cout << "Blad: dzielenie przez zero" << endl;
				}
				L[j][i] /= U[i][i];
			}
		}
		cout << "Macierz L po " << i + 1 << " iteracji:" << endl;
		pokazMacierz2(L, n);
		cout << "Macierz U po " << i + 1 << " iteracji:" << endl;
		pokazMacierz2(U, n);
	}
}

void rozwiazanieLz(double** L, double* b, double* z, int n) {			//rozwiazanie Lz=b (podstawienie w przód)
	for (int i = 0; i < n; i++) {
		z[i] = b[i];
		for (int j = 0; j < i; j++) {
			z[i] -= L[i][j] * z[j];
		}
	}
}

void rozwiazanieUx(double** U, double* z, double* x, int n) {			//rozwiazanie Ux=z (podstawienie w ty³)
	for (int i = n - 1; i >= 0; i--) {
		x[i] = z[i];
		for (int j = i + 1; j < n; j++) {
			x[i] -= U[i][j] * x[j];
		}
		x[i] /= U[i][i];
	}
}
