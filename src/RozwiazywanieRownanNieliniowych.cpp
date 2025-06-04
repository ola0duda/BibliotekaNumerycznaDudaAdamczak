#include "RozwiazywanieRownanNieliniowych.h"

#include <iostream>
#include <iomanip>

using namespace std;

const double EPSILON = 1e-6;        //epsilon to wartosc wyznaczajaca do ktorego momentu program ma wykonywaæ metodê na przedziale

//metoda bisekcji
Wyniki bisekcja(double (*f)(double), double a, double b, double epsilon) {
	Wyniki wynik;
	wynik.iter_count = 0;       //poczatkowo zerowa liczba iteracji
	wynik.wynik = NAN;

	double fa = f(a);
	double fb = f(b);

	if (isnan(fa) || isnan(fb) || fa * fb >= 0) {
		return wynik;
	}

	double c = a;
	for (int i = 0; i < MAX_ITER; ++i) {
		c = (a + b) / 2.0;

		//zapisywanie iteracji
		if (wynik.iter_count < MAX_ITER) {
			wynik.iteracje[wynik.iter_count++] = c;
		}
		else {
			break; //osi¹gniêty limit iteracji w tablicy
		}

		double fc = f(c);
		if (isnan(fc)) {
			return wynik;
		}
		if (fabs(fc) < epsilon || (b - a) / 2.0 < epsilon) {
			wynik.wynik = c;
			return wynik;
		}
		if (fa * fc < 0) {
			b = c;
		}
		else {
			a = c;
			fa = fc;
		}
	}

	wynik.wynik = c;
	return wynik;
}

//metoda Newtona
Wyniki newton(double (*f)(double), double (*df)(double), double x0_start, double epsilon) {
	Wyniki wynik;
	wynik.iter_count = 0;
	wynik.wynik = NAN;

	double x0 = x0_start;
	if (isnan(x0)) {
		return wynik;
	}

	if (wynik.iter_count < MAX_ITER) {
		wynik.iteracje[wynik.iter_count++] = x0;
	}
	else { //MAX_ITER == 0
		double fx_check = f(x0); //sprawdzenie, czy x0 jest ju¿ rozwi¹zaniem
		wynik.wynik = (!isnan(fx_check) && fabs(fx_check) < epsilon) ? x0 : NAN;
		return wynik;
	}

	double x1;
	for (int iter = 0; iter < MAX_ITER; ++iter) {       //rowniez sprawdzenie czy liczba iteracji nie przekracza maksymalnej ilosci
		double fx = f(x0);
		double dfx = df(x0);
		if (isnan(fx) || isnan(dfx) || isinf(fx) || isinf(dfx)) {
			return wynik;
		}
		if (fabs(dfx) < 1e-12) { //dzielenie przez bardzo ma³¹ liczbê lub zero
			return wynik;
		}
		x1 = x0 - fx / dfx;     //zgodnie z twierdzeniem Taylora
		if (isnan(x1)) { //jeœli wynik operacji jest NAN
			return wynik;
		}

		if (wynik.iter_count < MAX_ITER) {
			wynik.iteracje[wynik.iter_count++] = x1;
		}
		else {
			break;
		}

		if (fabs(x1 - x0) < epsilon) {      //warunek konczacy obliczenia
			wynik.wynik = x1;
			return wynik;
		}

		x0 = x1;
	}
	wynik.wynik = x0;
	return wynik;
}

//metoda siecznych
Wyniki sieczna(double (*f)(double), double x0_start, double x1_start, double epsilon) {
	Wyniki wynik;
	wynik.iter_count = 0;
	wynik.wynik = NAN;

	double x0 = x0_start;
	double x1 = x1_start;
	if (isnan(x0) || isnan(x1)) {
		return wynik;
	}

	if (wynik.iter_count < MAX_ITER) {
		wynik.iteracje[wynik.iter_count++] = x0;
	}
	else {
		return wynik;
	}
	if (wynik.iter_count < MAX_ITER) {
		wynik.iteracje[wynik.iter_count++] = x1;
	}
	else {
		return wynik;
	}

	double x2;
	for (int iter = 0; iter < MAX_ITER; ++iter) {
		double f0 = f(x0);
		double f1 = f(x1);

		if (isnan(f0) || isnan(f1) || isinf(f0) || isinf(f1)) {
			return wynik;
		}
		if (fabs(f1 - f0) < 1e-12) {
			return wynik;
		}
		x2 = x1 - f1 * (x1 - x0) / (f1 - f0);        //rownanie metod siecznych
		if (isnan(x2)) {
			return wynik;
		}

		if (wynik.iter_count < MAX_ITER) {
			wynik.iteracje[wynik.iter_count++] = x2;
		}
		else {
			break;
		}
		if (fabs(x2 - x1) < epsilon) {
			wynik.wynik = x2;
			return wynik;
		}

		x0 = x1;
		x1 = x2;
	}
	wynik.wynik = x1;
	return wynik;
}

void wypisz_iteracje(const char* metoda, Wyniki& wynik) {
	if (wynik.iter_count == 0 && isnan(wynik.wynik)) {
		return;
	}
	cout << "Kolejne przyblizenia miejsc zerowych - " << metoda << ":" << endl;
	for (int i = 0; i < wynik.iter_count; i++) {
		cout << "Iteracja " << i << ": " << fixed << setprecision(8) << wynik.iteracje[i] << endl;
	}
	cout << endl;
}

//funkcja do szukania bledu
double min_blad(const double wynik, const double* x_dokladne, int ile_dokladnych) {
	if (isnan(wynik)) {
		return NAN; //jeœli metoda nie znalaz³a wyniku, b³¹d jest NAN
	}
	if (ile_dokladnych == 0) {
		return NAN; //jeœli nie ma dok³adnych wartoœci do porównania
	}
	double min_blad = numeric_limits<double>::infinity();
	for (int i = 0; i < ile_dokladnych; i++) {
		double blad = fabs(wynik - x_dokladne[i]);
		if (blad < min_blad) {
			min_blad = blad;
		}
	}
	return min_blad;
}

void wynikiFunkcji(double (*f)(double), double (*df)(double), const char* metoda, double start, double koniec, const double* x_dokladne, int ile_dokladnych) {
	cout << endl;
	cout << "Funkcja: " << metoda << endl;

	double krok = 0.1;

	for (double a = start; a < koniec; a += krok) {
		double b = a + krok;
		if (b > koniec) {
			b = koniec;     //upewnienie sie, ze nie wyjdzie poza zakres
		}

		double fa = f(a);
		double fb = f(b);

		if (isnan(fa) || isnan(fb) || isinf(fa) || isinf(fb)) {     //pominiecie przedzialu jestli ktorys kraniec jest nan/inf
			continue;
		}

		//sprawdzenie, czy wystepuje zmiana znaku - jesli nie ma zmiany znaku, nie ma pierwiastka
		if (fa * fb >= 0) {
			continue;
		}

		//punkt startowy dla newtona
		double x0_newton = (a + b) / 2.0;

		Wyniki bis = bisekcja(f, a, b, EPSILON);
		Wyniki newt = newton(f, df, x0_newton, EPSILON);
		Wyniki siecz = sieczna(f, a, b, EPSILON);

		double blad_bis = min_blad(bis.wynik, x_dokladne, ile_dokladnych);
		double blad_newt = min_blad(newt.wynik, x_dokladne, ile_dokladnych);
		double blad_siecz = min_blad(siecz.wynik, x_dokladne, ile_dokladnych);

		cout << fixed << setprecision(13);
		if (blad_bis < 3 * krok) {            //ustawi³am pewn¹ tolerancjê b³êdu, poniewa¿ program znajdywa³ zbyt du¿o zbyt odleg³ych pierwiastków od prawid³owych wyników
			cout << endl;
			cout << "Przedzial: [" << a << ", " << b << "]" << endl;
			cout << "Bisekcja: x = " << bis.wynik << ", blad bezwzgledny = " << blad_bis << endl;
			wypisz_iteracje("Bisekcja", bis);
		}

		cout << fixed << setprecision(13);
		if (blad_newt < 3 * krok) {
			cout << endl;
			cout << "Przedzial: [" << a << ", " << b << "]" << endl;
			cout << "Newton: x = " << newt.wynik << ", blad bezwzgledny = " << blad_newt << endl;
			wypisz_iteracje("Newton", newt);
		}

		cout << fixed << setprecision(13);
		if (blad_siecz < 3 * krok) {
			cout << endl;
			cout << "Przedzial: [" << a << ", " << b << "]" << endl;
			cout << "Sieczne: x = " << siecz.wynik << ", blad bezwzgledny = " << blad_siecz << endl;
			wypisz_iteracje("Sieczne", siecz);
		}
	}
}

Wyniki falszywa_linia(double (*f)(double), double a, double b, double epsilon, const double* x_dokladne, int ile_dokladnych) {
	Wyniki wynik;
	wynik.iter_count = 0;
	wynik.wynik = NAN;

	double fa = f(a);
	double fb = f(b);

	if (isnan(fa) || isnan(fb) || fa * fb >= 0) {
		return wynik; //brak zmiany znaku lub niedefiniowana funkcja
	}

	double x = (a * fb - b * fa) / (fb - fa);  //pierwszy punkt przeciecia prostej z osia X
	double fx = f(x);

	while (wynik.iter_count < MAX_ITER && !isnan(fx) && fabs(fx) >= epsilon) {
		wynik.iteracje[wynik.iter_count] = x;
		wynik.bledy[wynik.iter_count] = min_blad(x, x_dokladne, ile_dokladnych);
		wynik.iter_count++;

		if (fa * fx < 0) {
			b = x;
			fb = fx;
		}
		else {
			a = x;
			fa = fx;
		}

		//punkt przeciecia w kolejnych iteracjach
		x = (a * fb - b * fa) / (fb - fa);
		fx = f(x);
	}

	if (!isnan(fx) && fabs(fx) < epsilon) {
		wynik.iteracje[wynik.iter_count] = x;
		wynik.bledy[wynik.iter_count] = min_blad(x, x_dokladne, ile_dokladnych);
		wynik.iter_count++;
		wynik.wynik = x;
	}

	return wynik;
}