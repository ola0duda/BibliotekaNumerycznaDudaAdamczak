double horner(double wartosc, double* wspolczynniki, int rozmiar); /**
Opis dzia�ania: Oblicza warto�� wielomianu metod� Hornera.
 Opis argument�w:
 wartosc � punkt, w kt�rym obliczamy warto�� wielomianu
 wspolczynniki � tablica wsp�czynnik�w wielomianu
 rozmiar � liczba wsp�czynnik�w
 
 Warto�ci zwracane: Warto�� wielomianu w punkcie wartosc
 Przyk�ad u�ycia:
 double wsp[] = {1, -2, 1}; // x^2 - 2x + 1
 double wynik = horner(3.0, wsp, 3); // wynik = 4
 */

double calka_prostokat(double* przedzial, double* ai, int stopien, int n); /**
Opis dzia�ania: Oblicza ca�k� oznaczon� wielomianu metod� prostok�t�w.
 Opis argument�w:
 przedzial � tablica 2-elementowa: [a, b]
 ai � wsp�czynniki wielomianu
 stopien � stopie� wielomianu
 n � liczba przedzia��w podzia�u
 
 Warto�ci zwracane: Przybli�ona warto�� ca�ki oznaczonej
 Przyk�ad u�ycia:
 double przedzial[] = {0, 1};
 double ai[] = {0, 0, 1}; // x^2
 double wynik = calka_prostokat(przedzial, ai, 2, 1000);
 */

double calka_trapez(double* przedzial, double* ai, int stopien, int n); /**
Opis dzia�ania: Oblicza ca�k� oznaczon� wielomianu metod� trapez�w.
 Opis argument�w:
 przedzial � tablica 2-elementowa: [a, b]
 ai � wsp�czynniki wielomianu
 stopien � stopie� wielomianu
 n � liczba przedzia��w podzia�u
 
 Warto�ci zwracane: Przybli�ona warto�� ca�ki.
 Przyk�ad u�ycia:
 double wynik = calka_trapez(przedzial, ai, 2, 1000);
 */

double calka_simpson(double* przedzial, double* ai, int stopien, int n); /**
Opis dzia�ania: Oblicza ca�k� oznaczon� wielomianu metod� Simpsona (n � liczba par przedzia��w).
 Opis argument�w:
 przedzial � tablica 2-elementowa: [a, b]
 ai � wsp�czynniki wielomianu
 stopien � stopie� wielomianu
 n � liczba przedzia��w podzia�u
 
 Warto�ci zwracane: Przybli�ona warto�� ca�ki
 Przyk�ad u�ycia:
 double wynik = calka_simpson(przedzial, ai, 2, 500); // 2n = 1000 punkt�w
 */

double calka_prostokat_funkcja(double* przedzial, double(*f)(double), int n); /**
Opis dzia�ania: Oblicza ca�k� funkcji metod� prostok�t�w.
 Opis argument�w:
 przedzial � tablica 2-elementowa: [a, b]
 f � wska�nik na funkcj� double f(double)
 n � liczba przedzia��w podzia�u

 Warto�ci zwracane: Przybli�ona warto�� ca�ki
 Przyk�ad u�ycia:
 double wynik = calka_prostokat_funkcja(przedzial, sin, 1000);
 */

double calka_trapez_funkcja(double* przedzial, double(*f)(double), int n); /**
Opis dzia�ania: Ca�kowanie funkcji metod� trapez�w.
 Opis argument�w:
 przedzial � tablica 2-elementowa: [a, b]
 f � wska�nik na funkcj� double f(double)
 n � liczba przedzia��w podzia�u
 
 Warto�ci zwracane: Przybli�ona warto�� ca�ki
 Przyk�ad u�ycia:
 double wynik = calka_trapez_funkcja(przedzial, exp, 1000);
 */

double calka_simpson_funkcja(double* przedzial, double(*f)(double), int n); /**
Opis dzia�ania: Ca�kowanie funkcji metod� Simpsona (n � liczba par przedzia��w).
 Opis argument�w:
 przedzial � tablica 2-elementowa: [a, b]
 f � wska�nik na funkcj� double f(double)
 n � liczba przedzia��w podzia�u
 
 Warto�ci zwracane: Przybli�ona warto�� ca�ki
 Przyk�ad u�ycia:
 double wynik = calka_simpson_funkcja(przedzial, cos, 500);
 */

void GL(int n, double x[], double w[]); /**
Opis dzia�ania: Inicjalizuje tablice x i w w�z�ami i wagami Gaussa-Legendre�a dla n w�z��w (2, 3, lub 4).
 Opis argument�w:
 n � liczba w�z��w
 x � tablica na w�z�y
 w � tablica na wagi
 
 Warto�ci zwracane: Dane przez x[] i w[] (przez referencj�).
 Przyk�ad u�ycia:
 double x[4], w[4];
GL(4, x, w); // x[] i w[] wype�nione
 */

double kwadraturaGL1(double (*f)(double), double a, double b, int n, int m); /**
Opis dzia�ania: Oblicza ca�k� funkcji f metod� Gaussa-Legendre�a (n w�z��w, m podprzedzia��w).
 Opis argument�w:
 f � wska�nik na funkcj�
 a, b � granice ca�kowania
 n � liczba w�z��w
 m � liczba podprzedzia��w (podzia� ca�ego przedzia�u)
 
 Warto�ci zwracane: Przybli�ona warto�� ca�ki.
 Przyk�ad u�ycia:
 double wynik = kwadraturaGL1(exp, 0.0, 1.0, 4, 1000);
 */

double wielomian(double x); /**
Opis dzia�ania: Zwraca warto�� wielomianu w punkcie x, korzystaj�c z globalnej tablicy ai[].
Uwaga: Funkcja pomocnicza, zale�na od globalnej tablicy ai[].
 Opis argument�w:
 x - punkt, w kt�rym chcemy obliczy� wielomian
 
 Warto�ci zwracane: warto�� wielomianu w punkcie x
 Przyk�ad u�ycia:
 ai = new double[5]{0, 1, 0, 0, 1}; // x + x^4
 double wynik = wielomian(2.0); // wynik = 2 + 16 = 18
 */