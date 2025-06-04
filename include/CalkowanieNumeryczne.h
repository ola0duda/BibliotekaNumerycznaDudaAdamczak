double horner(double wartosc, double* wspolczynniki, int rozmiar); /**
Opis dzia³ania: Oblicza wartoœæ wielomianu metod¹ Hornera.
 Opis argumentów:
 wartosc – punkt, w którym obliczamy wartoœæ wielomianu
 wspolczynniki – tablica wspó³czynników wielomianu
 rozmiar – liczba wspó³czynników
 
 Wartoœci zwracane: Wartoœæ wielomianu w punkcie wartosc
 Przyk³ad u¿ycia:
 double wsp[] = {1, -2, 1}; // x^2 - 2x + 1
 double wynik = horner(3.0, wsp, 3); // wynik = 4
 */

double calka_prostokat(double* przedzial, double* ai, int stopien, int n); /**
Opis dzia³ania: Oblicza ca³kê oznaczon¹ wielomianu metod¹ prostok¹tów.
 Opis argumentów:
 przedzial – tablica 2-elementowa: [a, b]
 ai – wspó³czynniki wielomianu
 stopien – stopieñ wielomianu
 n – liczba przedzia³ów podzia³u
 
 Wartoœci zwracane: Przybli¿ona wartoœæ ca³ki oznaczonej
 Przyk³ad u¿ycia:
 double przedzial[] = {0, 1};
 double ai[] = {0, 0, 1}; // x^2
 double wynik = calka_prostokat(przedzial, ai, 2, 1000);
 */

double calka_trapez(double* przedzial, double* ai, int stopien, int n); /**
Opis dzia³ania: Oblicza ca³kê oznaczon¹ wielomianu metod¹ trapezów.
 Opis argumentów:
 przedzial – tablica 2-elementowa: [a, b]
 ai – wspó³czynniki wielomianu
 stopien – stopieñ wielomianu
 n – liczba przedzia³ów podzia³u
 
 Wartoœci zwracane: Przybli¿ona wartoœæ ca³ki.
 Przyk³ad u¿ycia:
 double wynik = calka_trapez(przedzial, ai, 2, 1000);
 */

double calka_simpson(double* przedzial, double* ai, int stopien, int n); /**
Opis dzia³ania: Oblicza ca³kê oznaczon¹ wielomianu metod¹ Simpsona (n – liczba par przedzia³ów).
 Opis argumentów:
 przedzial – tablica 2-elementowa: [a, b]
 ai – wspó³czynniki wielomianu
 stopien – stopieñ wielomianu
 n – liczba przedzia³ów podzia³u
 
 Wartoœci zwracane: Przybli¿ona wartoœæ ca³ki
 Przyk³ad u¿ycia:
 double wynik = calka_simpson(przedzial, ai, 2, 500); // 2n = 1000 punktów
 */

double calka_prostokat_funkcja(double* przedzial, double(*f)(double), int n); /**
Opis dzia³ania: Oblicza ca³kê funkcji metod¹ prostok¹tów.
 Opis argumentów:
 przedzial – tablica 2-elementowa: [a, b]
 f – wskaŸnik na funkcjê double f(double)
 n – liczba przedzia³ów podzia³u

 Wartoœci zwracane: Przybli¿ona wartoœæ ca³ki
 Przyk³ad u¿ycia:
 double wynik = calka_prostokat_funkcja(przedzial, sin, 1000);
 */

double calka_trapez_funkcja(double* przedzial, double(*f)(double), int n); /**
Opis dzia³ania: Ca³kowanie funkcji metod¹ trapezów.
 Opis argumentów:
 przedzial – tablica 2-elementowa: [a, b]
 f – wskaŸnik na funkcjê double f(double)
 n – liczba przedzia³ów podzia³u
 
 Wartoœci zwracane: Przybli¿ona wartoœæ ca³ki
 Przyk³ad u¿ycia:
 double wynik = calka_trapez_funkcja(przedzial, exp, 1000);
 */

double calka_simpson_funkcja(double* przedzial, double(*f)(double), int n); /**
Opis dzia³ania: Ca³kowanie funkcji metod¹ Simpsona (n – liczba par przedzia³ów).
 Opis argumentów:
 przedzial – tablica 2-elementowa: [a, b]
 f – wskaŸnik na funkcjê double f(double)
 n – liczba przedzia³ów podzia³u
 
 Wartoœci zwracane: Przybli¿ona wartoœæ ca³ki
 Przyk³ad u¿ycia:
 double wynik = calka_simpson_funkcja(przedzial, cos, 500);
 */

void GL(int n, double x[], double w[]); /**
Opis dzia³ania: Inicjalizuje tablice x i w wêz³ami i wagami Gaussa-Legendre’a dla n wêz³ów (2, 3, lub 4).
 Opis argumentów:
 n – liczba wêz³ów
 x – tablica na wêz³y
 w – tablica na wagi
 
 Wartoœci zwracane: Dane przez x[] i w[] (przez referencjê).
 Przyk³ad u¿ycia:
 double x[4], w[4];
GL(4, x, w); // x[] i w[] wype³nione
 */

double kwadraturaGL1(double (*f)(double), double a, double b, int n, int m); /**
Opis dzia³ania: Oblicza ca³kê funkcji f metod¹ Gaussa-Legendre’a (n wêz³ów, m podprzedzia³ów).
 Opis argumentów:
 f – wskaŸnik na funkcjê
 a, b – granice ca³kowania
 n – liczba wêz³ów
 m – liczba podprzedzia³ów (podzia³ ca³ego przedzia³u)
 
 Wartoœci zwracane: Przybli¿ona wartoœæ ca³ki.
 Przyk³ad u¿ycia:
 double wynik = kwadraturaGL1(exp, 0.0, 1.0, 4, 1000);
 */

double wielomian(double x); /**
Opis dzia³ania: Zwraca wartoœæ wielomianu w punkcie x, korzystaj¹c z globalnej tablicy ai[].
Uwaga: Funkcja pomocnicza, zale¿na od globalnej tablicy ai[].
 Opis argumentów:
 x - punkt, w którym chcemy obliczyæ wielomian
 
 Wartoœci zwracane: wartoœæ wielomianu w punkcie x
 Przyk³ad u¿ycia:
 ai = new double[5]{0, 1, 0, 0, 1}; // x + x^4
 double wynik = wielomian(2.0); // wynik = 2 + 16 = 18
 */