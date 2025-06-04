#include <functional>


void trojkat(double** A, double* b, double* x, int n); /**
Opis dzia�ania: Rozwi�zuje uk�ad r�wna� liniowych Ax = b metod� podstawiania wstecznego, zak�adaj�c, �e macierz A jest tr�jk�tna g�rna.
 Opis argument�w:
 A - wska�nik do macierzy wsp�czynnik�w (n x n) tr�jk�tnej g�rnej
 b - wektor wyraz�w wolnych (rozmiar n)
 x - wektor wynikowy (rozmiar n), wype�niany rozwi�zaniami
 n - liczba r�wna� i niewiadomych
 
 Warto�ci zwracane: brak (x jest wype�niane wynikami)
 Przyk�ad u�ycia:
 trojkat(A, b, x, n);
 */

double baza(double x, int n); /**
Opis dzia�ania: Oblicza warto�� funkcji bazowej dla danego x i stopnia n. W tym przypadku jest to po prostu x^n (baza wielomianu pot�gowego).
 Opis argument�w:
 x - argument funkcji bazowej
 n - stopie� funkcji bazowej
 
 Warto�ci zwracane: warto�� x podniesionego do pot�gi n
 Przyk�ad u�ycia:
 double y = baza(2.0, 3); // zwr�ci 8.0
 */


double kwadraturaGL2(std::function<double(double)> f, double a, double b); /**
Opis dzia�ania: Oblicza ca�k� oznaczon� funkcji f na przedziale [a, b] przy u�yciu kwadratury Gaussa-Legendre'a z 4 w�z�ami, z podzia�em przedzia�u ca�kowania na 10000 podprzedzia��w.
 Opis argument�w:
 f - funkcja, kt�r� ca�kujemy
 a - pocz�tek przedzia�u ca�kowania
 b - koniec przedzia�u ca�kowania
 
 Warto�ci zwracane: przybli�ona warto�� ca�ki oznaczonej
 Przyk�ad u�ycia:
 double I = kwadraturaGL2([](double x){ return sin(x); }, 0.0, M_PI);
 */

void gauss2(double** A, double* b, double* x, int n); /**
Opis dzia�ania: Rozwi�zuje uk�ad r�wna� liniowych Ax = b metod� eliminacji Gaussa z cz�ciowym wyborem elementu g��wnego, a nast�pnie stosuje podstawienie wsteczne.
 Opis argument�w:
 A - wska�nik do macierzy wsp�czynnik�w (n x n)
 b - wektor wyraz�w wolnych (rozmiar n)
 x - wektor wynikowy (rozmiar n)
 n - liczba r�wna� i niewiadomych
 
 Warto�ci zwracane: brak (x zostaje wype�niony wynikami)
 Przyk�ad u�ycia:
 gauss2(A, b, x, n);
 */

void aproksymacja(double a, double b, int stopien, double (*funkcja)(double)); /**
Opis dzia�ania: Przeprowadza aproksymacj� funkcji `funkcja` na przedziale [a, b] za pomoc� wielomianu w bazie pot�gowej o zadanym stopniu. Oblicza wsp�czynniki wielomianu metod� rzutowania na przestrze� bazow�, korzystaj�c z metody najmniejszych kwadrat�w oraz kwadratury Gaussa-Legendre'a. Na ko�cu wypisuje wsp�czynniki oraz por�wnanie warto�ci aproksymacji i funkcji dla wybranych punkt�w z b��dem bezwzgl�dnym.
 Opis argument�w:
 a - pocz�tek przedzia�u aproksymacji
 b - koniec przedzia�u aproksymacji
 stopien - stopie� wielomianu aproksymacyjnego
 funkcja - wska�nik do funkcji aproksymowanej: double f(double)
 
 Warto�ci zwracane: brak (wyniki wypisywane s� na ekran)
 Przyk�ad u�ycia:
 aproksymacja(0.0, 1.0, 4, [](double x) { return exp(x); });
 */