#include <functional>


void trojkat(double** A, double* b, double* x, int n); /**
Opis dzia³ania: Rozwi¹zuje uk³ad równañ liniowych Ax = b metod¹ podstawiania wstecznego, zak³adaj¹c, ¿e macierz A jest trójk¹tna górna.
 Opis argumentów:
 A - wskaŸnik do macierzy wspó³czynników (n x n) trójk¹tnej górnej
 b - wektor wyrazów wolnych (rozmiar n)
 x - wektor wynikowy (rozmiar n), wype³niany rozwi¹zaniami
 n - liczba równañ i niewiadomych
 
 Wartoœci zwracane: brak (x jest wype³niane wynikami)
 Przyk³ad u¿ycia:
 trojkat(A, b, x, n);
 */

double baza(double x, int n); /**
Opis dzia³ania: Oblicza wartoœæ funkcji bazowej dla danego x i stopnia n. W tym przypadku jest to po prostu x^n (baza wielomianu potêgowego).
 Opis argumentów:
 x - argument funkcji bazowej
 n - stopieñ funkcji bazowej
 
 Wartoœci zwracane: wartoœæ x podniesionego do potêgi n
 Przyk³ad u¿ycia:
 double y = baza(2.0, 3); // zwróci 8.0
 */


double kwadraturaGL2(std::function<double(double)> f, double a, double b); /**
Opis dzia³ania: Oblicza ca³kê oznaczon¹ funkcji f na przedziale [a, b] przy u¿yciu kwadratury Gaussa-Legendre'a z 4 wêz³ami, z podzia³em przedzia³u ca³kowania na 10000 podprzedzia³ów.
 Opis argumentów:
 f - funkcja, któr¹ ca³kujemy
 a - pocz¹tek przedzia³u ca³kowania
 b - koniec przedzia³u ca³kowania
 
 Wartoœci zwracane: przybli¿ona wartoœæ ca³ki oznaczonej
 Przyk³ad u¿ycia:
 double I = kwadraturaGL2([](double x){ return sin(x); }, 0.0, M_PI);
 */

void gauss2(double** A, double* b, double* x, int n); /**
Opis dzia³ania: Rozwi¹zuje uk³ad równañ liniowych Ax = b metod¹ eliminacji Gaussa z czêœciowym wyborem elementu g³ównego, a nastêpnie stosuje podstawienie wsteczne.
 Opis argumentów:
 A - wskaŸnik do macierzy wspó³czynników (n x n)
 b - wektor wyrazów wolnych (rozmiar n)
 x - wektor wynikowy (rozmiar n)
 n - liczba równañ i niewiadomych
 
 Wartoœci zwracane: brak (x zostaje wype³niony wynikami)
 Przyk³ad u¿ycia:
 gauss2(A, b, x, n);
 */

void aproksymacja(double a, double b, int stopien, double (*funkcja)(double)); /**
Opis dzia³ania: Przeprowadza aproksymacjê funkcji `funkcja` na przedziale [a, b] za pomoc¹ wielomianu w bazie potêgowej o zadanym stopniu. Oblicza wspó³czynniki wielomianu metod¹ rzutowania na przestrzeñ bazow¹, korzystaj¹c z metody najmniejszych kwadratów oraz kwadratury Gaussa-Legendre'a. Na koñcu wypisuje wspó³czynniki oraz porównanie wartoœci aproksymacji i funkcji dla wybranych punktów z b³êdem bezwzglêdnym.
 Opis argumentów:
 a - pocz¹tek przedzia³u aproksymacji
 b - koniec przedzia³u aproksymacji
 stopien - stopieñ wielomianu aproksymacyjnego
 funkcja - wskaŸnik do funkcji aproksymowanej: double f(double)
 
 Wartoœci zwracane: brak (wyniki wypisywane s¹ na ekran)
 Przyk³ad u¿ycia:
 aproksymacja(0.0, 1.0, 4, [](double x) { return exp(x); });
 */