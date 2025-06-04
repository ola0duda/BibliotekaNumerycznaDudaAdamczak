double lagrange(double punkt, double rozmiar, double* tX, double* tY, int odlegloscWezly); /**
Opis dzia³ania: Oblicza wartoœæ funkcji interpoluj¹cej w punkcie punkt za pomoc¹ interpolacji Lagrange’a.
 Opis argumentów:
 punkt – wartoœæ x, dla której chcemy obliczyæ wartoœæ wielomianu,
 rozmiar – liczba elementów w tablicach tX i tY,
 tX – tablica wartoœci x wêz³ów interpolacji,
 tY – tablica wartoœci y funkcji w punktach tX,
 odlegloscWezly – krok indeksowania (zwykle 1; jeœli dane s¹ rzadziej rozmieszczone – wiêkszy).
 
 Wartoœci zwracane: Wartoœæ wielomianu interpolacyjnego w punkcie punkt.
 Przyk³ad u¿ycia:
 double x[] = {1, 2, 3};
 double y[] = {2, 4, 6};
 double wynik = lagrange(2.5, 3, x, y, 1);
 */

double naturalna(double wartosc, double* wspolczynniki, int rozmiar); /**
Opis dzia³ania: Oblicza wartoœæ wielomianu w postaci naturalnej (standardowej) w punkcie wartosc.
 Opis argumentów:
 wartosc – punkt x,
 wspolczynniki – wspó³czynniki wielomianu: a0 + a1*x + a2*x^2 + ...,
 rozmiar – liczba wspó³czynników.
 
 Wartoœci zwracane: Wartoœæ wielomianu w punkcie wartosc.
 Przyk³ad u¿ycia:
 double wsp[] = {1, 2, 3}; // 1 + 2x + 3x^2
 double wynik = naturalna(2.0, wsp, 3); // = 1 + 4 + 12 = 17
 */

double horner(double wartosc, double* wspolczynniki, int rozmiar); /**
Opis dzia³ania: Oblicza wartoœæ wielomianu metod¹ Hornera w punkcie wartosc.
 Opis argumentów:
 wartosc – punkt x,
 wspolczynniki – wspó³czynniki wielomianu: a0 + a1*x + a2*x^2 + ...,
 rozmiar – liczba wspó³czynników.
 
 Wartoœci zwracane: Wartoœæ wielomianu w punkcie wartosc.
 Przyk³ad u¿ycia:
 double wsp[] = {1, 2, 3}; // 1 + 2x + 3x^2
 double wynik = horner(2.0, wsp, 3); // = 3*2 + 2 = 8, 8*2 + 1 = 17
 */

void ilorazy(double** f, double* x, double* y, int n); /**
Opis dzia³ania: Wype³nia tablicê f ilorazami ró¿nicowymi potrzebnymi do interpolacji Newtona. f[i][j] to j-ty rz¹d ilorazu ró¿nicowego dla i-tego wêz³a.
 Opis argumentów:
 f – dwuwymiarowa tablica (macierz) ró¿nicowych ilorazów, modyfikowana w miejscu,
 x – tablica wêz³ów interpolacji x,
 y – tablica wartoœci funkcji w punktach x,
 n – liczba punktów.

Wartoœci zwracane: Nic (dzia³a przez modyfikacjê tablicy f)
 Przyk³ad u¿ycia:
double x[] = {1, 2, 3};
double y[] = {1, 4, 9};
double** f = new double*[3];
for (int i = 0; i < 3; ++i) f[i] = new double[3];
ilorazy(f, x, y, 3);
 */

double newton(double wartosc, double* a, double* x, int rozmiar); /**
Opis dzia³ania: Oblicza wartoœæ funkcji interpoluj¹cej w punkcie wartosc dla wielomianu Newtona z gotowymi wspó³czynnikami a.
 Opis argumentów:
 wartosc – punkt x, dla którego szukamy wartoœci,
 a – wspó³czynniki (ilorazy ró¿nicowe z funkcji ilorazy()),
 x – tablica wêz³ów interpolacji,
 rozmiar – liczba wêz³ów.
 
 Wartoœci zwracane: Wartoœæ wielomianu interpolacyjnego Newtona w wartosc.
 Przyk³ad u¿ycia:
double x[] = {1, 2, 3};
double y[] = {1, 4, 9};
double** f = new double*[3];
for (int i = 0; i < 3; ++i) f[i] = new double[3];
double wynik = newton(2.5, f[2], x, 3);
 */