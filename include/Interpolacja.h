double lagrange(double punkt, double rozmiar, double* tX, double* tY, int odlegloscWezly); /**
Opis dzia�ania: Oblicza warto�� funkcji interpoluj�cej w punkcie punkt za pomoc� interpolacji Lagrange�a.
 Opis argument�w:
 punkt � warto�� x, dla kt�rej chcemy obliczy� warto�� wielomianu,
 rozmiar � liczba element�w w tablicach tX i tY,
 tX � tablica warto�ci x w�z��w interpolacji,
 tY � tablica warto�ci y funkcji w punktach tX,
 odlegloscWezly � krok indeksowania (zwykle 1; je�li dane s� rzadziej rozmieszczone � wi�kszy).
 
 Warto�ci zwracane: Warto�� wielomianu interpolacyjnego w punkcie punkt.
 Przyk�ad u�ycia:
 double x[] = {1, 2, 3};
 double y[] = {2, 4, 6};
 double wynik = lagrange(2.5, 3, x, y, 1);
 */

double naturalna(double wartosc, double* wspolczynniki, int rozmiar); /**
Opis dzia�ania: Oblicza warto�� wielomianu w postaci naturalnej (standardowej) w punkcie wartosc.
 Opis argument�w:
 wartosc � punkt x,
 wspolczynniki � wsp�czynniki wielomianu: a0 + a1*x + a2*x^2 + ...,
 rozmiar � liczba wsp�czynnik�w.
 
 Warto�ci zwracane: Warto�� wielomianu w punkcie wartosc.
 Przyk�ad u�ycia:
 double wsp[] = {1, 2, 3}; // 1 + 2x + 3x^2
 double wynik = naturalna(2.0, wsp, 3); // = 1 + 4 + 12 = 17
 */

double horner(double wartosc, double* wspolczynniki, int rozmiar); /**
Opis dzia�ania: Oblicza warto�� wielomianu metod� Hornera w punkcie wartosc.
 Opis argument�w:
 wartosc � punkt x,
 wspolczynniki � wsp�czynniki wielomianu: a0 + a1*x + a2*x^2 + ...,
 rozmiar � liczba wsp�czynnik�w.
 
 Warto�ci zwracane: Warto�� wielomianu w punkcie wartosc.
 Przyk�ad u�ycia:
 double wsp[] = {1, 2, 3}; // 1 + 2x + 3x^2
 double wynik = horner(2.0, wsp, 3); // = 3*2 + 2 = 8, 8*2 + 1 = 17
 */

void ilorazy(double** f, double* x, double* y, int n); /**
Opis dzia�ania: Wype�nia tablic� f ilorazami r�nicowymi potrzebnymi do interpolacji Newtona. f[i][j] to j-ty rz�d ilorazu r�nicowego dla i-tego w�z�a.
 Opis argument�w:
 f � dwuwymiarowa tablica (macierz) r�nicowych iloraz�w, modyfikowana w miejscu,
 x � tablica w�z��w interpolacji x,
 y � tablica warto�ci funkcji w punktach x,
 n � liczba punkt�w.

Warto�ci zwracane: Nic (dzia�a przez modyfikacj� tablicy f)
 Przyk�ad u�ycia:
double x[] = {1, 2, 3};
double y[] = {1, 4, 9};
double** f = new double*[3];
for (int i = 0; i < 3; ++i) f[i] = new double[3];
ilorazy(f, x, y, 3);
 */

double newton(double wartosc, double* a, double* x, int rozmiar); /**
Opis dzia�ania: Oblicza warto�� funkcji interpoluj�cej w punkcie wartosc dla wielomianu Newtona z gotowymi wsp�czynnikami a.
 Opis argument�w:
 wartosc � punkt x, dla kt�rego szukamy warto�ci,
 a � wsp�czynniki (ilorazy r�nicowe z funkcji ilorazy()),
 x � tablica w�z��w interpolacji,
 rozmiar � liczba w�z��w.
 
 Warto�ci zwracane: Warto�� wielomianu interpolacyjnego Newtona w wartosc.
 Przyk�ad u�ycia:
double x[] = {1, 2, 3};
double y[] = {1, 4, 9};
double** f = new double*[3];
for (int i = 0; i < 3; ++i) f[i] = new double[3];
double wynik = newton(2.5, f[2], x, 3);
 */