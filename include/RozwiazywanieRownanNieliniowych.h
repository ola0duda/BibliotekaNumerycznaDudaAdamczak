const int MAX_ITER = 200; //maksymalna ilosc iteracji w funkcjach (warunek koncowy)

struct Wyniki {		//struktura do prostszego zapisywania i rozdzielania wyników na wartosc miejsca zerowego, iteracje, bledy
	double wynik;	//przybli¿on¹ wartoœæ miejsca zerowego
	double iteracje[MAX_ITER];		//tablicê kolejnych przybli¿eñ
	double bledy[MAX_ITER];		//tablicê b³êdów przybli¿eñ wzglêdem dok³adnych rozwi¹zañ,
	int iter_count;		//liczbê wykonanych iteracji
};

Wyniki bisekcja(double (*f)(double), double a, double b, double epsilon); /**
Opis dzia³ania: Metoda bisekcji do znajdowania miejsca zerowego funkcji. Dzia³a na przedziale [a, b], na którym funkcja przyjmuje wartoœci o przeciwnych znakach. Metoda dzieli przedzia³ na pó³ i wybiera podprzedzia³, w którym wystêpuje zmiana znaku, iteruj¹c a¿ do zbie¿noœci.
 Opis argumentów:
 f - WskaŸnik na funkcjê pojedynczej zmiennej typu double,
 a - Lewy kraniec przedzia³u.
 b - Prawy kraniec przedzia³u.
 epsilon - Dok³adnoœæ zbie¿noœci (tolerancja b³êdu).
 
 Wartoœci zwracane: struktura wyniki zawieraj¹ca przybli¿on¹ wartoœæ miejsca zerowego, tablicê kolejnych przybli¿eñ, liczbê wykonanych iteracji
 Przyk³ad u¿ycia:
 auto f = [](double x) { return x*x - 2; };
  Wyniki wynik = bisekcja(f, 0, 2, 1e-6);
  std::cout << "Pierwiastek: " << wynik.wynik << std::endl;
 */

Wyniki newton(double (*f)(double), double (*df)(double), double x0_start, double epsilon); /**
Opis dzia³ania: Metoda Newtona do znajdowania miejsca zerowego funkcji. Metoda wykorzystuje wzór iteracyjny x_{n+1} = x_n - f(x_n)/f'(x_n). Wymaga podania funkcji oraz jej pochodnej.
 Opis argumentów:
 f - WskaŸnik na funkcjê pojedynczej zmiennej typu double.
 df - WskaŸnik na pochodn¹ funkcji f.
 epsilon - Dok³adnoœæ zbie¿noœci (tolerancja b³êdu).
 
 Wartoœci zwracane: struktura wyniki zawieraj¹ca przybli¿on¹ wartoœæ miejsca zerowego, tablicê kolejnych przybli¿eñ, liczbê wykonanych iteracji
 Przyk³ad u¿ycia:
 auto f = [](double x) { return x*x - 2; };
  auto df = [](double x) { return 2*x; };
  Wyniki wynik = newton(f, df, 1.0, 1e-6);
  std::cout << "Pierwiastek: " << wynik.wynik << std::endl;
 */

Wyniki sieczna(double (*f)(double), double x0_start, double x1_start, double epsilon); /**
Opis dzia³ania: Metoda siecznych do znajdowania miejsca zerowego funkcji. Wykorzystuje przybli¿enia x0 i x1 oraz iteracyjnie wyznacza kolejne punkty na podstawie stycznych przechodz¹cych przez punkty (x0,f(x0)) i (x1,f(x1)). Nie wymaga znajomoœci pochodnej funkcji.
 Opis argumentów:
 f - WskaŸnik na funkcjê pojedynczej zmiennej typu double.
 x0_start - Pierwszy punkt startowy iteracji.
 x1_start - Drugi punkt startowy iteracji.
 epsilon - Dok³adnoœæ zbie¿noœci (tolerancja b³êdu).

 Wartoœci zwracane: struktura wyniki zawieraj¹ca przybli¿on¹ wartoœæ miejsca zerowego, tablicê kolejnych przybli¿eñ, liczbê wykonanych iteracji
 Przyk³ad u¿ycia:
 auto f = [](double x) { return x*x - 2; };
  Wyniki wynik = sieczna(f, 0.0, 2.0, 1e-6);
  std::cout << "Pierwiastek: " << wynik.wynik << std::endl;
 */

void wypisz_iteracje(const char* metoda, Wyniki& wynik); /**
Opis dzia³ania: Funkcja wypisuj¹ca kolejne przybli¿enia iteracji danej metody na standardowe wyjœcie.
 Opis argumentów:
 metoda - Nazwa metody (np. "Bisekcja", "Newton").
 wynik - Struktura Wyniki zawieraj¹ca iteracje i liczbê wykonanych iteracji.
 
 Wartoœci zwracane: brak (wypisuje kolejne przybli¿enia iteracji danej metody)
 Przyk³ad u¿ycia:
 Wyniki wynik = bisekcja(f, 0, 2, 1e-6);
  wypisz_iteracje("Bisekcja", wynik);
 */

double min_blad(const double wynik, const double* x_dokladne, int ile_dokladnych); /**
Opis dzia³ania: Funkcja oblicza minimalny b³¹d bezwzglêdny miêdzy wynikiem metody a zestawem znanych dok³adnych miejsc zerowych.
 Opis argumentów:
 wynik - Wynik obliczony przez metodê.
 x_dokladne - WskaŸnik na tablicê znanych dok³adnych rozwi¹zañ.
 ile_dokladnych - Liczba dok³adnych rozwi¹zañ w tablicy.
 
 Wartoœci zwracane: Minimalny b³¹d bezwzglêdny. Jeœli wynik jest NAN lub brak dok³adnych wartoœci, zwraca NAN.
 Przyk³ad u¿ycia:
 double dokladne[] = {1.414213562};
  double blad = min_blad(wynik.wynik, dokladne, 1);
  std::cout << "Minimalny b³¹d: " << blad << std::endl;
 */

void wynikiFunkcji(double (*f)(double), double (*df)(double), const char* metoda, double start, double koniec, const double* x_dokladne, int ile_dokladnych); /**
Opis dzia³ania: Funkcja przeprowadza testy na podanym przedziale, sprawdzaj¹c metody bisekcji, Newtona i siecznych. Dla ka¿dego przedzia³u z wykryt¹ zmian¹ znaku funkcji wypisuje znalezione przybli¿enia i ich b³êdy.
 Opis argumentów:
 f - WskaŸnik na funkcjê.
 df - WskaŸnik na pochodn¹ funkcji (dla metody Newtona).
 metoda - Nazwa testowanej funkcji (do wypisania).
 start - Pocz¹tek przedzia³u.
 koniec - Koniec przedzia³u.
 x_dokladne - WskaŸnik na tablicê dok³adnych rozwi¹zañ.
 ile_dokladnych - Liczba dok³adnych rozwi¹zañ.
 
 Wartoœci zwracane: brak (wypisuje przyblizenia i bledy)
 Przyk³ad u¿ycia:
 auto f = [](double x) { return x*x - 2; };
  auto df = [](double x) { return 2*x; };
  double dokladne[] = {1.414213562};
  wynikiFunkcji(f, df, "x^2-2", 0, 2, dokladne, 1);
 */

Wyniki falszywa_linia(double (*f)(double), double a, double b, double epsilon, const double* x_dokladne, int ile_dokladnych); /**
Opis dzia³ania: Metoda fa³szywej linii (regula falsi) do znajdowania miejsca zerowego funkcji. Podobna do bisekcji, ale zamiast œrodka przedzia³u u¿ywa punktu przeciêcia prostej ³¹cz¹cej koñce przedzia³u z osi¹ X.
 Opis argumentów:
 f - WskaŸnik na funkcjê pojedynczej zmiennej typu double.
 a - Lewy kraniec przedzia³u.
 b - Prawy kraniec przedzia³u.
 epsilon - Dok³adnoœæ zbie¿noœci (tolerancja b³êdu).
 x_dokladne - WskaŸnik na tablicê dok³adnych rozwi¹zañ.
 ile_dokladnych - Liczba dok³adnych rozwi¹zañ.
 
 Wartoœci zwracane: struktura wyniki zawieraj¹ca przybli¿on¹ wartoœæ miejsca zerowego, tablicê kolejnych przybli¿eñ, tablicê b³êdów przybli¿eñ wzglêdem dok³adnych rozwi¹zañ, liczbê wykonanych iteracji
 Przyk³ad u¿ycia:
 auto f = [](double x) { return x*x - 2; };
  double dokladne[] = {1.414213562};
  Wyniki wynik = falszywa_linia(f, 0, 2, 1e-6, dokladne, 1);
  std::cout << "Pierwiastek: " << wynik.wynik << std::endl;
 */