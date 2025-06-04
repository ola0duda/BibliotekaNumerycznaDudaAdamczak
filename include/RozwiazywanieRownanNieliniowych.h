const int MAX_ITER = 200; //maksymalna ilosc iteracji w funkcjach (warunek koncowy)

struct Wyniki {		//struktura do prostszego zapisywania i rozdzielania wynik�w na wartosc miejsca zerowego, iteracje, bledy
	double wynik;	//przybli�on� warto�� miejsca zerowego
	double iteracje[MAX_ITER];		//tablic� kolejnych przybli�e�
	double bledy[MAX_ITER];		//tablic� b��d�w przybli�e� wzgl�dem dok�adnych rozwi�za�,
	int iter_count;		//liczb� wykonanych iteracji
};

Wyniki bisekcja(double (*f)(double), double a, double b, double epsilon); /**
Opis dzia�ania: Metoda bisekcji do znajdowania miejsca zerowego funkcji. Dzia�a na przedziale [a, b], na kt�rym funkcja przyjmuje warto�ci o przeciwnych znakach. Metoda dzieli przedzia� na p� i wybiera podprzedzia�, w kt�rym wyst�puje zmiana znaku, iteruj�c a� do zbie�no�ci.
 Opis argument�w:
 f - Wska�nik na funkcj� pojedynczej zmiennej typu double,
 a - Lewy kraniec przedzia�u.
 b - Prawy kraniec przedzia�u.
 epsilon - Dok�adno�� zbie�no�ci (tolerancja b��du).
 
 Warto�ci zwracane: struktura wyniki zawieraj�ca przybli�on� warto�� miejsca zerowego, tablic� kolejnych przybli�e�, liczb� wykonanych iteracji
 Przyk�ad u�ycia:
 auto f = [](double x) { return x*x - 2; };
  Wyniki wynik = bisekcja(f, 0, 2, 1e-6);
  std::cout << "Pierwiastek: " << wynik.wynik << std::endl;
 */

Wyniki newton(double (*f)(double), double (*df)(double), double x0_start, double epsilon); /**
Opis dzia�ania: Metoda Newtona do znajdowania miejsca zerowego funkcji. Metoda wykorzystuje wz�r iteracyjny x_{n+1} = x_n - f(x_n)/f'(x_n). Wymaga podania funkcji oraz jej pochodnej.
 Opis argument�w:
 f - Wska�nik na funkcj� pojedynczej zmiennej typu double.
 df - Wska�nik na pochodn� funkcji f.
 epsilon - Dok�adno�� zbie�no�ci (tolerancja b��du).
 
 Warto�ci zwracane: struktura wyniki zawieraj�ca przybli�on� warto�� miejsca zerowego, tablic� kolejnych przybli�e�, liczb� wykonanych iteracji
 Przyk�ad u�ycia:
 auto f = [](double x) { return x*x - 2; };
  auto df = [](double x) { return 2*x; };
  Wyniki wynik = newton(f, df, 1.0, 1e-6);
  std::cout << "Pierwiastek: " << wynik.wynik << std::endl;
 */

Wyniki sieczna(double (*f)(double), double x0_start, double x1_start, double epsilon); /**
Opis dzia�ania: Metoda siecznych do znajdowania miejsca zerowego funkcji. Wykorzystuje przybli�enia x0 i x1 oraz iteracyjnie wyznacza kolejne punkty na podstawie stycznych przechodz�cych przez punkty (x0,f(x0)) i (x1,f(x1)). Nie wymaga znajomo�ci pochodnej funkcji.
 Opis argument�w:
 f - Wska�nik na funkcj� pojedynczej zmiennej typu double.
 x0_start - Pierwszy punkt startowy iteracji.
 x1_start - Drugi punkt startowy iteracji.
 epsilon - Dok�adno�� zbie�no�ci (tolerancja b��du).

 Warto�ci zwracane: struktura wyniki zawieraj�ca przybli�on� warto�� miejsca zerowego, tablic� kolejnych przybli�e�, liczb� wykonanych iteracji
 Przyk�ad u�ycia:
 auto f = [](double x) { return x*x - 2; };
  Wyniki wynik = sieczna(f, 0.0, 2.0, 1e-6);
  std::cout << "Pierwiastek: " << wynik.wynik << std::endl;
 */

void wypisz_iteracje(const char* metoda, Wyniki& wynik); /**
Opis dzia�ania: Funkcja wypisuj�ca kolejne przybli�enia iteracji danej metody na standardowe wyj�cie.
 Opis argument�w:
 metoda - Nazwa metody (np. "Bisekcja", "Newton").
 wynik - Struktura Wyniki zawieraj�ca iteracje i liczb� wykonanych iteracji.
 
 Warto�ci zwracane: brak (wypisuje kolejne przybli�enia iteracji danej metody)
 Przyk�ad u�ycia:
 Wyniki wynik = bisekcja(f, 0, 2, 1e-6);
  wypisz_iteracje("Bisekcja", wynik);
 */

double min_blad(const double wynik, const double* x_dokladne, int ile_dokladnych); /**
Opis dzia�ania: Funkcja oblicza minimalny b��d bezwzgl�dny mi�dzy wynikiem metody a zestawem znanych dok�adnych miejsc zerowych.
 Opis argument�w:
 wynik - Wynik obliczony przez metod�.
 x_dokladne - Wska�nik na tablic� znanych dok�adnych rozwi�za�.
 ile_dokladnych - Liczba dok�adnych rozwi�za� w tablicy.
 
 Warto�ci zwracane: Minimalny b��d bezwzgl�dny. Je�li wynik jest NAN lub brak dok�adnych warto�ci, zwraca NAN.
 Przyk�ad u�ycia:
 double dokladne[] = {1.414213562};
  double blad = min_blad(wynik.wynik, dokladne, 1);
  std::cout << "Minimalny b��d: " << blad << std::endl;
 */

void wynikiFunkcji(double (*f)(double), double (*df)(double), const char* metoda, double start, double koniec, const double* x_dokladne, int ile_dokladnych); /**
Opis dzia�ania: Funkcja przeprowadza testy na podanym przedziale, sprawdzaj�c metody bisekcji, Newtona i siecznych. Dla ka�dego przedzia�u z wykryt� zmian� znaku funkcji wypisuje znalezione przybli�enia i ich b��dy.
 Opis argument�w:
 f - Wska�nik na funkcj�.
 df - Wska�nik na pochodn� funkcji (dla metody Newtona).
 metoda - Nazwa testowanej funkcji (do wypisania).
 start - Pocz�tek przedzia�u.
 koniec - Koniec przedzia�u.
 x_dokladne - Wska�nik na tablic� dok�adnych rozwi�za�.
 ile_dokladnych - Liczba dok�adnych rozwi�za�.
 
 Warto�ci zwracane: brak (wypisuje przyblizenia i bledy)
 Przyk�ad u�ycia:
 auto f = [](double x) { return x*x - 2; };
  auto df = [](double x) { return 2*x; };
  double dokladne[] = {1.414213562};
  wynikiFunkcji(f, df, "x^2-2", 0, 2, dokladne, 1);
 */

Wyniki falszywa_linia(double (*f)(double), double a, double b, double epsilon, const double* x_dokladne, int ile_dokladnych); /**
Opis dzia�ania: Metoda fa�szywej linii (regula falsi) do znajdowania miejsca zerowego funkcji. Podobna do bisekcji, ale zamiast �rodka przedzia�u u�ywa punktu przeci�cia prostej ��cz�cej ko�ce przedzia�u z osi� X.
 Opis argument�w:
 f - Wska�nik na funkcj� pojedynczej zmiennej typu double.
 a - Lewy kraniec przedzia�u.
 b - Prawy kraniec przedzia�u.
 epsilon - Dok�adno�� zbie�no�ci (tolerancja b��du).
 x_dokladne - Wska�nik na tablic� dok�adnych rozwi�za�.
 ile_dokladnych - Liczba dok�adnych rozwi�za�.
 
 Warto�ci zwracane: struktura wyniki zawieraj�ca przybli�on� warto�� miejsca zerowego, tablic� kolejnych przybli�e�, tablic� b��d�w przybli�e� wzgl�dem dok�adnych rozwi�za�, liczb� wykonanych iteracji
 Przyk�ad u�ycia:
 auto f = [](double x) { return x*x - 2; };
  double dokladne[] = {1.414213562};
  Wyniki wynik = falszywa_linia(f, 0, 2, 1e-6, dokladne, 1);
  std::cout << "Pierwiastek: " << wynik.wynik << std::endl;
 */