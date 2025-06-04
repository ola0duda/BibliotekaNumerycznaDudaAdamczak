double euler(int N, double T0_val, double t_tab[], double T_tab[], int max_rozmiar_tab, double (*rownanie)(double), double CZAS); /**
Opis dzia�ania: Funkcja rozwi�zuj�ca r�wnanie r�niczkowe metod� Eulera.
 Opis argument�w:
N - liczba krok�w ca�kowania
T0_val - warto�� pocz�tkowa funkcji T(0)
t_tab[] - tablica na kolejne warto�ci czasu
T_tab[] - tablica na kolejne warto�ci funkcji T
max_rozmiar_tab - maksymalny rozmiar tablic t_tab i T_tab (niewykorzystywany w tej implementacji)
double (*rownanie)(double) - wska�nik na funkcj� opisuj�c� pochodn� (np. dT/dt = f(T))
CZAS - ca�kowity czas symulacji
 
 Warto�ci zwracane: Liczba wype�nionych element�w w tablicy (czyli N + 1)
 Przyk�ad u�ycia:
 double wynik = euler(100, 20.0, t, T, 101, funkcjaPochodna, 1000.0);
 */

int heun(int N, double T0, double t[], double T[], double (*rownanie)(double), double CZAS); /**
Opis dzia�ania: Funkcja rozwi�zuj�ca r�wnanie r�niczkowe metod� Heuna (�redniego nachylenia).
 Opis argument�w:
N - liczba krok�w ca�kowania
T0 - warto�� pocz�tkowa funkcji T(0)
t[] - tablica na kolejne warto�ci czasu
T[] - tablica na kolejne warto�ci funkcji T
double (*rownanie)(double) - funkcja opisuj�ca pochodn�
CZAS - ca�kowity czas symulacji
 
 Warto�ci zwracane: Liczba wype�nionych element�w (N + 1)
 Przyk�ad u�ycia:
 int rozmiar = heun(100, 20.0, t, T, funkcjaPochodna, 1000.0);
 */

int midpoint(int N, double T0, double t[], double T[], double (*rownanie)(double), double CZAS); /**
Opis dzia�ania: Funkcja rozwi�zuj�ca r�wnanie r�niczkowe metod� punktu �rodkowego.
 Opis argument�w:
N - liczba krok�w ca�kowania
T0 - warto�� pocz�tkowa funkcji T(0)
t[] - tablica na kolejne warto�ci czasu
T[] - tablica na kolejne warto�ci funkcji T
double (*rownanie)(double) - wska�nik na funkcj� opisuj�c� pochodn�
CZAS - ca�kowity czas symulacji
 
 Warto�ci zwracane: Liczba wype�nionych element�w (N + 1)
 Przyk�ad u�ycia:
 int rozmiar = midpoint(100, 20.0, t, T, funkcjaPochodna, 1000.0);
 */

int runge_kutta(int N, double T0, double t[], double T[], double (*rownanie)(double), double CZAS); /**
Opis dzia�ania: Funkcja rozwi�zuj�ca r�wnanie r�niczkowe metod� Rungego-Kutty 4. rz�du.
 Opis argument�w:
N - liczba krok�w ca�kowania
T0 - warto�� pocz�tkowa funkcji T(0)
t[] - tablica na kolejne warto�ci czasu
T[] - tablica na kolejne warto�ci funkcji T
double (*rownanie)(double) - wska�nik na funkcj� opisuj�c� pochodn�
CZAS - ca�kowity czas symulacji
 
 Warto�ci zwracane: Liczba wype�nionych element�w (N + 1)
 Przyk�ad u�ycia:
 int rozmiar = runge_kutta(100, 20.0, t, T, funkcjaPochodna, 1000.0);
 */