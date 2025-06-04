double euler(int N, double T0_val, double t_tab[], double T_tab[], int max_rozmiar_tab, double (*rownanie)(double), double CZAS); /**
Opis dzia³ania: Funkcja rozwi¹zuj¹ca równanie ró¿niczkowe metod¹ Eulera.
 Opis argumentów:
N - liczba kroków ca³kowania
T0_val - wartoœæ pocz¹tkowa funkcji T(0)
t_tab[] - tablica na kolejne wartoœci czasu
T_tab[] - tablica na kolejne wartoœci funkcji T
max_rozmiar_tab - maksymalny rozmiar tablic t_tab i T_tab (niewykorzystywany w tej implementacji)
double (*rownanie)(double) - wskaŸnik na funkcjê opisuj¹c¹ pochodn¹ (np. dT/dt = f(T))
CZAS - ca³kowity czas symulacji
 
 Wartoœci zwracane: Liczba wype³nionych elementów w tablicy (czyli N + 1)
 Przyk³ad u¿ycia:
 double wynik = euler(100, 20.0, t, T, 101, funkcjaPochodna, 1000.0);
 */

int heun(int N, double T0, double t[], double T[], double (*rownanie)(double), double CZAS); /**
Opis dzia³ania: Funkcja rozwi¹zuj¹ca równanie ró¿niczkowe metod¹ Heuna (œredniego nachylenia).
 Opis argumentów:
N - liczba kroków ca³kowania
T0 - wartoœæ pocz¹tkowa funkcji T(0)
t[] - tablica na kolejne wartoœci czasu
T[] - tablica na kolejne wartoœci funkcji T
double (*rownanie)(double) - funkcja opisuj¹ca pochodn¹
CZAS - ca³kowity czas symulacji
 
 Wartoœci zwracane: Liczba wype³nionych elementów (N + 1)
 Przyk³ad u¿ycia:
 int rozmiar = heun(100, 20.0, t, T, funkcjaPochodna, 1000.0);
 */

int midpoint(int N, double T0, double t[], double T[], double (*rownanie)(double), double CZAS); /**
Opis dzia³ania: Funkcja rozwi¹zuj¹ca równanie ró¿niczkowe metod¹ punktu œrodkowego.
 Opis argumentów:
N - liczba kroków ca³kowania
T0 - wartoœæ pocz¹tkowa funkcji T(0)
t[] - tablica na kolejne wartoœci czasu
T[] - tablica na kolejne wartoœci funkcji T
double (*rownanie)(double) - wskaŸnik na funkcjê opisuj¹c¹ pochodn¹
CZAS - ca³kowity czas symulacji
 
 Wartoœci zwracane: Liczba wype³nionych elementów (N + 1)
 Przyk³ad u¿ycia:
 int rozmiar = midpoint(100, 20.0, t, T, funkcjaPochodna, 1000.0);
 */

int runge_kutta(int N, double T0, double t[], double T[], double (*rownanie)(double), double CZAS); /**
Opis dzia³ania: Funkcja rozwi¹zuj¹ca równanie ró¿niczkowe metod¹ Rungego-Kutty 4. rzêdu.
 Opis argumentów:
N - liczba kroków ca³kowania
T0 - wartoœæ pocz¹tkowa funkcji T(0)
t[] - tablica na kolejne wartoœci czasu
T[] - tablica na kolejne wartoœci funkcji T
double (*rownanie)(double) - wskaŸnik na funkcjê opisuj¹c¹ pochodn¹
CZAS - ca³kowity czas symulacji
 
 Wartoœci zwracane: Liczba wype³nionych elementów (N + 1)
 Przyk³ad u¿ycia:
 int rozmiar = runge_kutta(100, 20.0, t, T, funkcjaPochodna, 1000.0);
 */