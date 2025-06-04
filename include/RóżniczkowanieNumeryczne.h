double df_numeryczne(double x, double (*f)(double)); /**
Opis dzia³ania: Oblicza pochodn¹ funkcji f w punkcie x metod¹ ró¿nic centralnych.
 Opis argumentów:
x - Punkt, w którym liczona jest pochodna.
f - WskaŸnik do funkcji f(x), dla której liczona jest pochodna.
 
 Wartoœci zwracane: Przybli¿ona wartoœæ pochodnej funkcji f w punkcie x. Zwraca NaN, jeœli obliczenie jest niemo¿liwe (np. nieokreœlona wartoœæ funkcji).
 Metoda u¿ywa wzoru: f'(x) = (f(x + h) - f(x - h)) / (2h), gdzie h = 1e-6
 Przyk³ad u¿ycia:
 double funkcja(double x) {
      return x * x;
 }
 double pochodna = df_numeryczne(2.0, funkcja);  // Wynik = 4.0
 */