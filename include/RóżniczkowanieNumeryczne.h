double df_numeryczne(double x, double (*f)(double)); /**
Opis dzia�ania: Oblicza pochodn� funkcji f w punkcie x metod� r�nic centralnych.
 Opis argument�w:
x - Punkt, w kt�rym liczona jest pochodna.
f - Wska�nik do funkcji f(x), dla kt�rej liczona jest pochodna.
 
 Warto�ci zwracane: Przybli�ona warto�� pochodnej funkcji f w punkcie x. Zwraca NaN, je�li obliczenie jest niemo�liwe (np. nieokre�lona warto�� funkcji).
 Metoda u�ywa wzoru: f'(x) = (f(x + h) - f(x - h)) / (2h), gdzie h = 1e-6
 Przyk�ad u�ycia:
 double funkcja(double x) {
      return x * x;
 }
 double pochodna = df_numeryczne(2.0, funkcja);  // Wynik = 4.0
 */