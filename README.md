#BibliotekaNumerycznaDudaAdamczak

##Opis biblioteki

Nasza biblioteka numeryczna to biblioteka statyczna napisana w języku C++ zawierająca implementacje podstawowych metod obliczeń numerycznych, w tym:
- rozwiązywanie układów równań liniowych,
- interpolację - interpolacja Lagrange'a, wielomian interpolacyjny Newtona
- aproksymację
- całkowanie numeryczne - metoda prostokątów, metoda trapezów, metoda Simpsona, metoda Gaussa-Legendre'a
- różniczkowanie numeryczne - metoda różnic centralnych
- rozwiązywanie równań nieliniowych - metoda Newtona, siecznych, bisekcji, fałszywej linii
- rozwiązywanie równań różniczkowych - metoda Eulera, metoda Heuna, metoda punktu środkowego, metoda Rungego-Kutty 4. rzędu

Projekt jest zorganizowany w sposób przejrzysty i modułowy, struktura:
BibliotekaNumerycznaDudaAdamczak/
├──include/ # pliki nagłówkowe (*.h)
├──src/ # pliki źródłowe (*.cpp)
├──tests/ # testy jednostkowe
├──examples/

Biblioteka jest kompatybilna z Microsoft Visual Studio.
Projekt edukacyjny przeznaczony do nauki.

----------------------------------------------------------

##Instrukcja instalacji

1.Sklonuj repozytorium:
    git clone https://github.com/ola0duda/BibliotekaNumerycznaDudaAdamczak.git

2.Otwórz projekt w Visual Studio:
    - Otwórz plik '.sln' znajdujący się w katalogu głównym.
    - Upewnij się, że aktywna konfiguracja to 'Release' lub 'Debug', w zależności od potrzeb.

3. Zbuduj bibliotekę:
    - W Visual Studio wybierz 'Build' -> 'Build Solution' (lub naciśnij 'Ctrl + Shift + B').

4. Dodaj bibliotekę do swojego projektu:
    - W 'C/C++ -> General -> Additional Include Directiories' i 'Linker -> General -> Additional Include Directories' dołącz katalog 'include/' do ścieżek nagłówków.
    - W 'Linker -> General -> Additional Library Directories' dodaj katalog 'src/' do ścieżek bibliotek.
    - Dołącz bibliotekę dodając nazwę final.lib w 'Linker -> Input -> Additional Dependencies'.

----------------------------------------------------------

##Przykłady użycia

###Przykład 1: Różniczkowanie numeryczne
double funkcja(double x) {
    return x * x;
}

double wynik = df_numeryczne(2.0, funkcja);
//wynik = 4.0

###Przykład 2: Rozwiązywanie równań różniczkowych metodą Eulera
double rownanie(double x) {
    return -0.3 * x;
}

double y[100], z[100];
int n = euler(10, 100.0, y, z, 100, rownanie, 100.0);

###Przykład 3: Rozwiązywanie układu równań metodą Gaussa
double* b = new double[n];
double** A = new double*[n];
for (int i = 0; i < n; ++i) {
    A[i] = new double[n];
    //uzupełnij A[i][j] i b[i]
}
gauss(A, b, n);
trojkat(A, b, x, n) //do uzyskania rozwiązań
