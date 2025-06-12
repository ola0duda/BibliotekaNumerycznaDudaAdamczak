#include <iostream>
#include <cassert>
#include <cmath>
#include <vector>
#include <iomanip>
#include <functional>
#include <string>

// Definicja stałej PI dla zapewnienia kompatybilności, jeśli nie jest dostępna w <cmath>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// --- Nagłówki z biblioteki numerycznej ---
// Zakładamy, że projekt jest skonfigurowany tak, aby kompilator mógł znaleźć te pliki.
#include "../include/Aproksymacja.h"
#include "../include/CalkowanieNumeryczne.h"
#include "../include/Interpolacja.h"
#include "../include/RozwiazywanieRownanNieliniowych.h"
#include "../include/RozwiazywanieRownanRozniczkowych.h"
#include "../include/RozwiazywanieUkladowRownanLiniowych.h"
#include "../include/RóżniczkowanieNumeryczne.h"

// --- Narzędzia testowe ---

// Funkcja pomocnicza do porównywania wartości zmiennoprzecinkowych z zadaną tolerancją.
void check_close(double expected, double actual, const std::string& test_name, double epsilon = 1e-6) {
    if (std::isnan(expected) && std::isnan(actual)) return; // Obie wartości są NAN - test poprawny
    if (std::isnan(expected) || std::isnan(actual) || fabs(expected - actual) > epsilon) {
        std::cerr << "TEST FAILED: " << test_name << ". Oczekiwano " << expected << ", otrzymano " << actual << std::endl;
        exit(1);
    }
}

// Funkcja pomocnicza do tworzenia macierzy 2D.
double** create_matrix(int rows, int cols) {
    double** matrix = new double* [rows];
    for (int i = 0; i < rows; ++i) {
        matrix[i] = new double[cols](); // Inicjalizacja zerami
    }
    return matrix;
}

// Funkcja pomocnicza do zwalniania pamięci macierzy 2D.
void delete_matrix(double** matrix, int rows) {
    for (int i = 0; i < rows; ++i) {
        delete[] matrix[i];
    }
    delete[] matrix;
}

// --- Testy dla każdego modułu ---

// ===========================================
// Testy modułu: RóżniczkowanieNumeryczne
// ===========================================
double func_square(double x) { return x * x; }       // f'(x) = 2x
double func_const(double x) { return 10.0; }        // f'(x) = 0
double func_sin(double x) { return sin(x); }         // f'(x) = cos(x)

void test_rozniczkowanie() {
    std::cout << "Testowanie RóżniczkowanieNumeryczne..." << std::endl;

    // Test 1: Pochodna funkcji kwadratowej
    check_close(4.0, df_numeryczne(2.0, func_square), "df_numeryczne(x^2 dla x=2)");
    // Test 2: Pochodna funkcji stałej
    check_close(0.0, df_numeryczne(5.0, func_const), "df_numeryczne(const dla x=5)");
    // Test 3: Pochodna sinusa
    check_close(cos(M_PI / 3.0), df_numeryczne(M_PI / 3.0, func_sin), "df_numeryczne(sin(x) dla x=pi/3)");
    // Test 4: Przypadek błędny (dzielenie przez zero wewnątrz funkcji f)
    check_close(NAN, df_numeryczne(-2.0, func_square), "df_numeryczne - blad");

    std::cout << "PASSED!" << std::endl;
}

// ====================================================
// Testy modułu: RozwiazywanieRownanRozniczkowych
// ====================================================
// Równanie 1: y' = y, y(0) = 1. Rozwiązanie analityczne: y(t) = e^t
double ode_growth(double y) { return y; }
// Równanie 2: y' = -y, y(0) = 1. Rozwiązanie analityczne: y(t) = e^(-t)
double ode_decay(double y) { return -y; }

void test_rownania_rozniczkowe() {
    std::cout << "Testowanie RozwiazywanieRownanRozniczkowych..." << std::endl;
    const int N = 100;
    const double T0 = 1.0;
    const double CZAS = 1.0;
    double t[N + 1], y[N + 1];

    // Test 1: Wzrost wykładniczy
    double exact_growth = exp(1.0);
    euler(N, T0, t, y, N + 1, ode_growth, CZAS);
    // POPRAWKA: Zwiększona tolerancja dla metody Eulera ze względu na jej niski rząd dokładności.
    check_close(exact_growth, y[N], "euler (wzrost)", 1.5e-2);
    runge_kutta(N, T0, t, y, ode_growth, CZAS);
    check_close(exact_growth, y[N], "runge_kutta (wzrost)", 1e-7);

    // Test 2: Zanik wykładniczy
    double exact_decay = exp(-1.0);
    heun(N, T0, t, y, ode_decay, CZAS);
    check_close(exact_decay, y[N], "heun (zanik)", 1e-4);
    midpoint(N, T0, t, y, ode_decay, CZAS);
    check_close(exact_decay, y[N], "midpoint (zanik)", 1e-4);

    std::cout << "PASSED!" << std::endl;
}


// =================================================
// Testy modułu: RozwiazywanieRownanNieliniowych
// =================================================
double func_poly(double x) { return x * x - 9.0; } // Pierwiastki: 3, -3
double func_poly_deriv(double x) { return 2.0 * x; }
double func_no_root(double x) { return x * x + 9.0; } // Brak pierwiastków rzeczywistych

void test_rownania_nieliniowe() {
    std::cout << "Testowanie RozwiazywanieRownanNieliniowych..." << std::endl;
    const double epsilon = 1e-7;
    double exact_roots[] = { 3.0 };

    // Test 1: Poprawne znajdowanie pierwiastka
    Wyniki res_bisekcja = bisekcja(func_poly, 0, 5, epsilon);
    check_close(3.0, res_bisekcja.wynik, "bisekcja (poprawny)");

    Wyniki res_newton = newton(func_poly, func_poly_deriv, 5, epsilon);
    check_close(3.0, res_newton.wynik, "newton (poprawny)");

    // Test 2: Przypadki błędne
    Wyniki res_bisekcja_fail = bisekcja(func_no_root, 0, 5, epsilon);
    check_close(NAN, res_bisekcja_fail.wynik, "bisekcja (brak pierwiastka)");

    // Pochodna równa zero w punkcie startowym x=0 dla func_poly
    Wyniki res_newton_fail = newton(func_poly, func_poly_deriv, 0, epsilon);
    check_close(NAN, res_newton_fail.wynik, "newton (pochodna zero)");

    // Test 3: Inne metody
    Wyniki res_sieczna = sieczna(func_poly, 1, 5, epsilon);
    check_close(3.0, res_sieczna.wynik, "sieczna");

    Wyniki res_falszywa_linia = falszywa_linia(func_poly, 0, 5, epsilon, exact_roots, 1);
    check_close(3.0, res_falszywa_linia.wynik, "falszywa_linia");

    // Test 4: Funkcja licząca błąd
    check_close(0.0001, min_blad(3.0001, exact_roots, 1), "min_blad", 1e-5);
    check_close(NAN, min_blad(NAN, exact_roots, 1), "min_blad (wynik NAN)");

    std::cout << "PASSED!" << std::endl;
}

// ==============================
// Testy modułu: Interpolacja
// ==============================
void test_interpolacja() {
    std::cout << "Testowanie Interpolacja..." << std::endl;

    // Zestaw 1: y = 2x, punkty (1,2), (2,4), (3,6)
    double x1[] = { 1, 2, 3 };
    double y1[] = { 2, 4, 6 };
    check_close(5.0, lagrange(2.5, 3, x1, y1, 1), "lagrange (liniowa)");

    // Zestaw 2: y = x^2, punkty (1,1), (3,9), (4,16)
    double x2[] = { 1, 3, 4 };
    double y2[] = { 1, 9, 16 };
    check_close(4.0, lagrange(2.0, 3, x2, y2, 1), "lagrange (kwadratowa)");

    // Testy schematu Hornera i postaci naturalnej.
    // Konwencja: tablica wspolczynnikow {a0, a1, a2, ...}
    // odpowiada wielomianowi P(x) = a0 + a1*x + a2*x^2 + ...
    double coeffs[] = { 1, 2, 3 }; // Reprezentuje wielomian P(x) = 1 + 2x + 3x^2

    // Test dla x=2. P(2) = 1 + 2(2) + 3(2^2) = 1 + 4 + 12 = 17
    check_close(17.0, horner(2.0, coeffs, 3), "horner (x=2)");
    // Test dla x=0. P(0) = 1 + 2(0) + 3(0^2) = 1
    check_close(1.0, horner(0.0, coeffs, 3), "horner (x=0)", 1e-9);

    // Testy interpolacji Newtona
    int n_newton = 3;
    double** f = create_matrix(n_newton, n_newton);
    ilorazy(f, x2, y2, n_newton);
    double* newton_coeffs = new double[n_newton];
    for (int i = 0; i < n_newton; ++i) newton_coeffs[i] = f[i][i];

    check_close(4.0, newton(2.0, newton_coeffs, x2, n_newton), "newton (x=2)");
    check_close(25.0, newton(5.0, newton_coeffs, x2, n_newton), "newton (x=5)");

    delete[] newton_coeffs;
    delete_matrix(f, n_newton);
    std::cout << "PASSED!" << std::endl;
}


// =========================================
// Testy modułu: CalkowanieNumeryczne
// =========================================
// Zmienna globalna `ai` wymagana przez funkcję `wielomian` z biblioteki.
// W testach jest to obsłużone lokalnie.

void test_calkowanie() {
    std::cout << "Testowanie CalkowanieNumeryczne..." << std::endl;

    // Test 1: Całka z x^2 na [0, 1] = 1/3
    double interval[] = { 0.0, 1.0 };
    double poly_coeffs[] = { 0, 0, 1 }; // x^2
    double exact1 = 1.0 / 3.0;
    check_close(exact1, calka_prostokat(interval, poly_coeffs, 2, 10000), "calka_prostokat (x^2)", 1e-4);
    check_close(exact1, calka_trapez(interval, poly_coeffs, 2, 10000), "calka_trapez (x^2)", 1e-6);
    check_close(exact1, calka_simpson(interval, poly_coeffs, 2, 5000), "calka_simpson (x^2)", 1e-9);

    // Test 2: Całka z sin(x) na [-pi, pi] = 0
    double interval_sym[] = { -M_PI, M_PI };
    double exact2 = 0.0;
    check_close(exact2, calka_prostokat_funkcja(interval_sym, func_sin, 1000), "prostokat (sin(x) na [-pi,pi])", 1e-3);
    check_close(exact2, calka_trapez_funkcja(interval_sym, func_sin, 1000), "trapez (sin(x) na [-pi,pi])", 1e-5);
    check_close(exact2, calka_simpson_funkcja(interval_sym, func_sin, 500), "simpson (sin(x) na [-pi,pi])", 1e-7);

    // Test 3: Kwadratura Gaussa-Legendre'a
    auto func_x_cubed = [](double x) { return x * x * x; };
    check_close(0.0, kwadraturaGL1(func_x_cubed, -1.0, 1.0, 4, 100), "kwadraturaGL1 (x^3)");
    auto func_exp = [](double x) { return exp(x); };
    check_close(exp(1.0) - exp(0.0), kwadraturaGL1(func_exp, 0.0, 1.0, 4, 1000), "kwadraturaGL1 (exp(x))");

    std::cout << "PASSED!" << std::endl;
}

// ==========================================================
// Testy modułu: RozwiazywanieUkladowRownanLiniowych
// ==========================================================
void test_uklady_rownan() {
    std::cout << "Testowanie RozwiazywanieUkladowRownanLiniowych..." << std::endl;

    // Test 1: Układ 2x2
    const int n2 = 2;
    double** A2 = create_matrix(n2, n2);
    double b2[] = { 9, 2 };
    double x2[n2];
    A2[0][0] = 3; A2[0][1] = 2;
    A2[1][0] = 1; A2[1][1] = 2;
    // Rozwiązanie: x = 3.5, y = -0.75
    double** L2 = create_matrix(n2, n2);
    double** U2 = create_matrix(n2, n2);
    rozkladLU(A2, L2, U2, b2, n2);
    double z2[n2], x_lu2[n2];
    rozwiazanieLz(L2, b2, z2, n2);
    rozwiazanieUx(U2, z2, x_lu2, n2);
    check_close(3.5, x_lu2[0], "LU 2x2 - x[0]");
    check_close(-0.75, x_lu2[1], "LU 2x2 - x[1]");

    // Test 2: Układ 3x3
    const int n3 = 3;
    double** A3 = create_matrix(n3, n3);
    // POPRAWKA: Użycie poprawnego wektora b dla oczekiwanego rozwiązania.
    double b3[] = { 9, -13, -5 };
    A3[0][0] = 2; A3[0][1] = 1; A3[0][2] = -1;
    A3[1][0] = -3; A3[1][1] = -1; A3[1][2] = 2;
    A3[2][0] = -2; A3[2][1] = 1; A3[2][2] = 2;
    // Rozwiązanie: x=2, y=3, z=-2
    double** L3 = create_matrix(n3, n3);
    double** U3 = create_matrix(n3, n3);
    double b3_copy[] = { 9, -13, -5 };
    rozkladLU(A3, L3, U3, b3_copy, n3);
    double z3[n3], x_lu3[n3];
    rozwiazanieLz(L3, b3_copy, z3, n3);
    rozwiazanieUx(U3, z3, x_lu3, n3);
    check_close(2.0, x_lu3[0], "LU 3x3 - x[0]");
    check_close(3.0, x_lu3[1], "LU 3x3 - x[1]");
    check_close(-2.0, x_lu3[2], "LU 3x3 - x[2]");

    delete_matrix(A2, n2); delete_matrix(L2, n2); delete_matrix(U2, n2);
    delete_matrix(A3, n3); delete_matrix(L3, n3); delete_matrix(U3, n3);
    std::cout << "PASSED!" << std::endl;
}

// ==============================
// Testy modułu: Aproksymacja
// ==============================
void test_aproksymacja() {
    std::cout << "Testowanie Aproksymacja..." << std::endl;

    // Test 1: Funkcja bazowa
    check_close(8.0, baza(2.0, 3), "baza (dodatnia podstawa)");
    check_close(1.0, baza(123.0, 0), "baza (wykładnik zero)");
    check_close(-27.0, baza(-3.0, 3), "baza (ujemna podstawa)");

    // Test 2: Kwadratura GL2
    auto func_x_cubed = [](double x) { return x * x * x; };
    check_close(0.0, kwadraturaGL2(func_x_cubed, -1.0, 1.0), "kwadraturaGL2 (nieparzysta)");
    auto func_exp = [](double x) { return exp(x); };
    check_close(exp(1.0) - exp(0.0), kwadraturaGL2(func_exp, 0.0, 1.0), "kwadraturaGL2 (exp(x))");

    // Test 3: Eliminacja Gaussa z rozwiązaniem (gauss2)
    const int n = 2;
    double** A = create_matrix(n, n);
    double b[] = { 9, 2 };
    double x[n];
    A[0][0] = 3; A[0][1] = 2;
    A[1][0] = 1; A[1][1] = 2;
    gauss2(A, b, x, n);
    check_close(3.5, x[0], "gauss2 - x[0]");
    check_close(-0.75, x[1], "gauss2 - x[1]");

    // Uwaga: Główna funkcja `aproksymacja` nie jest testowana bezpośrednio,
    // ponieważ jej wynikiem jest wydruk na standardowe wyjście.
    // Testowane są jednak jej kluczowe komponenty (`kwadraturaGL2`, `gauss2`).

    delete_matrix(A, n);
    std::cout << "PASSED!" << std::endl;
}

// =============================================
// Główna funkcja uruchamiająca wszystkie testy
// =============================================
int main() {
    // Ustawienie polskiego locale dla komunikatów C++ (opcjonalne)
    setlocale(LC_ALL, "pl_PL.UTF-8");

    try {
        test_rozniczkowanie();
        test_rownania_rozniczkowe();
        test_rownania_nieliniowe();
        test_interpolacja();
        test_calkowanie();
        test_uklady_rownan();
        test_aproksymacja();

        std::cout << "\n---------------------------------" << std::endl;
        std::cout << "WSZYSTKIE TESTY ZAKONCZONE SUKCESEM!" << std::endl;
        std::cout << "---------------------------------" << std::endl;

    }
    catch (const std::exception& e) {
        std::cerr << "Wystapil krytyczny blad: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
