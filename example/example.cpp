#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

// --- Nagłówki z biblioteki numerycznej ---
#include "../include/RóżniczkowanieNumeryczne.h"
#include "../include/Interpolacja.h"
#include "../include/CalkowanieNumeryczne.h"
#include "../include/RozwiazywanieRownanNieliniowych.h"
#include "../include/RozwiazywanieUkladowRownanLiniowych.h"
#include "../include/RozwiazywanieRownanRozniczkowych.h"
#include "../include/Aproksymacja.h"

// Definicja stałej PI dla zapewnienia kompatybilności
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// --- Funkcje pomocnicze do demonstracji ---

// Funkcja do wypisywania wektora
void print_vector(const double* vec, int size, const std::string& name) {
    std::cout << name << " = [ ";
    for (int i = 0; i < size; ++i) {
        std::cout << vec[i] << (i == size - 1 ? " " : ", ");
    }
    std::cout << "]" << std::endl;
}

// Funkcja do wypisywania macierzy
void print_matrix(double** M, int n, const std::string& name) {
    std::cout << "Macierz " << name << ":" << std::endl;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            std::cout << std::setw(8) << std::fixed << std::setprecision(4) << M[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

// Funkcja pomocnicza do tworzenia macierzy 2D
double** create_matrix(int rows, int cols) {
    double** matrix = new double* [rows];
    for (int i = 0; i < rows; ++i) {
        matrix[i] = new double[cols]();
    }
    return matrix;
}

// Funkcja pomocnicza do zwalniania pamięci macierzy 2D
void delete_matrix(double** matrix, int rows) {
    for (int i = 0; i < rows; ++i) {
        delete[] matrix[i];
    }
    delete[] matrix;
}


int main() {
    setlocale(LC_ALL, "pl_PL.UTF-8");
    std::cout << std::fixed << std::setprecision(6);

    std::cout << "=====================================================" << std::endl;
    std::cout << "==      PRZYKŁADY UŻYCIA BIBLIOTEKI NUMERYCZNEJ    ==" << std::endl;
    std::cout << "=====================================================" << std::endl << std::endl;

    //--------------------------------------------------------------------
    // 1. RÓŻNICZKOWANIE NUMERYCZNE
    //--------------------------------------------------------------------
    std::cout << "--- 1. Różniczkowanie Numeryczne ---" << std::endl;
    auto funkcja_do_rozniczkowania = [](double x) { return x * x * x - 2 * x; }; // f(x) = x^3 - 2x
    double punkt = 2.0;
    double pochodna_analityczna = 3 * punkt * punkt - 2; // f'(x) = 3x^2 - 2 => f'(2) = 10
    double pochodna_numeryczna = df_numeryczne(punkt, funkcja_do_rozniczkowania);
    std::cout << "Funkcja: f(x) = x^3 - 2x" << std::endl;
    std::cout << "Pochodna w punkcie x = " << punkt << ":" << std::endl;
    std::cout << "  - Wynik analityczny: " << pochodna_analityczna << std::endl;
    std::cout << "  - Wynik numeryczny:  " << pochodna_numeryczna << std::endl << std::endl;

    //--------------------------------------------------------------------
    // 2. INTERPOLACJA
    //--------------------------------------------------------------------
    std::cout << "--- 2. Interpolacja ---" << std::endl;
    double x_nodes[] = { 0.0, M_PI / 2, M_PI, 3 * M_PI / 2 };
    double y_nodes[] = { 0.0, 1.0, 0.0, -1.0 }; // Węzły dla funkcji sin(x)
    int n_nodes = 4;
    double punkt_interpolacji = M_PI / 4; // sin(pi/4) = sqrt(2)/2 ~= 0.707106
    std::cout << "Interpolacja dla funkcji f(x) = sin(x) w punkcie x = PI/4" << std::endl;

    // a) Interpolacja Lagrange'a
    double wynik_lagrange = lagrange(punkt_interpolacji, n_nodes, x_nodes, y_nodes, 1);
    std::cout << "  - Metoda Lagrange'a: f(" << punkt_interpolacji << ") ~= " << wynik_lagrange << std::endl;

    // b) Interpolacja Newtona
    double** f_ilorazy = create_matrix(n_nodes, n_nodes);
    ilorazy(f_ilorazy, x_nodes, y_nodes, n_nodes);
    double* newton_coeffs = new double[n_nodes];
    for (int i = 0; i < n_nodes; ++i) newton_coeffs[i] = f_ilorazy[i][i];
    double wynik_newton = newton(punkt_interpolacji, newton_coeffs, x_nodes, n_nodes);
    std::cout << "  - Metoda Newtona:    f(" << punkt_interpolacji << ") ~= " << wynik_newton << std::endl;
    std::cout << "  - Wartość dokładna:  f(" << punkt_interpolacji << ") =  " << sin(punkt_interpolacji) << std::endl << std::endl;
    delete_matrix(f_ilorazy, n_nodes);
    delete[] newton_coeffs;

    //--------------------------------------------------------------------
    // 3. CAŁKOWANIE NUMERYCZNE
    //--------------------------------------------------------------------
    std::cout << "--- 3. Całkowanie Numeryczne ---" << std::endl;
    auto funkcja_do_calkowania = [](double x) { return x * x; }; // f(x) = x^2
    double przedzial[] = { 0.0, 1.0 };
    int n_calk = 1000;
    std::cout << "Całka z f(x) = x^2 na przedziale [0, 1]. Wynik analityczny: 1/3 ~= 0.333333" << std::endl;
    double calka_p = calka_prostokat_funkcja(przedzial, funkcja_do_calkowania, n_calk);
    double calka_t = calka_trapez_funkcja(przedzial, funkcja_do_calkowania, n_calk);
    double calka_s = calka_simpson_funkcja(przedzial, funkcja_do_calkowania, n_calk / 2);
    double calka_gl = kwadraturaGL1(funkcja_do_calkowania, przedzial[0], przedzial[1], 4, 100);
    std::cout << "  - Metoda prostokątów: " << calka_p << std::endl;
    std::cout << "  - Metoda trapezów:    " << calka_t << std::endl;
    std::cout << "  - Metoda Simpsona:    " << calka_s << std::endl;
    std::cout << "  - Kwadratura G-L:     " << calka_gl << std::endl << std::endl;

    //--------------------------------------------------------------------
    // 4. ROZWIĄZYWANIE RÓWNAŃ NIELINIOWYCH
    //--------------------------------------------------------------------
    std::cout << "--- 4. Rozwiązywanie Równań Nieliniowych ---" << std::endl;
    auto f_nielin = [](double x) { return x * x - 2; }; // f(x) = x^2 - 2, pierwiastek = sqrt(2) ~= 1.414213
    auto df_nielin = [](double x) { return 2 * x; }; // f'(x) = 2x
    std::cout << "Szukanie pierwiastka funkcji f(x) = x^2 - 2" << std::endl;
    Wyniki res_bis = bisekcja(f_nielin, 1.0, 2.0, 1e-7);
    Wyniki res_new = newton(f_nielin, df_nielin, 1.0, 1e-7);
    Wyniki res_sie = sieczna(f_nielin, 1.0, 2.0, 1e-7);
    std::cout << "  - Metoda bisekcji: x = " << res_bis.wynik << std::endl;
    std::cout << "  - Metoda Newtona:  x = " << res_new.wynik << std::endl;
    std::cout << "  - Metoda siecznych:x = " << res_sie.wynik << std::endl;
    std::cout << "  - Wartość dokładna:  x = " << sqrt(2.0) << std::endl << std::endl;

    //--------------------------------------------------------------------
    // 5. ROZWIĄZYWANIE UKŁADÓW RÓWNAŃ LINIOWYCH
    //--------------------------------------------------------------------
    std::cout << "--- 5. Rozwiązywanie Układów Równań Liniowych ---" << std::endl;
    const int n_rown = 3;
    double** A = create_matrix(n_rown, n_rown);
    double b[] = { 9, -13, -5 };
    A[0][0] = 2; A[0][1] = 1; A[0][2] = -1;
    A[1][0] = -3; A[1][1] = -1; A[1][2] = 2;
    A[2][0] = -2; A[2][1] = 1; A[2][2] = 2;
    std::cout << "Rozwiązywanie układu 3x3. Oczekiwane rozwiązanie: x = [2, 3, -2]" << std::endl;

    // a) Metoda dekompozycji LU
    double** L = create_matrix(n_rown, n_rown);
    double** U = create_matrix(n_rown, n_rown);
    double b_lu[] = { 9, -13, -5 };
    double z[n_rown], x_lu[n_rown];
    rozkladLU(A, L, U, b_lu, n_rown); // Ta funkcja wypisuje swoje kroki
    rozwiazanieLz(L, b_lu, z, n_rown);
    rozwiazanieUx(U, z, x_lu, n_rown);
    print_vector(x_lu, n_rown, "  - Rozwiązanie (metoda LU)");
    std::cout << std::endl;
    delete_matrix(A, n_rown); delete_matrix(L, n_rown); delete_matrix(U, n_rown);

    //--------------------------------------------------------------------
    // 6. ROZWIĄZYWANIE RÓWNAŃ RÓŻNICZKOWYCH
    //--------------------------------------------------------------------
    std::cout << "--- 6. Rozwiązywanie Równań Różniczkowych ---" << std::endl;
    auto rownanie_rozniczkowe = [](double T) { return -0.1 * (T - 20.0); }; // Chłodzenie obiektu: T'(t) = -k(T-T_otoczenia)
    double T0 = 100.0; // Temperatura początkowa
    double CZAS_R = 10.0; // Czas symulacji
    const int N_R = 50;
    double t_tab[N_R + 1], T_tab_euler[N_R + 1], T_tab_rk[N_R + 1];
    std::cout << "Symulacja chłodzenia obiektu z T(0)=100 do T_otoczenia=20" << std::endl;

    euler(N_R, T0, t_tab, T_tab_euler, N_R + 1, rownanie_rozniczkowe, CZAS_R);
    runge_kutta(N_R, T0, t_tab, T_tab_rk, rownanie_rozniczkowe, CZAS_R);
    // Rozwiązanie analityczne: T(t) = 20 + (100-20) * exp(-0.1*t)
    double T_analityczne = 20 + 80 * exp(-0.1 * CZAS_R);

    std::cout << "Temperatura po " << CZAS_R << " sekundach:" << std::endl;
    std::cout << "  - Metoda Eulera:       T(" << CZAS_R << ") = " << T_tab_euler[N_R] << std::endl;
    std::cout << "  - Metoda Rungego-Kutty: T(" << CZAS_R << ") = " << T_tab_rk[N_R] << std::endl;
    std::cout << "  - Wynik analityczny:   T(" << CZAS_R << ") = " << T_analityczne << std::endl << std::endl;

    //--------------------------------------------------------------------
    // 7. APROKSYMACJA
    //--------------------------------------------------------------------
    std::cout << "--- 7. Aproksymacja ---" << std::endl;
    std::cout << "Aproksymacja funkcji f(x) = exp(x) na przedziale [0, 1] wielomianem stopnia 3." << std::endl;
    std::cout << "Wyniki (współczynniki i błędy) zostaną wypisane przez funkcję aproksymacja:" << std::endl;
    // Funkcja `aproksymacja` sama wypisuje szczegółowe wyniki
    aproksymacja(0.0, 1.0, 4, exp);
    std::cout << std::endl;

    std::cout << "=====================================================" << std::endl;
    std::cout << "== ZAKOŃCZONO DEMONSTRACJĘ                         ==" << std::endl;
    std::cout << "=====================================================" << std::endl;

    return 0;
}
