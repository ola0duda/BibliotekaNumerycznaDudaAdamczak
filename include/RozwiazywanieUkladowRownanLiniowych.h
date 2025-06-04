void pokazMacierz1(double** A, double* b, int n); /**
Opis dzia³ania: Wyœwietla rozszerzon¹ macierz wspó³czynników A oraz wektora wyrazów wolnych b.
 Opis argumentów:
A - Macierz wspó³czynników (n x n).
b - Wektor wyrazów wolnych (n).
n - Rozmiar uk³adu równañ.
 
 Wartoœci zwracane: Brak (wyœwietla macierz)
 Przyk³ad u¿ycia:
 pokazMacierz1(A, b, n);
 */

void gauss(double** A, double* b, int n); /**
Opis dzia³ania: Rozwi¹zuje uk³ad równañ metod¹ eliminacji Gaussa z czêœciowym wyborem elementu g³ównego (pivoting).
 Opis argumentów:
A - Macierz wspó³czynników (n x n), modyfikowana w trakcie obliczeñ.
b - Wektor wyrazów wolnych (n), modyfikowany w trakcie obliczeñ.
n - Rozmiar uk³adu równañ.
 
 Wartoœci zwracane: Modyfikuje macierz A do postaci trójk¹tnej górnej oraz przekszta³ca wektor b.
 Przyk³ad u¿ycia:
gauss(A, b, n);
 */

void trojkat(double** A, double* b, double* x, int n); /**
Opis dzia³ania: Rozwi¹zuje uk³ad równañ w postaci trójk¹tnej górnej (Ax = b) przez podstawianie wsteczne.
 Opis argumentów:
A - Macierz trójk¹tna górna (n x n).
b - Wektor prawej strony równania (n).
x - Wektor rozwi¹zania (n) – wynik dzia³ania funkcji.
n - Rozmiar uk³adu.
 
 Wartoœci zwracane: Brak (wyniki zapisuje w x)
 Przyk³ad u¿ycia:
 trojkat(A, b, x, n);
 */

void pokazMacierz2(double** A, int n); /**
Opis dzia³ania: Wyœwietla macierz A o wymiarze n x n.
 Opis argumentów:
A - Macierz do wyœwietlenia.
n - Rozmiar macierzy.
 
 Wartoœci zwracane: Brak (wyswietla macierz)
 Przyk³ad u¿ycia:
pokazMacierz2(L, n);
 */

void rozkladLU(double** A, double** L, double** U, double* b, int n); /**
Opis dzia³ania: Wykonuje rozk³ad macierzy A na macierze doln¹ L i górn¹ U (LU-dekompozycja) z pivotingiem.
 Opis argumentów:
A - Macierz wspó³czynników (n x n), modyfikowana podczas rozk³adu.
L - Macierz dolna (n x n) – wynik dzia³ania funkcji.
U - Macierz górna (n x n) – wynik dzia³ania funkcji.
b - Wektor wyrazów wolnych, modyfikowany przy pivotingu.
n - Rozmiar uk³adu.
 
 Wartoœci zwracane: Brak (Funkcja dzieli macierz A na iloczyn dwóch macierzy: L i U tak, ¿e A = L * U)
 Przyk³ad u¿ycia:
 rozkladLU(A, L, U, b, n);
 */

void rozwiazanieLz(double** L, double* b, double* z, int n); /**
Opis dzia³ania: Rozwi¹zuje uk³ad równañ Lz = b (macierz L – dolna trójk¹tna) metod¹ podstawiania w przód.
 Opis argumentów:
L - Macierz dolna (n x n).
b - Wektor prawej strony równania (n).
z - Wektor rozwi¹zania poœredniego (n) – wynik dzia³ania funkcji.
n - Rozmiar uk³adu.
 
 Wartoœci zwracane: Brak (wyniki zapisuje w z)
 Przyk³ad u¿ycia:
 rozwiazanieLz(L, b, z, n);
 */

void rozwiazanieUx(double** U, double* z, double* x, int n); /**
Opis dzia³ania: Rozwi¹zuje uk³ad równañ Ux = z (macierz U – górna trójk¹tna) metod¹ podstawiania wstecz.
 Opis argumentów:
U - Macierz górna (n x n).
z - Wektor poœredni (n), wynik rozwi¹zania Lz = b.
x - Wektor rozwi¹zania koñcowego (n) – wynik dzia³ania funkcji.
n - Rozmiar uk³adu.
 
 Wartoœci zwracane: Brak (wyniki zapisuje w x)
 Przyk³ad u¿ycia:
 rozwiazanieUx(U, z, x, n);
 */