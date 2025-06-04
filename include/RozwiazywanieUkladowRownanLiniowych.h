void pokazMacierz1(double** A, double* b, int n); /**
Opis dzia�ania: Wy�wietla rozszerzon� macierz wsp�czynnik�w A oraz wektora wyraz�w wolnych b.
 Opis argument�w:
A - Macierz wsp�czynnik�w (n x n).
b - Wektor wyraz�w wolnych (n).
n - Rozmiar uk�adu r�wna�.
 
 Warto�ci zwracane: Brak (wy�wietla macierz)
 Przyk�ad u�ycia:
 pokazMacierz1(A, b, n);
 */

void gauss(double** A, double* b, int n); /**
Opis dzia�ania: Rozwi�zuje uk�ad r�wna� metod� eliminacji Gaussa z cz�ciowym wyborem elementu g��wnego (pivoting).
 Opis argument�w:
A - Macierz wsp�czynnik�w (n x n), modyfikowana w trakcie oblicze�.
b - Wektor wyraz�w wolnych (n), modyfikowany w trakcie oblicze�.
n - Rozmiar uk�adu r�wna�.
 
 Warto�ci zwracane: Modyfikuje macierz A do postaci tr�jk�tnej g�rnej oraz przekszta�ca wektor b.
 Przyk�ad u�ycia:
gauss(A, b, n);
 */

void trojkat(double** A, double* b, double* x, int n); /**
Opis dzia�ania: Rozwi�zuje uk�ad r�wna� w postaci tr�jk�tnej g�rnej (Ax = b) przez podstawianie wsteczne.
 Opis argument�w:
A - Macierz tr�jk�tna g�rna (n x n).
b - Wektor prawej strony r�wnania (n).
x - Wektor rozwi�zania (n) � wynik dzia�ania funkcji.
n - Rozmiar uk�adu.
 
 Warto�ci zwracane: Brak (wyniki zapisuje w x)
 Przyk�ad u�ycia:
 trojkat(A, b, x, n);
 */

void pokazMacierz2(double** A, int n); /**
Opis dzia�ania: Wy�wietla macierz A o wymiarze n x n.
 Opis argument�w:
A - Macierz do wy�wietlenia.
n - Rozmiar macierzy.
 
 Warto�ci zwracane: Brak (wyswietla macierz)
 Przyk�ad u�ycia:
pokazMacierz2(L, n);
 */

void rozkladLU(double** A, double** L, double** U, double* b, int n); /**
Opis dzia�ania: Wykonuje rozk�ad macierzy A na macierze doln� L i g�rn� U (LU-dekompozycja) z pivotingiem.
 Opis argument�w:
A - Macierz wsp�czynnik�w (n x n), modyfikowana podczas rozk�adu.
L - Macierz dolna (n x n) � wynik dzia�ania funkcji.
U - Macierz g�rna (n x n) � wynik dzia�ania funkcji.
b - Wektor wyraz�w wolnych, modyfikowany przy pivotingu.
n - Rozmiar uk�adu.
 
 Warto�ci zwracane: Brak (Funkcja dzieli macierz A na iloczyn dw�ch macierzy: L i U tak, �e A = L * U)
 Przyk�ad u�ycia:
 rozkladLU(A, L, U, b, n);
 */

void rozwiazanieLz(double** L, double* b, double* z, int n); /**
Opis dzia�ania: Rozwi�zuje uk�ad r�wna� Lz = b (macierz L � dolna tr�jk�tna) metod� podstawiania w prz�d.
 Opis argument�w:
L - Macierz dolna (n x n).
b - Wektor prawej strony r�wnania (n).
z - Wektor rozwi�zania po�redniego (n) � wynik dzia�ania funkcji.
n - Rozmiar uk�adu.
 
 Warto�ci zwracane: Brak (wyniki zapisuje w z)
 Przyk�ad u�ycia:
 rozwiazanieLz(L, b, z, n);
 */

void rozwiazanieUx(double** U, double* z, double* x, int n); /**
Opis dzia�ania: Rozwi�zuje uk�ad r�wna� Ux = z (macierz U � g�rna tr�jk�tna) metod� podstawiania wstecz.
 Opis argument�w:
U - Macierz g�rna (n x n).
z - Wektor po�redni (n), wynik rozwi�zania Lz = b.
x - Wektor rozwi�zania ko�cowego (n) � wynik dzia�ania funkcji.
n - Rozmiar uk�adu.
 
 Warto�ci zwracane: Brak (wyniki zapisuje w x)
 Przyk�ad u�ycia:
 rozwiazanieUx(U, z, x, n);
 */