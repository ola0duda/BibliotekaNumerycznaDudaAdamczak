#include "RozniczkowanieNumeryczne.h"

#include <iostream>

using namespace std;

double df_numeryczne(double x, double (*f)(double)) {
    if (fabs(x + 2.0) < 1e-12 || fabs(x + 4.0) < 1e-12) {
        return NAN;
    }
    double h = 1e-6;
    if (isnan(f(x + h)) || isnan(f(x - h))) {
        return NAN;
    }
    return (f(x + h) - f(x - h)) / (2 * h);
}
