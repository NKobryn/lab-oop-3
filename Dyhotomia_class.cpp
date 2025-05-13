#include "Dyhotomia_class.h"
#include <cmath>
#include <limits>

Dyhotomia_class::Dyhotomia_class() {}

Dyhotomia_class::~Dyhotomia_class() {}

void Dyhotomia_class::setVolumes(double vol_a, double vol_b) {
    a = vol_a;
    b = vol_b;
}

void Dyhotomia_class::setTolerance(double vol_eps) {
    eps = vol_eps;
}

double Dyhotomia_class::f(double x) {
    return x + sqrt(x) + cbrt(x) - 2.5;
}

int Dyhotomia_class::count(double &x) {
    double fa = f(a);
    double fb = f(b);

    if (fa * fb >= 0) return 1;

    double c;
    while ((b - a) >= eps) {
        c = (a + b) / 2.0;
        double fc = f(c);

        if (fa * fc < 0) {
            b = c;
            fb = fc;
        } else {
            a = c;
            fa = fc;
        }
    }

    x = (a + b) / 2.0;
    return 0;
}

Newton_class::Newton_class() {}

Newton_class::~Newton_class() {}

void Newton_class::setVolumes(double vol_a, double vol_b) {
    a = vol_a;
    b = vol_b;
}

void Newton_class::setTolerance(double vol_eps) {
    eps = vol_eps;
}

double Newton_class::f(double x) {
    return x + sqrt(x) + cbrt(x) - 2.5;
}

double Newton_class::df(double x) {
    return 1.0 + 1.0/(2.0*sqrt(x)) + 1.0/(3.0*pow(cbrt(x), 2));
}

int Newton_class::count_newton(double &x) {
    x = (a + b) / 2.0;
    int max_iter = 100;

    for (int i = 0; i < max_iter; ++i) {
        double fx = f(x);
        double dfx = df(x);

        if (fabs(dfx) < 1e-12) return 2;

        double x_new = x - fx / dfx;

        if (x_new < a || x_new > b) {
            x_new = (a + b) / 2.0;
        }

        if (fabs(x_new - x) < eps) {
            x = x_new;
            return 0;
        }
        x = x_new;
    }

    return 1;
}