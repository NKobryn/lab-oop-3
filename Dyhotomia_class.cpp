#include "Dyhotomia_class.h"
#include <cmath>
#include <stdexcept>

Dyhotomia_class::Dyhotomia_class() : a(0), b(0), eps(1e-6) {}
Dyhotomia_class::~Dyhotomia_class() {}

void Dyhotomia_class::setVolumes(double vol_a, double vol_b) {
    a = vol_a;
    b = vol_b;
}

void Dyhotomia_class::setTolerance(double vol_eps) {
    eps = vol_eps;
}

double Dyhotomia_class::func(double x) const {
    return x + sqrt(x) + cbrt(x) - 2.5;
}

double Dyhotomia_class::derivative(double x) const {
    return 1.0 + (1.0 / (2.0 * sqrt(x))) + (1.0 / (3.0 * cbrt(x * x)));
}

double Dyhotomia_class::solveByDichotomy() {
    double left = a;
    double right = b;

    while ((right - left) / 2 > eps) {
        double middle = (left + right) / 2.0;
        if (func(middle) == 0.0) {
            return middle;
        }
        if (func(left) * func(middle) < 0) {
            right = middle;
        } else {
            left = middle;
        }
    }

    return (left + right) / 2.0;
}

double Dyhotomia_class::solveByNewton() const {
    double x = (a + b) / 2;
    double delta = eps;

    while (fabs(func(x)) > eps) {
        double fx = func(x);
        double dfx = (func(x + delta) - func(x - delta)) / (2 * delta);

        if (fabs(dfx) < 1e-12) {
            std::cerr << "Похідна дуже мала. Метод Ньютона може не працювати коректно.\n";
            return x;
        }

        double x_new = x - fx / dfx;

        x = x_new;
    }

    return x;
}
