#include <iostream>
#include <iomanip>
#include "Dyhotomia_class.h"

int main() {
    Dyhotomia_class solver;

    solver.setVolumes(0.4, 1.0);
    solver.setTolerance(1e-6);

    double rootDichotomy = solver.solveByDichotomy();
    double rootNewton = solver.solveByNewton();

    std::cout << "Дихотомія: " << std::setprecision(5) << rootDichotomy << std::endl;
    std::cout << "Ньютона: " << std::setprecision(5) << rootNewton << std::endl;

    return 0;
}
