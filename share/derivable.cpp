#include "derivable.h"

#include <cmath>

Derivable operator+(const Derivable &f1, const Derivable &f2) {
    return Derivable(f1.val + f2.val, f1.deriv + f2.deriv);
}

Derivable operator-(const Derivable &f1, const Derivable &f2) {
    return Derivable(f1.val - f2.val, f1.deriv - f2.deriv);
}

Derivable operator*(const Derivable &f1, const Derivable &f2) {
    return Derivable(f1.val * f2.val, f1.deriv * f2.val + f1.val * f2.deriv);
}

Derivable operator/(const Derivable &f1, const Derivable &f2) {
    return Derivable(f1.val / f2.val, (f1.deriv * f2.val - f1.val * f2.deriv) / f2.val / f2.val);
}

Derivable cos(Derivable f) {
    return Derivable(cos(f.val), -sin(f.val) * f.deriv);
}

Derivable sin(Derivable f) {
    return Derivable(sin(f.val), cos(f.val) * f.deriv);
}

Derivable::Derivable(double _val, double _deriv) : val(_val), deriv(_deriv) {
}

Derivable::Derivable(double c) : val(c), deriv(0) {
}

Derivable Derivable::IndependendVariable(double x) {
    return Derivable(x, 1);
}

double Derivable::getValue() const {
    return val;
}

double Derivable::getDerivative() const {
    return deriv;
}