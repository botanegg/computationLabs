#ifndef SHARE_DERIVABLE_H
#define SHARE_DERIVABLE_H

class Derivable {
    double val, deriv;

public:
    Derivable(double _val, double _deriv);

    Derivable(double c);

    static Derivable IndependendVariable(double x);

    double getValue() const;

    double getDerivative() const;

    friend Derivable operator+(const Derivable &f1, const Derivable &f2);

    friend Derivable operator-(const Derivable &f1, const Derivable &f2);

    friend Derivable operator*(const Derivable &f1, const Derivable &f2);

    friend Derivable operator/(const Derivable &f1, const Derivable &f2);

    friend Derivable cos(Derivable f);

    friend Derivable sin(Derivable f);
};

#endif //SHARE_DERIVABLE_H