#include <iostream>
#include <cmath>
#include <vector>

#include "share/vector.h"
#include "share/matrix.h"
#include "share/utils.h"


double c1 = 1;
double c2 = 2;
double d1 = 1;
double d2 = 2;

double c = 3;
double d = 3;

double a = 2;
double b = 3;

size_t n = 10;

using namespace std;

double p(double x) {
    return x / (x + 1);
}

double q(double x) {
    return sin(x);
}

double f(double x) {
    return exp(-x);
}

void progon() {


    double h = (b - a) / n;

    Matrix matrix;
    matrix.setDimension(n + 1, n + 1);

    matrix[0][0] = c1 - c2 / h;
    matrix[0][1] = c2 / h;

    for (size_t i = 1; i < n; i++) {
        double xi = a + h * i;
        matrix[i][i - 1] = 1 / (h * h) - p(xi) / (2 * h);
        matrix[i][i] = q(xi) - 2 / (h * h);
        matrix[i][i + 1] = 1 / (h * h) + p(xi) / (2 * h);
    }

    matrix[n][n - 1] = -d2 / h;
    matrix[n][n] = d1 + d2 / h;

    Vector vector;
    vector.setDimension(n + 1);

    vector[0] = c;

    for (size_t i = 1; i < n; i++) {
        double xi = a + h * i;
        vector[i] = f(xi);
    }

    vector[n] = d;

    Vector ans = Utils::solveProgon(matrix, vector);
    Utils::printVector(ans);
   // Utils::printMatrix(matrix);

}

double f1(double x, double y1, double y2) {
    return y2;
}

double f2(double x, double y1, double y2) {
    return f(x) - p(x) * y2 - q(x) * y1;
}

vector<pair<double, double>> rungeKuttaSystemSolution(double (*fxu1)(double x, double y1, double y2),
                                                      double (*fxu2)(double x, double y1, double y2), double a,
                                                      double b, double n, double eta, double psi) {
    vector<pair<double, double>> res;
    double h = (b - a) / n;
    res.reserve(n + 1);

    res.push_back(make_pair(eta, psi)); // first elem;

    for (int i = 0; i < n; i++) {
        double yn = res.back().first;
        double yn2 = res.back().second;
        double xn = a + i * h;

        double k11 = fxu1(xn, yn, yn2);
        double k12 = fxu2(xn, yn, yn2);

        double k21 = fxu1(xn + h / 2, yn + (h * k11) / 2, yn2 + (h * k12) / 2);
        double k22 = fxu2(xn + h / 2, yn + (h * k11) / 2, yn2 + (h * k12) / 2);

        double k31 = fxu1(xn + h / 2, yn + (h * k21) / 2, yn2 + (h * k22) / 2);
        double k32 = fxu2(xn + h / 2, yn + (h * k21) / 2, yn2 + (h * k22) / 2);

        double k41 = fxu1(xn + h, yn + (h * k31), yn2 + (h * k32));
        double k42 = fxu2(xn + h, yn + (h * k31), yn2 + (h * k32));

        res.push_back(
                make_pair(yn + (k11 + 2 * k21 + 2 * k31 + k41) * h / 6, yn2 + (k12 + 2 * k22 + 2 * k32 + k42) * h / 6));


//        double k1 = fxu(xn, yn);
//        double k2 = fxu(xn + h / 4, yn + (h * k1) / 4);
//        double k3 = fxu(xn + h / 2, yn + (h * k2) / 2);
//        double k4 = fxu(xn + h, yn + (h * (k1 - 2*k2 + 2*k3)));
//
//        res.push_back(yn + (k1  + 4 * k3 + k4) * h / 6);
    }
    return res;
}

int main() {
    progon();

    double eta1 = 1;
    double beta1 = (c - c1 * eta1) / c2;

    auto res1 = rungeKuttaSystemSolution(f1, f2, a, b, n, eta1, beta1);
    double psi1 = d1 * res1[res1.size() - 1].first + d2 * res1[res1.size() - 1].second - d;

    double eta2 = 2;
    double beta2 = (c - c1 * eta2) / c2;

    auto res2 = rungeKuttaSystemSolution(f1, f2, a, b, n, eta2, beta2);
    double psi2 = d1 * res2[res2.size() - 1].first + d2 * res2[res2.size() - 1].second - d;

    double eta = eta2 - ((eta2 - eta1) / (psi2 - psi1)) * psi2;
    double beta = (c - c1 * eta) / c2;

    auto res = rungeKuttaSystemSolution(f1, f2, a, b, n, eta, beta);

    for (auto &el:res) {
        cout << el.first << ' ';
    }
    cout << endl;

    return 0;
}
