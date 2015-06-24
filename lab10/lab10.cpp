// chislolab10.cpp: определяет точку входа для консольного приложения.
//

#include <iostream>
#include <iomanip>
#include <cmath>

#include "share/matrix.h"

using namespace std;

double u(double x, double t) {
    return (t * t + t + 1) * (x * x * x + x + 1);
}

void progonka(double *a, double *b, double *c, double *F, int n, double *res) {
    double m;
    for (int i = 0; i < n; i++) {
        m = a[i] / c[i - 1];
        c[i] = c[i] - m * b[i - 1];
        F[i] = F[i] - m * F[i - 1];
    }
    res[n - 1] = F[n - 1] / c[n - 1];
    for (int i = n - 2; i >= 0; i--)
        res[i] = (F[i] - b[i] * res[i + 1]) / c[i];
}

double f(double x, double t) ///when k=1
{
    return (-1) * 6 * t * t * x + 2 * t * (x * x * x - 2 * x + 1) - x * x * x - 5 * x - 1;
}

void clearSystem(double A, double B, double T, double T0, double h, double tau) {
    size_t m = (size_t) ((B - A) / h + 1);
    size_t n = (size_t) ((T - T0) / h + 1);

    Matrix Matr = Matrix::getE(n,m);

    for (int i = 0; i < m; i++)
        Matr[0][i] = u(A + h * i, 0);
    for (int i = 1; i < n; i++) {
        Matr[i][0] = u(A, tau * i);
        for (int j = 1; j < m - 1; j++)
            Matr[i][j] = Matr[i - 1][j] +
                         ((Matr[i - 1][j - 1] - 2 * Matr[i - 1][j] + Matr[i - 1][j + 1] + Matr[i - 1][j + 1]) /
                          (h * h) + f(A + j * h, i * tau)) * tau;
        Matr[i][m - 1] = u(B, tau * i);
    }
    cout << setw(10) << "U:" << " " << setw(10) << "solution:" << endl;
    for (int i = 0; i < m - 1; i++)
        cout << setw(8) << u(A + i * h, T) << " " << setw(10) << Matr[n - 1][i] << endl;
}

void unclearSystem(double A, double B, double T, double T0, double h, double tau) {
    size_t m = (size_t) ((B - A) / h + 1);
    size_t n = (size_t) ((T - T0) / h + 1);

    Matrix Matr = Matrix::getE(n,m);

    for (int i = 0; i < m; i++)
        Matr[0][i] = u(A + h * i, 0);


    for (int i = 1; i < n; i++) {
        Matr[i][0] = u(A, tau * i);
        Matr[i][m - 1] = u(B, tau * i);
    }
    for (int i = 1; i < n; i++) {
        double *a = new double[m - 2];
        double *b = new double[m - 2];
        double *c = new double[m - 2];
        double *F = new double[m - 2];
        double *res = new double[n - 2];
        b[0] = -tau / (h * h);
        c[0] = 1 + (2 * tau) / (h * h);
        F[0] = Matr[i - 1][1] + (tau * Matr[i][0]) / (h * h) + f(A + h, i * tau) * tau;

        for (int j = 1; j < m - 3; ++j) {
            a[j] = -tau / (h * h);
            c[j] = 1 + (2 * tau) / (h * h);
            b[j] = -tau / (h * h);
            F[j] = Matr[i - 1][j + 1] + f(A + (j + 1) * h, i * tau) * tau;
        }

        a[m - 3] = -tau / (h * h);
        c[m - 3] = 1 + (2 * tau) / (h * h);
        F[m - 3] = Matr[i - 1][m - 2] + (tau * Matr[i][m - 1]) / (h * h) + f(B - h, i * tau) * tau;
        progonka(a, b, c, F, m - 2, res);
        for (int j = 0; j < m - 2; j++)
            Matr[i][j + 1] = res[j];
        delete[]a;
        delete[]b;
        delete[]c;
        delete[]F;
        delete[]res;
    }
    cout << setw(10) << "U:" << " " << setw(10) << "solution" << endl;
    for (int i = 0; i < m; i++)
        cout << setw(10) << u(A + i * h, T) << " " << setw(10) << Matr[n - 1][i] << endl;
}

int main() {
    setlocale(LC_ALL, "");
    double h = 0.1;
    double A = 0, B = 1;
    double T0 = 0, T = 1;
    double tau = 0.0001;
    cout << "Clear system: " << endl;
    clearSystem(A, B, T, T0, h, tau);
    h = 0.1;
    tau = 0.1;
    cout << "Unclear system: " << endl;
    h = 0.1;
    tau = 0.1;
    unclearSystem(A, B, T, T0, h, tau);

    return 0;
}
