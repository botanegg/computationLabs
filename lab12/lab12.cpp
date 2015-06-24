#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>

#include "share/matrix.h"

using namespace std;

double u(double x, double y, double alpha) {
    return ((x + alpha) * (x + alpha)) / (1 + y);
}

double f(double x, double y, double alpha) {
    return (2 / (1 + y)) * (1 + ((x + alpha) / (1 + y)) * ((x + alpha) / (1 + y)));
}

void simpleIterMethod(double h, double alpha, double eps, double x1, double x0, double y1, double y0) {
    size_t n = (size_t) ((x1 - x0) / h + 1);
    Matrix solutionMatrix = Matrix::getE(n, n);

    //initial rules
    for (size_t i = 0; i < n; i++) {
        solutionMatrix[0][i] = u(x0 + i * h, y0, alpha);
        solutionMatrix[n - 1][i] = u(x0 + i * h, y1, alpha);
    }

    for (size_t i = 1; i < n - 1; i++) {
        solutionMatrix[i][0] = u(x0, y0 + h * i, alpha);
        for (size_t j = 1; j < n; j++)
            solutionMatrix[i][j] = 0;
        solutionMatrix[i][n - 1] = u(x1, y0 + h * i, alpha);
    }


    Matrix tmpMatrix;
    double EPS = 1000000;
    int iter = 1;
    while (EPS > eps) {
        iter++;
        tmpMatrix = solutionMatrix;
		for (size_t i = 1; i < n - 1; ++i)
		    for (size_t j = 1; j < n - 1; ++j)
                tmpMatrix[i][j] = (solutionMatrix[i - 1][j] + solutionMatrix[i + 1][j] + solutionMatrix[i][j - 1] +
                                   solutionMatrix[i][j + 1]) / 4 - (h * h * f(x0 + h * j, y0 + h * i, alpha)) / 4;
        EPS = 0;
		for (size_t i = 0; i < n; ++i)
		    for (size_t j = 0; j < n; ++j)
                EPS = max(EPS, abs(tmpMatrix[i][j] - solutionMatrix[i][j]));

        solutionMatrix = tmpMatrix;
    }

    cout << "iterations = " << iter << ";  Epsilon = " << EPS << endl;
    EPS = 0;
	for (size_t i = 0; i < n; ++i)
	    for (size_t j = 0; j < n; ++j)
            EPS = max(EPS, abs(solutionMatrix[i][j] - u(x0 + j * h, y0 + i * h, alpha)));
    cout << EPS << endl << endl;
}

void GaussZelde(double h, double alpha, double eps, double x1, double x0, double y1, double y0) {
    size_t n = (size_t) ((x1 - x0) / h + 1);
    Matrix solutionMatrix = Matrix::getE(n, n);
	for (size_t i = 0; i < n; i++) {
        solutionMatrix[0][i] = u(x0 + i * h, y0, alpha);
        solutionMatrix[n - 1][i] = u(x0 + i * h, y1, alpha);
    }
	for (size_t i = 1; i < n - 1; i++) {
        solutionMatrix[i][0] = u(x0, y0 + h * i, alpha);
		for (size_t j = 1; j < n; j++)
            solutionMatrix[i][j] = 0;
        solutionMatrix[i][n - 1] = u(x1, y0 + h * i, alpha);
    }

    Matrix tmpMatrix;
    double EPS = 1000000;
    int iter = 1;
    while (EPS > eps) {
        iter++;

        tmpMatrix = solutionMatrix;
		for (size_t i = 1; i < n - 1; ++i)
		    for (size_t j = 1; j < n - 1; ++j)
                solutionMatrix[i][j] = (solutionMatrix[i - 1][j] + solutionMatrix[i + 1][j] + solutionMatrix[i][j - 1] +
                                        solutionMatrix[i][j + 1]) / 4 - (h * h * f(x0 + h * j, y0 + h * i, alpha)) / 4;
        EPS = 0;
		for (size_t i = 0; i < n; ++i)
		    for (size_t j = 0; j < n; ++j)
                EPS = max(EPS, abs(tmpMatrix[i][j] - solutionMatrix[i][j]));

    }
    cout << "iterations = " << iter << ";  Epsilon = " << EPS << endl;
    EPS = 0;
	for (size_t i = 0; i < n; ++i)
	    for (size_t j = 0; j < n; ++j)
            EPS = max(EPS, abs(solutionMatrix[i][j] - u(x0 + j * h, y0 + i * h, alpha)));
    cout << EPS << endl << endl;
}


int main() {
    setlocale(LC_ALL, "");
    double x0 = 0;
    double x1 = 1;
    double y0 = 0;
    double y1 = 1;
    double alpha = 1;
    double h = 0.1;
    double eps = 0.0001;
    cout << "Simple Iteration method: " << endl;
    simpleIterMethod(h, alpha, eps, x1, x0, y1, y0);
    cout << "Gauss-Zelde method: " << endl;
    GaussZelde(h, alpha, eps, x1, x0, y1, y0);

    return 0;
}
