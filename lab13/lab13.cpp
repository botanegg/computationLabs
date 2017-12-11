#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <share/utils.h>

#include "share/matrix.h"


using namespace std;

double u(double x, double y, double alpha) {
	return ((x + alpha) * (x + alpha)) / (1 + y);
}

double F(double x) {
	return M_PI * M_PI * sin(M_PI * x) + sin(M_PI * x);
}


double getElemY(size_t n, size_t i, double x) {
	double delta = 1. / (n + 1);
	double pos = (i + 1) * delta;
	double dist = pos - x;

	if (abs(dist) >= delta) return 0;

	if (dist < 0) return 1 + dist / delta;
	else return 1 - dist / delta;
}

double getElemDiffY(size_t n, size_t i, double x) {
	double delta = 1. / (n + 1);
	double pos = (i + 1) * delta;
	double dist = pos - x;

	if (abs(dist) > delta) return 0;
	return dist >= 0 ? 1 / delta : -1 / delta;
}

double calcAij(size_t n, ssize_t i, ssize_t j) {
	double delta = 1. / (n + 1);
/*	ssize_t firstNumber = min(i, j);
	ssize_t secondNumber = max(i, j);

	double startPoint = (firstNumber + 1 - 1) * delta;
	double endPoint = (secondNumber + 1 + 1) * delta;

	double integral = 0;
	const double STEP = delta / 100;
	while (startPoint + STEP < endPoint) {
		double first = getElemY(n, firstNumber, startPoint);
		double second = getElemY(n, secondNumber, startPoint);

		double firstDiff = getElemDiffY(n, firstNumber, startPoint);
		double secondDiff = getElemDiffY(n, secondNumber, startPoint);

		integral += STEP * first * second;
		integral += STEP * firstDiff * secondDiff;

		startPoint += STEP;
	}
*/
	if (i == j) 
		integral = 2 * delta / 3 + 2 / delta;
	else 
		integral = delta / 6 - 1 / delta;
	
	return integral;
}


double calcBi(size_t n, ssize_t i) {
	double delta = 1. / (n + 1);


	double startPoint = ((i + 1) - 1) * delta;
	double endPoint = ((i + 1) + 1) * delta;

	double integral = 0;
	const double STEP = delta / 100;
	while (startPoint + STEP < endPoint) {
		double phi = getElemY(n, i, startPoint);
		double func = F(startPoint);

		integral += STEP * phi * func;

		startPoint += STEP;
	}

	return integral;
}


Matrix calculateMatrix(size_t n) {
	Matrix m = Matrix::getE(n, n);
	for (ssize_t i = 0; i < n; i++) {
		for (ssize_t j = 0; j < n; j++) {
			if (abs(i - j) >= 2)
				m[i][j] = 0;
			else
				m[i][j] = calcAij(n, i, j);
		}
	}

	Utils::printMatrix(m);
	cout << endl;

	return m;
}

Vector calculateVector(size_t n) {
	Vector v = Vector::get0(n);
	for (ssize_t i = 0; i < n; i++) {
		v[i] = calcBi(n, i);
	}

	Utils::printVector(v);
	cout << endl;

	return v;
}


int main() {
	setlocale(LC_ALL, "");

	size_t n = 100;
	auto mat = calculateMatrix(n);
	auto vec = calculateVector(n);

	auto solve = Utils::solveProgon(mat, vec);

	Utils::printVector(solve);
	cout << endl;
	return 0;
}
