// VAR 6 u = sin(x^2)
// u' = (x^2 + 1)u + f(x)
// u(0) = 0;
// => f(x) = 2x*cos(x^2)+(-1-x^2)sin(x^2)
//
// u' = (x^2 + 1)u + 2x*cos(x^2)+(-1-x^2)sin(x^2)
// u(0) = 0;

#include <iostream>
#include <cmath>
#include <vector>

double runge_eps = 0;
double adams_eps = 0;
int n = 20;
double a = 0;
double b = 1;

using namespace std;

double fxu(double x, double u) {
    if (x == 0) return 0;
    return (x * x + 1) * u + 2 * x * cos(x * x) + (-1 - x * x) * sin(x * x);
}

double ans(double x) {
    return sin(x * x);
}

double test_fxu(double x, double u) {
    if (x == -1.5) return 2 * exp(x);
    return 2 * x * u + 2 * exp(x) * (1 - 2 * x);
//    if (x == 0) return 0;
//    return (x * x + 1) * u + cos(x) + (-1 - x * x) * sin(x);
}

double test_ans(double x) {
    return 2 * exp(x);
}

vector<double> rungeKuttaSolution(double (*fxu)(double x, double y), double a, double b, double n) {
    vector<double> res;
    double h = (b - a) / n;
    res.reserve(n + 1);

    res.push_back(fxu(a, 1)); // first elem;

    for (int i = 0; i < n; i++) {
        double yn = res.back();
        double xn = a + i * h;

        double k1 = fxu(xn, yn);
        double k2 = fxu(xn + h / 2, yn + (h * k1) / 2);
        double k3 = fxu(xn + h / 2, yn + (h * k2) / 2);
        double k4 = fxu(xn + h, yn + (h * k3));

        res.push_back(yn + (k1 + 2 * k2 + 2 * k3 + k4) * h / 6);

//        double k1 = fxu(xn, yn);
//        double k2 = fxu(xn + h / 4, yn + (h * k1) / 4);
//        double k3 = fxu(xn + h / 2, yn + (h * k2) / 2);
//        double k4 = fxu(xn + h, yn + (h * (k1 - 2*k2 + 2*k3)));
//
//        res.push_back(yn + (k1  + 4 * k3 + k4) * h / 6);
    }
    return res;
}

vector<double> adamsSolution(double (*fxu)(double x, double y), double a, double b, double n, double y0,
                             double y1, double y2, double y3) {
    vector<double> res;
    double h = (b - a) / n;
    res.reserve(n + 1);

    res.push_back(y0);
    res.push_back(y1);
    res.push_back(y2);
    res.push_back(y3);


    for (int i = 3; i < n; i++) {
        double yn = *(res.end() - 1);
        double yn_1 = *(res.end() - 2);
        double yn_2 = *(res.end() - 3);
        double yn_3 = *(res.end() - 4);

        double xn = a + i * h;
        double xn_1 = a + (i - 1) * h;
        double xn_2 = a + (i - 2) * h;
        double xn_3 = a + (i - 3) * h;

        res.push_back(
                yn + h * (55 * fxu(xn, yn) - 59 * fxu(xn_1, yn_1) + 37 * fxu(xn_2, yn_2) - 9 * fxu(xn_3, yn_3)) / 24);

//        double k1 = fxu(xn, yn);
//        double k2 = fxu(xn + h / 4, yn + (h * k1) / 4);
//        double k3 = fxu(xn + h / 2, yn + (h * k2) / 2);
//        double k4 = fxu(xn + h, yn + (h * (k1 - 2*k2 + 2*k3)));
//
//        res.push_back(yn + (k1  + 4 * k3 + k4) * h / 6);
    }
    return res;
}

vector<double> realSolution(double (*fx)(double x), double a, double b, double n) {
    vector<double> res;
    double h = (b - a) / n;
    res.reserve(n + 1);

    for (int i = 0; i <= n; i++) {
        double xn = a + i * h;
        res.push_back(fx(xn));
    }
    return res;
}


vector<double> runge() {
    auto res1 = rungeKuttaSolution(test_fxu, a, b, n);
    auto res2 = rungeKuttaSolution(test_fxu, a, b, 2 * n);


    for (int i = 0; i <= n; i++) {
        runge_eps = max(runge_eps, abs(res1[i] - res2[2 * i]));
    }
    runge_eps /= 15;

//    for (auto item : res1) {
//        cout << item << endl;
//    }
//    cout << "eps: " << runge_eps << endl;

    return res1;
}

vector<double> adams(vector<double> runge) {
    auto res1 = adamsSolution(test_fxu, a, b, n, runge[0], runge[1], runge[2], runge[3]);
    auto res2 = adamsSolution(test_fxu, a, b, 2 * n, runge[0], runge[1], runge[2], runge[3]);
    for (int i = 0; i <= n; i++) {
        adams_eps = max(adams_eps, abs(res1[i] - res2[2 * i]));
    }
    adams_eps /= 15;

//    for (auto item : res1) {
//        cout << item << endl;
//    }
//    cout << "eps: " << adams_eps << endl;

    return res1;
}

vector<double> real() {
    auto res1 = realSolution(test_ans, a, b, n);
    return res1;
}

int main() {
    auto runge_res = runge();
    auto adams_res = adams(runge_res);
    auto real_res = real();

    cout.precision(12);
    cout.setf(std::ios::fixed, std::ios::floatfield);
    for (int i = 0; i <= n; i++) {
        cout << a + i * (b - a) / n << " \t" << real_res[i] << " \t" << runge_res[i] << " \t" << adams_res[i] <<
        " \t" <<
        abs(runge_res[i] - real_res[i]) << " \t" << abs(adams_res[i] - real_res[i]) << endl;
    }
    cout << endl;
    cout << "runge eps " << runge_eps << endl;
    cout << "adams eps " << adams_eps << endl;
    return 0;
}
