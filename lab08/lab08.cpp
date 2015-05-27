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
size_t n = 20;
double a = -1.5;
double b = 0.5;

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

vector<double> rungeKuttaSolution(double (*fxu)(double x, double y), double a, double b, size_t n) {
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

vector<double> adamsSolution(double (*fxu)(double x, double y), double a, double b, size_t n, double y0,
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

vector<double> realSolution(double (*fx)(double x), double a, double b, size_t n) {
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


vector<pair<double, double>> rungeKuttaSystemSolution(double (*fxu1)(double x, double y1, double y2),
                                                      double (*fxu2)(double x, double y1, double y2), double a,
                                                      double b, size_t n) {
    vector<pair<double, double>> res;
    double h = (b - a) / n;
    res.reserve(n + 1);

    res.push_back(make_pair(fxu1(a, 1, 1), fxu2(a, 1, 1))); // first elem;

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


vector<pair<double, double>> adamsSystemSolution(double (*fxu1)(double x, double y1, double y2),
                                                 double (*fxu2)(double x, double y1, double y2), double a, double b,
                                                 size_t n, vector<pair<double, double>> kutta) {
    vector<pair<double, double>> res;
    double h = (b - a) / n;
    res.reserve(n + 1);

    res.push_back(kutta[0]); // first elem;
    res.push_back(kutta[1]);
    res.push_back(kutta[2]);
    res.push_back(kutta[3]);

    for (int i = 3; i < n; i++) {
        double yn = (res.end() - 1)->first;
        double yn_1 = (res.end() - 2)->first;
        double yn_2 = (res.end() - 3)->first;
        double yn_3 = (res.end() - 4)->first;

        double _2yn = (res.end() - 1)->second;
        double _2yn_1 = (res.end() - 2)->second;
        double _2yn_2 = (res.end() - 3)->second;
        double _2yn_3 = (res.end() - 4)->second;

        double xn = a + i * h;
        double xn_1 = a + (i - 1) * h;
        double xn_2 = a + (i - 2) * h;
        double xn_3 = a + (i - 3) * h;

        res.push_back(
                make_pair(yn + h * (55 * fxu1(xn, yn, _2yn) - 59 * fxu1(xn_1, yn_1, _2yn_1) +
                                    37 * fxu1(xn_2, yn_2, _2yn_2) - 9 * fxu1(xn_3, yn_3, _2yn_3)) / 24,
                          _2yn + h * (55 * fxu2(xn, yn, _2yn) - 59 * fxu2(xn_1, yn_1, _2yn_1) +
                                      37 * fxu2(xn_2, yn_2, _2yn_2) - 9 * fxu2(xn_3, yn_3, _2yn_3)) / 24));


//        double k1 = fxu(xn, yn);
//        double k2 = fxu(xn + h / 4, yn + (h * k1) / 4);
//        double k3 = fxu(xn + h / 2, yn + (h * k2) / 2);
//        double k4 = fxu(xn + h, yn + (h * (k1 - 2*k2 + 2*k3)));
//
//        res.push_back(yn + (k1  + 4 * k3 + k4) * h / 6);
    }
    return res;
}

double test_fxu1(double x, double u, double u2) {
    if (x == 0) return 0.6;
    return x * cos(u + u2);
//    if (x == 0) return 0;
//    return (x * x + 1) * u + cos(x) + (-1 - x * x) * sin(x);
}

double test_fxu2(double x, double u, double u2) {
    if (x == 0) return 2;
    return sin(u - u2);
//    if (x == 0) return 0;
//    return (x * x + 1) * u + cos(x) + (-1 - x * x) * sin(x);
}

int main() {
    auto runge_res = runge();
    auto adams_res = adams(runge_res);
    auto real_res = real();

    cout.precision(8);
    cout.setf(std::ios::fixed, std::ios::floatfield);
    for (int i = 0; i <= n; i++) {
        cout << a + i * (b - a) / n << " " << real_res[i] << " " << runge_res[i] << " " << adams_res[i] <<
        " " <<
        abs(runge_res[i] - real_res[i]) << " " << abs(adams_res[i] - real_res[i]) << endl;
    }
    cout << endl;
    cout << "runge eps " << runge_eps << endl;
    cout << "adams eps " << adams_eps << endl;

    auto res = rungeKuttaSystemSolution(test_fxu1, test_fxu2, 0, 2, 20);
    for (int i = 0; i <= 20; i++) {
        cout << res[i].first << endl;

    }
    cout << endl << endl;
    for (int i = 0; i <= 20; i++) {
        cout << res[i].second << endl;

    }

    cout << endl << endl;

    auto res2 = adamsSystemSolution(test_fxu1, test_fxu2, 0, 2, 20, res);
    for (int i = 0; i <= 20; i++) {
        cout << res2[i].first << endl;

    }
    cout << endl << endl;
    for (int i = 0; i <= 20; i++) {
        cout << res2[i].second << endl;

    }
    return 0;
}
