#include <iostream>
#include <cmath>
#include <vector>
#include <float.h>


using namespace std;

//General
double f(double x);
double a = 0;
double b = 3.1415926;

//2
double formulaTrapeze(double a, double b, int m);
double formulaSimpson(double a, double b, int m);
void ruleRunge();

//3
double formulaGauss(double a, double b, int m);


int main() {
    ruleRunge();
    cout << "Gauss m = 2: " << formulaGauss(a, b, 2) << endl;
    cout << "Gauss m = 4: " << formulaGauss(a, b, 4) << endl;
    cout << "Gauss m = 5: " << formulaGauss(a, b, 5) << endl;
    cout << "Gauss m = 8: " << formulaGauss(a, b, 8) << endl;

    return 0;
}


double f(double x) {
    return exp(cos(x)) * cos(2 * x);
}

void ruleRunge() {
    cout << "I. Runge's rule" << endl;

    double epsMinusFive = 0.00001;
    double epsMinusSeven = 0.0000001;


    double eps = DBL_MAX;
    int m = 1;
    double zhDivTwo, zh;


    //e10-5
    while (eps > epsMinusFive) {
        zh = formulaTrapeze(a, b, m);
        zhDivTwo = formulaTrapeze(a, b, m * 2);
        eps = abs((zhDivTwo - zh) / 3);
        ++m;
    }

    cout << "formulaTrapeze" << endl;
    cout << "Res: " << zhDivTwo << endl;
    cout << "m: " << m * 2 << endl;
    cout << "h: " << (b - a) / (m * 2) << endl;
    cout << "Eps: " << eps << endl << endl;


    //e10-7
    eps = DBL_MAX;
    m = 1;

    while (eps > epsMinusSeven) {
        zh = formulaTrapeze(a, b, m);
        zhDivTwo = formulaTrapeze(a, b, m * 2);
        eps = abs((zhDivTwo - zh) / 3);
        ++m;
    }

    cout << "formulaTrapeze" << endl;
    cout << "Res: " << zhDivTwo << endl;
    cout << "m: " << m * 2 << endl;
    cout << "h: " << (b - a) / (m * 2) << endl;
    cout << "Eps: " << eps << endl << endl;

    //Simpson
    eps = DBL_MAX;
    m = 1;

    while (eps > epsMinusSeven) {
        zh = formulaSimpson(a, b, m);
        zhDivTwo = formulaSimpson(a, b, m * 2);
        eps = abs((zhDivTwo - zh) / 15);
        ++m;
    }

    cout << "formulaSimpson" << endl;
    cout << "Res: " << zhDivTwo << endl;
    cout << "m: " << m * 4 << endl;
    cout << "h: " << (b - a) / (m * 2) << endl;
    cout << "Eps: " << eps << endl << endl;
}


double formulaTrapeze(double a, double b, int m) {
    double res = 0;
    double h = (b - a) / m;

    for (int i = 0; i <= m; ++i)
        res += f(a + i * h);
    res = (2 * res - f(a) - f(a + h * m)) * (h / 2);


    return res;
}

double formulaSimpson(double a, double b, int m) {
    double res = 0;
    double h = (b - a) / (2 * m);

    res += f(a) + f(a + h * 2 * m);

    for (int i = 1; i < 2 * m; ++i) {
        if (i % 2) res += 4 * f(a + i * h); else
            res += 2 * f(a + i * h);
    }
    res *= (b - a) / (6 * m);
    return res;
}

double formulaGauss(double a, double b, int m) {
    vector<vector<double>> t, ag;
    t.resize(9);
    ag.resize(9);


    //2
    t[2].resize(2);
    t[2][0] = -0.577350;
    t[2][1] = 0.577350;

    ag[2].resize(2);
    ag[2][0] = ag[2][1] = 1;

    //4
    t[4].resize(4);
    t[4][0] = -0.86114;
    t[4][1] = -0.33998;
    t[4][2] = 0.33998;
    t[4][3] = 0.86114;

    ag[4].resize(4);
    ag[4][0] = ag[4][3] = 0.34785;
    ag[4][1] = ag[4][2] = 0.65215;

    //5
    t[5].resize(5);
    t[5][0] = -0.90618;
    t[5][1] = -0.538469;
    t[5][2] = 0;
    t[5][3] = 0.538469;
    t[5][4] = 0.90618;

    ag[5].resize(5);
    ag[5][0] = ag[5][4] = 0.23693;
    ag[5][1] = ag[5][3] = 0.47863;
    ag[5][2] = 0.56889;


    //8
    t[8].resize(8);
    t[8][0] = -0.96028986;
    t[8][1] = -0.79666648;
    t[8][2] = -0.52553242;
    t[8][3] = -0.18343464;

    t[8][7] = 0.96028986;
    t[8][6] = 0.79666648;
    t[8][5] = 0.52553242;
    t[8][4] = 0.18343464;

    ag[8].resize(8);
    ag[8][0] = ag[8][7] = 0.10122854;
    ag[8][1] = ag[8][6] = 0.22238103;
    ag[8][2] = ag[8][5] = 0.31370664;
    ag[8][3] = ag[8][4] = 0.36268378;


    double res = 0;
    for (int i = 0; i < m; ++i) {
        res += ag[m][i] * f((a + b) / 2 + (b - a) * t[m][i] / 2);

    }
    res *= (b - a) / 2;

    return res;

}