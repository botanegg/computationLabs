#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

int main() {
    double **matrix = NULL;
    double **orig_matrix = NULL;
    fstream file("input.txt");
    int size = 0;
    file >> size;

    int size_d1 = size - 1;
    int size_i1 = size + 1;


    matrix = new double *[size];
    orig_matrix = new double *[size];

    for (int i = 0; i < size; ++i) {
        matrix[i] = new double[size_i1];
        orig_matrix[i] = new double[size_i1];
    }

    //read matrix[][] + orig_matrix fill
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            file >> matrix[i][j];
            orig_matrix[i][j] = matrix[i][j];
        }
        file >> matrix[i][size];
        orig_matrix[i][size] = matrix[i][size];
    }

    file.close();


    /* COMPUTE */
    //compute matrix[][] LU method
    //forward

    //!TODO write it

    //compute x[]
    //backward

    double *x = NULL;
    x = new double[size];

    //!TODO write it

    //compute residual

    double *residual = NULL;
    residual = new double[size];

    for (int i = 0; i < size; ++i) {
        double sum = 0;
        for (int j = 0; j < size; ++j) {
            sum += orig_matrix[i][j] * x[j];
        }
        residual[i] = sum - orig_matrix[i][size];
    }


    /* WRITE */
    //write matrix[][]
    cout << "write matrix[][]" << endl;
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            cout << matrix[i][j] << ' ';
        }
        cout << matrix[i][size] << endl;
    }
    cout << endl;

    //write x[]
    cout << "write x[]" << endl;
    for (int i = 0; i < size; ++i) {
        cout << x[i] << endl;
    }
    cout << endl;

    //write residual[]
    cout << "write residual[]" << endl;
    for (int i = 0; i < size; ++i) {
        cout << residual[i] << endl;
    }
    cout << endl;


    /* CLEAN */
    //clean matrix array
    for (int i = 0; i < size; ++i) {
        delete[] matrix[i];
    }

    //clean orig_matrix array
    for (int i = 0; i < size; ++i) {
        delete[] orig_matrix[i];
    }

    delete[] matrix;
    delete[] orig_matrix;
    delete[] x;
    delete[] residual;

    return 0;
}