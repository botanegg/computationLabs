#include <iostream>
#include <fstream>
#include <memory.h>

void deleteArr(double **arr, int size);

void writeVec(double *vec, int size);

void writeMatrix(double **U, int size);

using namespace std;

int main() {
    double **matrix = NULL;
    double **L = NULL;
    double **U = NULL;
    fstream file("input.txt");
    int size = 0;
    file >> size;

    int size_d1 = size - 1;
    int size_i1 = size + 1;


    matrix = new double *[size];
    L = new double *[size];
    U = new double *[size];

    for (int i = 0; i < size; ++i) {
        matrix[i] = new double[size_i1];

        L[i] = new double[size];
        U[i] = new double[size];
        memset(L[i], 0, sizeof(double) * size);
        memset(U[i], 0, sizeof(double) * size);
    }

    //read matrix[][] + orig_matrix fill
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            file >> matrix[i][j];
        }
        file >> matrix[i][size];
    }

    file.close();


    /* COMPUTE */
    //compute LU method
    //forward

    //compute L U matrix
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            double sum = 0;
            for (int k = 0; k < i; ++k) {
                sum += L[i][k] * U[k][j];
            }

            if (i <= j) {
                //compute U matrix
                U[i][j] = matrix[i][j] - sum;
            } else {
                //compute L matrix
                L[i][j] = (matrix[i][j] - sum) / U[j][j];
            }

            if (i == j) L[i][j] = 1;
        }
    }

    //Ly = b
    //Ux = y

    //compute y[]
    //forward L

    double *y = NULL;
    y = new double[size];

    for (int i = 0; i < size; ++i) {
        double sum = 0;
        for (int j = 0; j < i; ++j)
            sum += L[i][j] * y[j];
        double L_ii_mul_x = matrix[i][size] - sum;
        y[i] = L_ii_mul_x;
    }

    //compute x[]
    //backward U

    double *x = NULL;
    x = new double[size];

    for (int i = size_d1; i >= 0; --i) {
        double sum = 0;
        for (int j = i + 1; j < size; ++j)
            sum += U[i][j] * x[j];
        double U_ii_mul_x = y[i] - sum;
        x[i] = U_ii_mul_x / U[i][i];
    }

    //compute residual

    double *residual = NULL;
    residual = new double[size];

    for (int i = 0; i < size; ++i) {
        double sum = 0;
        for (int j = 0; j < size; ++j) {
            sum += matrix[i][j] * x[j];
        }
        residual[i] = sum - matrix[i][size];
    }


    /* WRITE */
    //write matrix[][]
    cout << "write matrix[][]" << endl;
    writeMatrix(matrix, size);

    //write L[][]
    cout << "write L[][]" << endl;
    writeMatrix(L, size);

    //write U[][]
    cout << "write U[][]" << endl;
    writeMatrix(U, size);

    //write y[]
    cout << "write y[]" << endl;
    writeVec(y, size);

    //write x[]
    cout << "write x[]" << endl;
    writeVec(x, size);

    //write residual[]
    cout << "write residual[]" << endl;
    writeVec(residual, size);


    /* CLEAN */
    //clean matrix array
    deleteArr(matrix, size);

    //clean U array
    deleteArr(U, size);

    //clean L array
    deleteArr(L, size);

    delete[] matrix;
    delete[] U;
    delete[] L;
    delete[] y;
    delete[] x;
    delete[] residual;

    return 0;
}

void writeMatrix(double **U, int size) {
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            cout << U[i][j] << ' ';
        }
        cout << endl;
    }
    cout << endl;
}

void writeVec(double *vec, int size) {
    for (int i = 0; i < size; ++i) {
        cout << vec[i] << endl;
    }
    cout << endl;
}

void deleteArr(double **arr, int size) {
    for (int i = 0; i < size; ++i) {
        delete[] arr[i];
    }
}
