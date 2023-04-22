//
// Created by Тимур Жексимбаев on 22.04.2023.
//
#include<iomanip>
#include<iostream>

using namespace std;

class Matrix {
public:
    int n, m;
    int d = 1.0;
    double matrix[100][100];

    Matrix(int row, int column) {
        this->n = row;
        this->m = column;
    }

    Matrix operator+(Matrix v) {
        Matrix sumMatrix(n, m);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                sumMatrix.matrix[i][j] = matrix[i][j] + v.matrix[i][j];
            }
        }
        return sumMatrix;
    }

    Matrix operator-(Matrix v) {
        Matrix subMatrix(n, m);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                subMatrix.matrix[i][j] = matrix[i][j] - v.matrix[i][j];
            }
        }
        return subMatrix;
    }

    Matrix operator*(Matrix v) {
        Matrix mulMatrix(n, v.m);
        for (int k = 0; k < n; k++) {
            for (int i = 0; i < v.m; i++) {
                double sum = 0.0;
                for (int j = 0; j < m; j++) {
                    sum += matrix[k][j] * v.matrix[j][i];
                }
                mulMatrix.matrix[k][i] = sum;
            }
        }
        return mulMatrix;
    }

    friend istream &operator>>(istream &input, Matrix &v) {
        for (int i = 0; i < v.n; i++) {
            for (int j = 0; j < v.m; j++) {
                input >> v.matrix[i][j];
            }
        }
        return input;
    }

    friend ostream &operator<<(ostream &output, Matrix &v) {
        for (int i = 0; i < v.n; i++) {
            for (int j = 0; j < v.m; j++) {
                output << fixed << setprecision(4) << v.matrix[i][j] << ' ';
            }
            output << '\n';
        }
        return output;
    }

    Matrix &operator=(Matrix v) {

        Matrix assignMatrix(n, m);
        for (int i = 0; i < v.n; i++) {
            for (int j = 0; j < v.m; j++) {
                assignMatrix.matrix[i][j] = v.matrix[i][j];
            }
        }
        return assignMatrix;
    }

    Matrix Transpose() {

        Matrix transposeMatrix(m, n);

        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                transposeMatrix.matrix[i][j] = matrix[j][i];
            }
        }
        return transposeMatrix;
    }

};


int main() {
    /// Input data, for example
    // 2 3
    // 1 2 3
    // 6 3 8

    int n, m;
    cin >> n >> m;
    Matrix A(n, m);
    cin >> A;


    return 0;
}
