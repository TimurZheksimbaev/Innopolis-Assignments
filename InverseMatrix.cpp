//
// Created by Тимур Жексимбаев on 22.04.2023.
//

#include<iomanip>
#include<iostream>
#include<cfloat>
#include<cmath>

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

    double Norm() {
        double norm = 0.0;
        for (int i = 0; i < n; i++) {
            norm += pow(matrix[i][0], 2);
        }
        return sqrt(norm);
    }

};

class SquareMatrix : public Matrix {
public:

    SquareMatrix(int n) : Matrix(n, n) {}

    SquareMatrix &operator=(Matrix v) {

        SquareMatrix assignMatrix(n);
        for (int i = 0; i < v.n; i++) {
            for (int j = 0; j < v.n; j++) {
                matrix[i][j] = v.matrix[i][j];
            }
        }
        return assignMatrix;
    }

    /// (((((((((((((((((((((((((((((((((((((((((( INVERSE )))))))))))))))))))))))))))))))))))))))))))))))))))))))))

    SquareMatrix Inverse() {
        /// Filling Identity Matrix
        SquareMatrix id(n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (i == j) id.matrix[i][j] = 1.0;
                else id.matrix[i][j] = 0.0;
            }
        }

        ///  Pivoting and swapping rows
        for (int t = 0; t < n - 1; t++) {
            int I = t;
            double mx = matrix[t][t];
            for (int i = t + 1; i < n; i++) {
                if (abs(matrix[i][t]) > mx) {
                    mx = abs(matrix[i][t]);
                    I = i;
                }
            }
            if (I != t) {
                for (int j = 0; j < n; j++) {
                    swap(matrix[t][j], matrix[I][j]);
                    swap(id.matrix[t][j], id.matrix[I][j]);
                }
            }

            /// Upper Triangle Form
            for (int i = t + 1; i < n; i++) {
                double T = -matrix[i][t] / matrix[t][t];
                for (int j = 0; j < n; j++) {
                    matrix[i][j] += T * matrix[t][j];
                    id.matrix[i][j] += T * id.matrix[t][j];
                }
            }
        }

        /// Lower Triangle Form
        for (int t = n - 1; t >= 0; t--) {
            for (int i = t - 1; i >= 0; i--) {
                double T = -matrix[i][t] / matrix[t][t];
                for (int j = 0; j < n; j++) {
                    matrix[i][j] += T * matrix[t][j];
                    id.matrix[i][j] += T * id.matrix[t][j];
                }
            }
        }

        /// Diagonal Normalization
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                id.matrix[i][j] /= matrix[i][i];
            }
            matrix[i][i] = 1.0;
        }

        return id;
    }

    void printAug(Matrix identityMatrix) {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (abs(matrix[i][j]) < 0.01)matrix[i][j] = 0.0;
                cout << fixed << setprecision(2) << matrix[i][j] << " ";
            }

            for (int j = 0; j < n; j++) {
                if (abs(matrix[i][j]) < 0.01)matrix[i][j] = 0.0;
                cout << fixed << setprecision(2) << identityMatrix.matrix[i][j] << " ";
            }
            cout << '\n';

        }
    }

    int pivotAug(int pi, int pj) {
        int pivot = pi;
        double mx = matrix[pi][pj];
        for (int i = pi + 1; i < n; i++) {
            if (abs(matrix[i][pj]) > mx) {
                mx = abs(matrix[i][pj]);
                pivot = i;
            }
        }
        return pivot;
    }

    void Inverse(Matrix identityMatrix) { /// With printing steps

        int step = 1; // step of eliminations and permutations

        cout << "step #0: Augmented Matrix\n";
        printAug(identityMatrix);
        cout << "Direct way:\n";
        for (int k = 0; k < n - 1; k++) {
            int p = pivotAug(k, k);

            if (p != k) {
                cout << "step #" << step++ << ": permutation\n";
                for (int j = 0; j < n; j++) {
                    swap(matrix[k][j], matrix[p][j]);
                    swap(identityMatrix.matrix[k][j], identityMatrix.matrix[p][j]);
                }
                printAug(identityMatrix);
            }

            for (int i = k + 1; i < n; i++) {
                if (matrix[i][k] == 0.0) continue; /// If we don't need elimination.
                cout << "step #" << step++ << ": elimination\n";

                double multiplier = matrix[i][k] / matrix[k][k];
                for (int j = 0; j < n; j++) {
                    matrix[i][j] -= multiplier * matrix[k][j];
                    identityMatrix.matrix[i][j] -= multiplier * identityMatrix.matrix[k][j];
                }
                printAug(identityMatrix);
            }
        }

        cout << "Way back:\n";

        for (int t = n - 1; t >= 0; t--) {
            for (int i = t - 1; i >= 0; i--) {
                cout << "step #" << step++ << ": elimination\n";
                double multiplier = matrix[i][t] / matrix[t][t];
                for (int j = 0; j < n; j++) {
                    matrix[i][j] -= multiplier * matrix[t][j];
                    identityMatrix.matrix[i][j] -= multiplier * identityMatrix.matrix[t][j];
                }
                printAug(identityMatrix);
            }
        }

        cout << "Diagonal normalization:\n";
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                identityMatrix.matrix[i][j] /= matrix[i][i];
            }
            matrix[i][i] = 1.0;
        }
        printAug(identityMatrix);
        cout << "result:\n";
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                cout << fixed << setprecision(2) << identityMatrix.matrix[i][j] << " ";
            }
            cout << "\n";
        }
    }
};


int main() {

    return 0;
}
