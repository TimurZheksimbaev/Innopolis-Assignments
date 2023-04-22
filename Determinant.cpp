//Timur Zheksimbaev
//t.zheksimbaev@innopolis.university
#include<iostream>
#include<iomanip>

using namespace std;


class Matrix {

public:
    double matrix[100][100];
    int n, m;

    Matrix(int row, int column) {
        this->n = row;
        this->m = column;
    }

    // for adding matrices
    Matrix operator+(Matrix v) {

        Matrix sumMatrix(n, m);

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                sumMatrix.matrix[i][j] = matrix[i][j] + v.matrix[i][j];
            }
        }
        return sumMatrix;
    }

    // for subtracting matrices
    Matrix operator-(Matrix v) {

        Matrix subMatrix(n, m);

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                subMatrix.matrix[i][j] = matrix[i][j] - v.matrix[i][j];
            }
        }
        return subMatrix;
    }

    // for multiplying matrices
    Matrix operator*(Matrix v) {

        Matrix mulMatrix(n, v.m);

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < v.m; j++) {
                for (int k = 0; k < m; k++) {
                    mulMatrix.matrix[i][j] += matrix[i][k] * v.matrix[k][j];
                }
            }
        }
        return mulMatrix;
    }

    // for transposing matrices
    Matrix Transpose() {

        Matrix transposeMatrix(m, n);

        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                transposeMatrix.matrix[i][j] = matrix[j][i];
            }
        }
        return transposeMatrix;
    }

    // for assigning matrices
    Matrix &operator=(Matrix v) {

        Matrix assignMatrix(n, m);
        for (int i = 0; i < v.n; i++) {
            for (int j = 0; j < v.m; j++) {
                matrix[i][j] = v.matrix[i][j];
            }
        }
        return assignMatrix;
    }

    // for printing matrices
    friend ostream &operator<<(ostream &output, Matrix &v) {

        for (int i = 0; i < v.n; i++) {
            for (int j = 0; j < v.m; j++) {
                if (j == v.m-1) cout << v.matrix[i][j];
                else output << fixed << setprecision(2) << v.matrix[i][j] << "\n";
            }
        }
        return output;
    }

    // for filling matrices
    friend istream &operator>>(istream &input, Matrix &v) {
        for (int i = 0; i < v.n; i++) {
            for (int j = 0; j < v.m; j++) {
                input >> v.matrix[i][j];
            }
        }
        return input;
    }
};


// I created this class because in the task there are square matrices,
// and we cannot find determinant of non-square matrix
// It inherits class Matrix
class SquareMatrix : public Matrix {
public:

    SquareMatrix(int n) : Matrix(n, n) {}

    //function for printing matrix after each operation
    void print() {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (j == n-1) cout << matrix[i][j];
                else cout << fixed << setprecision(2) << matrix[i][j] << " ";
            }
            cout << '\n';
        }
    }

    //function for computing pivots by maximum
    int pivot(int ii, int jj) {
        int pivotIndex = ii;
        double mx = abs(matrix[ii][jj]);

        for (int i = ii + 1; i < n; i++)
            if (abs(matrix[i][jj]) > mx)
                mx = abs(matrix[i][jj]), pivotIndex = i;
        return pivotIndex;
    }

    // function for calculating determinant
    void Determinant() {

        double det = 1.0; //determinant of a matrix
        int numberOfPermutations = 0; // if we exchange rows in a matrix, then determinant changes its sign

        // swapping rows and upper triangle form
        for (int k = 0; k < n-1; k++) {

            int p = pivot(k, k); //finding pivot

            //swapping
            if (p != k) {
                numberOfPermutations++;
                cout << "step: permutation\n";
                for (int j = 0; j < n; j++) {
                    swap(matrix[k][j], matrix[p][j]);
                }
                print();
            }

            // reducing to upper triangle form
            for (int i = k + 1; i < n; i++) {
                cout << "step: elimination\n";

                double multiplier = matrix[i][k] / matrix[k][k];
                for (int j = 0; j < n; j++) {
                    matrix[i][j] -= multiplier * matrix[k][j];
                }
                print();
            }
        }

        //calculating determinant
        for (int i = 0; i < n; i++) {
            det *= matrix[i][i];
        }

        // when we exchange rows determinant changes its sign, so if we exchange rows odd number of times,
        // then we need to multiply by -1.0
        if (numberOfPermutations % 2 != 0) det *= -1.0;

        cout << "result:\n" << det;
    }
};


int main() {
    int n;
    cout << "Input the dimensions of a matrix:\nn: ";
    cin >> n;

    SquareMatrix A(n);
    cout << "Input n*n elements of a matrix:\n";
    cin >> A;

    A.Determinant();
    return 0;
}