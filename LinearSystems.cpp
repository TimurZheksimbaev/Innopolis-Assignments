
//Timur Zheksimbaev
// t.zheksimbaev@innopolis.university
#include<iostream>
#include<vector>
#include<iomanip>

using namespace std;

// Class Matrix from previous Assigment 2.
class Matrix {

public:
    double matrix[100][100];
    int n, m;

    Matrix(int row, int column) {
        this->n = row;
        this->m = column;
    }

    void print() {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n - 1; j++)
                cout << fixed << setprecision(2) << matrix[i][j] << " ";

            cout << fixed << setprecision(2) << matrix[i][n - 1] << '\n';
        }
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
    }

    // for printing matrices
    friend ostream &operator<<(ostream &output, Matrix &v) {

        for (int i = 0; i < v.n; i++) {
            for (int j = 0; j < v.m; j++) {
                output << fixed << setprecision(2) << v.matrix[i][j] << "\n";
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


// I created this class for storing square matrices, because in task there are square matrices.
// And Gauss-Jordan cannot be applied to non-square matrices.
// It inherits class Matrix
class SquareMatrix : public Matrix {
public:

    SquareMatrix(int n) : Matrix(n, n) {}

    /// (((((((((((((((((((((((((((((((((((((((( LINEAR SYSTEMS ))))))))))))))))))))))))))))))))))))))))))))))))))))

    //function for computing pivots
    int pivot(int ii, int jj) {
        int pivotIndex = ii;
        double mx = abs(matrix[ii][jj]);

        for (int i = ii + 1; i < n; i++)
            if (abs(matrix[i][jj]) > mx)
                mx = abs(matrix[i][jj]), pivotIndex = i;
        return pivotIndex;
    }

    // function for solving system of linear equations
    void linearSystems(Matrix b) {
        bool NO = false, INF = false;

        //swapping and upper triangle
        for (int k = 0; k < n; k++) {

            //finding pivot for swapping rows
            int p = pivot(k, k);

            //checking and swapping rows
            if (p != k) {
                for (int j = 0; j < n; j++) {
                    swap(matrix[k][j], matrix[p][j]); // in matrix A
                    swap(b.matrix[k][j], b.matrix[p][j]); //in matrix b
                }
            }

            //reducing to upper triangle
            // [ * * * ]
            // [ 0 * * ]
            // [ 0 0 * ]
            for (int i = k + 1; i < n; i++) {
                double multiplier = matrix[i][k] / matrix[k][k];

                for (int j = 0; j < n; j++) {
                    matrix[i][j] -= multiplier * matrix[k][j];
                    b.matrix[i][j] -= multiplier * b.matrix[k][j];
                }
            }
        }

        //checking for no solutions and infinite solutions
        for (int i = 0; i < n; i++) {
            int zeros = 0; //counter for zeros
            for (int j = 0; j < n; j++) {
                if (matrix[i][j] == 0) {
                    zeros++;
                }
            }
            // if a whole row consists of zeros, then
            if (zeros == n) {
                //if component b corresponding to that row is zero too, then infinite solutions
                if (b.matrix[i][0] == 0) INF = true;
                    //if component b corresponding to that row is not zero then no solutions
                else NO = true;
            }
        }


        // reducing to lower triangle
        // [ * 0 0 ]
        // [ * * 0 ]
        // [ * * * ]
        for (int k = n - 1; k >= 0; k--) {
            for (int i = k - 1; i >= 0; i--) {
                double multiplier = matrix[i][k] / matrix[k][k];
                for (int j = 0; j < n; j++) {
                    matrix[i][j] -= multiplier * matrix[k][j];
                    b.matrix[i][j] -= multiplier * b.matrix[k][j];
                }
            }
        }

        //diagonal normalization and zeros
        // [ 1 0 0 ]
        // [ 0 1 0 ]
        // [ 0 0 1 ]
        for (int i = 0; i < n; i++) {
            b.matrix[i][0] /= matrix[i][i];
            if (abs(b.matrix[i][0]) < 0.01) b.matrix[i][0] = 0.0; // if numbers are very small
        }

        //outputting results
        if (NO) cout << "NO";
        else if (INF) cout << "INF";
        else cout << b;

    }
};

// Class for storing column vectors, inherits class Matrix
class ColumnVector : public Matrix {
public:
    ColumnVector(int n) : Matrix(n, 1) {}
};

int main() {
    int n; //dimensions of matrix
    cin >> n;

    SquareMatrix A(n); //matrix A
    cin >> A;

    ColumnVector c(n); //matrix b
    cin >> c;

    A.linearSystems(c); //calling function, see Class SquareMatrix above ^

    return 0;
}




