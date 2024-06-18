/*
#include <iostream>
#include <vector>
#include <cmath>

const double EPSILON = 1e-10;
using namespace std;

void printMatrix(const vector<vector<double>>& A) {
    for (const auto& row : A) {
        for (double elem : row) {
            cout << elem << "\t";
        }
        cout << endl;
    }
    cout << endl;
}

bool gaussElimination(vector<vector<double>>& A, vector<double>& x) {
    int n = A.size();
    for (int i = 0; i < n; i++) {
        int maxRow = i;
        for (int j = i + 1; j < n; j++) {
            if (abs(A[j][i]) > abs(A[maxRow][i])) {
                maxRow = j;
            }
        }
        if (abs(A[maxRow][i]) < EPSILON) {
            return false; // Матрица вырожденная или бесконечное число решений
        }
        swap(A[i], A[maxRow]);
        for (int j = i + 1; j < n; j++) {
            double ratio = A[j][i] / A[i][i];
            for (int k = i; k <= n; k++) {
                A[j][k] -= ratio * A[i][k];
            }
        }
    }

    for (int i = n - 1; i >= 0; i--) {
        x[i] = A[i][n];
        for (int j = i + 1; j < n; j++) {
            x[i] -= A[i][j] * x[j];
        }
        x[i] /= A[i][i];
    }

    return true;
}

int main() {
    vector<vector<double>> A = {
        { 1,  2, -1, -7, -23},
        { 8,  0, -9, -3,  39},
        { 2, -3,  7,  1,  -7},
        { 1, -5, -6,  8,  30}
    };

    const int n = A.size();
    vector<double> x(n, 0);

    if (gaussElimination(A, x)) {
        cout << "Решение системы:\n";
        for (int i = 0; i < n; i++) {
            cout << "x[" << i << "] = " << x[i] << endl;
        }
    } else {
        cout << "Система уравнений вырожденная или имеет бесконечное количество решений.\n";
    }

    return 0;
}
*/

#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

// Функция для вычисления определителя матрицы 3x3
double determinant3x3(vector<vector<double>>& matrix) {
    return matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1])
         - matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0])
         + matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]);
}

// Метод Крамера для решения системы линейных уравнений
vector<double> solveCramer(vector<vector<double>>& A, vector<double>& b) {
    int n = A.size();
    vector<double> x(n, 0);

    double detA = determinant3x3(A);

    if (detA == 0.0) {
        cout << "Система уравнений вырожденная (определитель матрицы равен нулю)." << endl;
        return x;
    }

    vector<vector<double>> tempA = A;

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            tempA[j][i] = b[j];
        }
        x[i] = determinant3x3(tempA) / detA;
        tempA = A;
    }

    return x;
}

// Функция для проверки диагонального доминирования
bool isDiagonallyDominant(vector<vector<double>>& A) {
    int n = A.size();
    for (int i = 0; i < n; i++) {
        double diag = abs(A[i][i]);
        double sum = 0.0;
        for (int j = 0; j < n; j++) {
            if (i != j) {
                sum += abs(A[i][j]);
            }
        }
        if (diag <= sum) {
            return false;
        }
    }
    return true;
}

// Метод Якоби для решения системы линейных уравнений
vector<double> solveJacobi(vector<vector<double>>& A, vector<double>& b, int maxIter = 100, double tolerance = 1e-6) {
    int n = A.size();
    vector<double> x(n, 0);
    vector<double> xNew(n, 0);
    int iter = 0;
    double error = tolerance + 1.0;

    if (!isDiagonallyDominant(A)) {
        cout << "Матрица не является диагонально доминирующей." << endl;
        return x;
    }

    while (error > tolerance && iter < maxIter) {
        for (int i = 0; i < n; i++) {
            double sum = 0.0;
            for (int j = 0; j < n; j++) {
                if (j != i) {
                    sum += A[i][j] * x[j];
                }
            }
            xNew[i] = (b[i] - sum) / A[i][i];
        }

        error = 0.0;
        for (int i = 0; i < n; i++) {
            error += abs(xNew[i] - x[i]);
            x[i] = xNew[i];
        }

        iter++;
    }

    if (iter == maxIter) {
        cout << "Достигнуто максимальное число итераций." << endl;
    } else {
        cout << "Сходимость достигнута на итерации " << iter << endl;
    }

    return x;
}

int main() {
    vector<vector<double>> A = {
        { 2, -2,  5 },
        { -2,  3,  6 },
        { -10, 12, -4 }
    };

    vector<double> b = { 2, -1, 2 };

    // Выбор метода решения
    char method;
    cout << "Выберите метод решения (C - Крамер, J - Якоби): ";
    cin >> method;

    vector<double> x;
    if (method == 'C' || method == 'c') {
        x = solveCramer(A, b);
        cout << "Решение системы методом Крамера:" << endl;
    } else if (method == 'J' || method == 'j') {
        x = solveJacobi(A, b);
        cout << "Решение системы методом Якоби:" << endl;
    } else {
        cout << "Некорректный выбор метода." << endl;
        return 1;
    }

    for (int i = 0; i < x.size(); i++) {
        cout << "x[" << i << "] = " << x[i] << endl;
    }

    return 0;
}
