import numpy as np

'''
def check_convergence(A):
    n = len(A)
    for i in range(n):
        row_sum = sum(abs(A[i, j]) for j in range(n) if i != j)
        if abs(A[i, i]) <= row_sum:
            return False
    return True

def seidel(A, b, tol=0.01, max_iter=1000):
    n = len(A)
    x = np.zeros(n)
    x_new = np.zeros(n)
    iter_count = 0
    error = tol + 1

    if not check_convergence(A):
        print("Метод Зейделя не сходится для данной матрицы A.")
        return None

    while error > tol and iter_count < max_iter:
        for i in range(n):
            x_new[i] = b[i]
            for j in range(n):
                if j != i:
                    x_new[i] -= A[i, j] * x[j]
            x_new[i] /= A[i, i]

        error = np.max(np.abs(x_new - x))
        x = np.copy(x_new)
        iter_count += 1

    if iter_count == max_iter:
        print("Достигнуто максимальное число итераций.")
    else:
        print(f"Сходимость достигнута на итерации {iter_count}")

    return x

A = np.array([
    [23, -6, -5, 9],
    [8, 22, -2, 5],
    [7, -6, 18, -1],
    [3, 5, 5, -19]
], dtype=float)

b = np.array([232, -82, 202, -57], dtype=float)

solution = seidel(A, b)

if solution is not None:
    print("Решение системы методом Зейделя:")
    for i, sol in enumerate(solution):
        print(f"x[{i}] = {sol:.3f}")
'''
def thomas(A, b):
    n = len(b)
    c = np.zeros(n-1)
    d = np.zeros(n)
    x = np.zeros(n)

    c[0] = A[0, 1] / A[0, 0]
    d[0] = b[0] / A[0, 0]

    for i in range(1, n-1):
        c[i] = A[i, i+1] / (A[i, i] - A[i, i-1] * c[i-1])
    for i in range(1, n):
        d[i] = (b[i] - A[i, i-1] * d[i-1]) / (A[i, i] - A[i, i-1] * c[i-1])

    x[n-1] = d[n-1]
    for i in range(n-2, -1, -1):
        x[i] = d[i] - c[i] * x[i+1]

    return x

A = np.array([
    [6, -5, 0, 0, 0],
    [-6, 16, 9, 0, 0],
    [0, 9, -17, -3, 0],
    [0, 0, 8, 22, -8],
    [0, 0, 0, 6, -13]
], dtype=float)

b = np.array([-58, 161, -114, -90, -55], dtype=float)

solution = thomas(A, b)

print("Решение системы методом прогонки:")
for i, sol in enumerate(solution):
    print(f"x[{i}] = {sol:.3f}")
