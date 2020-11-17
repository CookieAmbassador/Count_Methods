import numpy as np


def krylov(A):
    print("Метод Крилова")
    print("")
    # print("Input size of the matrix ")
    # n = int(input())
    n = len(A)
    # print("Input matrix's elements ")
    # A = np.zeros((n, n))
    # for i in range(n):
    #   for j in range(n):
    #      A[i][j] = float(input())
    print("")
    Y = np.zeros((n + 1, n))
    Y[0][0] = 1
    while True:
        for i in range(1, len(Y)):
            Y[i] = A.dot(Y[i - 1])
        Y_n = -Y[n]
        Y_new = np.rot90(Y[:n], 3)
        if np.linalg.det(Y_new) != 0:
            break
        Y[0][0] += 1
    P = np.linalg.solve(Y_new, Y_n)
    P = np.insert(P, 0, 1)
    Lambdas = np.roots(P)
    print("Власні значення")
    print(Lambdas)
    Q = np.zeros((n, n))
    Q[0] = 1
    for j in range(1, n):
        Q[j] = Lambdas * Q[j - 1] + P[j]
    X = np.zeros((n, n))
    for i in range(0, n):
        X[i] = Y[n - 1]
        for j in range(1, n):
            X[i] += Q[j][i] * Y[n - j - 1]

    for i in range(0, n):
        X[i] = X[i] / np.max(np.abs(X[i]))
    print("")
    print("Власні вектори")
    w, v = np.linalg.eig(A)
    print(X)
    print("")
    print("Нев'язка знайденних векторів")
    for i in range(0, n):
        print(A.dot(X[i]) - Lambdas[i] * X[i])
        print("")
    print("")
    # print("Нев'язка векторів, знайдених стандартним пакетом")
    # for i in range(0, n):
    #   print(np.linalg.matrix_power(A, i + 1).dot(w[i]) - Lambdas[i] * w[i])
    #  print("")
    print("____")


def simple_vector_iteration(A):
    print("Метод простої векторної ітерації")
    print("")
    n = len(A)
    current_lambda = 0
    previous_lambda = 0
    iterations_counter = 0
    eps = 0.001
    Y0 = np.zeros((n,))
    Y0[0] = 1
    Y0[1] = 1
    Y0[2] = 2
    while True:
        previous_y = Y0
        Y0 = A.dot(Y0)
        current_lambda = np.sum(Y0 / previous_y) / n
        iterations_counter += 1
        if abs(current_lambda - previous_lambda) < eps:
            break
        previous_lambda = current_lambda
    print("Ітерацій знадобилось ")
    print(iterations_counter)
    print("")
    print("Найбільше власне значення")
    print(current_lambda)
    print("")
    Y = Y0 / np.max(np.abs(Y0))
    print("Відповідний йому власний вектор")
    print(Y)
    print("")
    print("Вектор нев'язки ")
    print(A.dot(Y) - current_lambda * Y)
    print("____")


def spectral_number(A):
    A = A.T.dot(A)
    w, v = np.linalg.eig(A)
    spectral_num = np.sqrt(np.max(w) / np.min(w))
    print("")
    print("Спектральне число обумовленості")
    print(spectral_num)


A = np.array([[2.8, 1, 1.8], [1, 3.3, 1.8], [1.8, 1.8, 3.8]])
print("Початкова матриця:")
print("")
print(A)
print("")
w, v = np.linalg.eig(A)
print("Еталонні результати:")
print("Власні значення")
print(w)
print("")
print("Власні вектори")
print(v)
print("")
print("Нев'язка")
for i in range(0, len(v)):
    print(A.dot(v[i]) - w[i] * v[i])
    print("")
print("")

krylov(A)

simple_vector_iteration(A)

spectral_number(A)
