import math
import numpy as np
import matplotlib.pyplot as plt
import functions as f
import first_order as first


alpha = (3 - 5 ** 0.5) / 2


def mult_matrix_and_vector(A, b) -> list:
    return A[0][0] * b[0] + A[0][1] * b[1], A[1][0] * b[0] + A[1][1] * b[1]


def mult_vectors(v1: list, v2: list) -> float:
    return v1[0] * v2[0] + v1[1] * v2[1]


def newton_method(func, grad_func, hess_func, eps):
    x_k = -1
    y_k = -1
    alpha_0 = 1
    delta = 0.5 
    lambd = 0.5
    numb_of_iter = 0

    while True:
        grad = grad_func(x_k, y_k)
        norm_grad = f.norm(grad[0], grad[1])
        if norm_grad < eps:
            return x_k, y_k, numb_of_iter 
        
        print(f'{norm_grad}')
        d_k = np.array(mult_matrix_and_vector(np.linalg.inv(hess_func(x_k, y_k)), grad))
        alpha_k = alpha_0

        F = func(x_k - alpha_k * d_k[0], y_k - alpha_k * d_k[1])
        f_k = func(x_k, y_k)
        while F - f_k > -delta * alpha_k * mult_vectors(grad, d_k):
            alpha_k *= lambd

            F = func(x_k - alpha_k * d_k[0], y_k - alpha_k * d_k[1])
            
        x_k -= alpha_k * d_k[0] 
        y_k -= alpha_k * d_k[1]
        numb_of_iter += 1


def plot_and_solve():
    x, y = np.meshgrid(np.linspace(-4, 1, 1000), np.linspace(-1.5, 1, 1000))
    # x, y = np.meshgrid(np.linspace(-5, 5, 1000), np.linspace(-5, 5, 1000))

    z = f.f(x, y)
    plt.figure()
    plt.contour(x, y, z, 50)
    result = newton_method(f.f, f.grad_f, f.hess_f, 10e-4)
    plt.scatter(result[0], result[1])
    plt.show()

    print(f'Newton: [{result[0]}, {result[1]}]')

    result1 = first.method_of_steepest_descent(f.f, f.grad_f, 10e-4)
    print(f'Gradient: [{result1[0]}, {result1[1]}]')


def compare_grad_and_newton():
    epsilons = [1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8]
    result_1, result_2 = [], []

    for eps in epsilons:
        result_1.append(first.method_of_steepest_descent(f.f, f.grad_f, eps)[2])
        result_2.append(newton_method(f.f, f.grad_f, f.hess_f, eps)[2])

    plt.figure()
    plt.semilogx(epsilons, result_1, linewidth=2, label="method of steepest descent")
    plt.semilogx(epsilons, result_2, linewidth=2, label="newton method")
    plt.legend()
    plt.grid()
    plt.xlabel("eps")
    plt.ylabel("number of iterations")
    plt.title("Зависимость числа итераций от точности")
    plt.show()
