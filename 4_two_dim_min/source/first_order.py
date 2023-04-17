import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d as ax, Axes3D
import functions as f


x_true = [-1.635575672731708, -0.24490976121576505]
f_true = f.f(x_true[0], x_true[1])
m = 2
M = 51
coef = m * (1 + m / M)
alpha = (3 - 5 ** 0.5) / 2


# def condition(a, b, f_1, f_2, grad, eps):
#     norm_grad_sq = (grad[0]**2 + grad[1]**2)
#     if  <= eps:
#     # if (abs(b - a) < eps):
#         return True
#     else:
#         return False


def golden_section_method(func, eps: float, grad: list, x: list):
    n = 0
    a = 0
    b = 1
    alpha_1 = -10e10
    alpha_2 = -10e10
    func_alpha_1 = -10e10
    func_alpha_2 = 10e10
    alpha_2_flag = 1  # флаг, что alpha_2 нужно пересчитать
    alpha_1_flag = 1  # флаг, что alpha_1 нужно пересчитать
    N = 100
    while N != 0:
        if abs(b - a) < eps:
            return (a + b) / 2, n
        if alpha_1_flag == 1:
            alpha_1 = a + alpha * (b - a)
            func_alpha_1 = func(x[0] - alpha_1 * grad[0], x[1] - alpha_1 * grad[1])
            n += 1

        if alpha_2_flag == 1:
            alpha_2 = b - alpha * (b - a)
            func_alpha_2 = func(x[0] - alpha_2 * grad[0], x[1] - alpha_2 * grad[1])
            n += 1

        if func_alpha_1 > func_alpha_2:
            a = alpha_1
            alpha_1 = alpha_2
            func_alpha_1 = func_alpha_2
            alpha_1_flag = 0
            alpha_2_flag = 1

        else:
            b = alpha_2
            alpha_2 = alpha_1
            func_alpha_2 = func_alpha_1
            alpha_1_flag = 1
            alpha_2_flag = 0

        N = N - 1


def dichotomy_method(func, eps: float, grad: list, x: list):
    n = 0
    a = 0
    b = 1
    N = 100
    while N != 0:
        if abs(b - a) < eps:
            return (a + b) / 2, n
        delta = (b - a) * 0.001
        alpha_1 = (a + b) / 2 - delta
        alpha_2 = (a + b) / 2 + delta
        func_alpha_1 = func(x[0] - alpha_1 * grad[0], x[1] - alpha_1 * grad[1])
        func_alpha_2 = func(x[0] - alpha_2 * grad[0], x[1] - alpha_2 * grad[1])
        n += 2
        if func_alpha_1 > func_alpha_2:
            a = alpha_1
        if func_alpha_1 <= func_alpha_2:
            b = alpha_2

        N = N - 1


def method_of_steepest_descent(func, grad_func, eps, method=golden_section_method):
    x_k = -1
    y_k = -1
    num_of_func_call = 0
    numb_of_iter = 0
    norms_grad = []
    alpha = 0  # скорость сходимости
    grad = grad_func(x_k, y_k)
    while True:
        norm_grad = f.norm(grad[0], grad[1])
        # norms_grad.append(norm_grad)
        if norm_grad < eps:
            return x_k, y_k, num_of_func_call, numb_of_iter, norms_grad, alpha
        result_find_alpha_k = method(func, eps, grad, [x_k, y_k])
        alpha_k = result_find_alpha_k[0]
        num_of_func_call += result_find_alpha_k[1]
        x_k -= alpha_k * grad[0]
        y_k -= alpha_k * grad[1]

        # x_k_1 = x_k - alpha_k * grad[0]
        # y_k_1 = y_k - alpha_k * grad[1]
        # alpha_iter = f.norm(x_k_1 - x_true[0], y_k_1 - x_true[1]) / f.norm(x_k - x_true[0], y_k - x_true[1])
        # if alpha_iter > alpha:
        #     alpha = alpha_iter

        # q = 1 - eps * alpha_k * coef
        # func_k = func(x_k, y_k)
        # x_k = x_k_1
        # y_k = y_k_1
        # func_k_plus_1 = func(x_k, y_k)
        # print(f"{func_k_plus_1 - f_true} <= {q * (func_k - f_true)}, {(func_k_plus_1 - f_true) <= q * (func_k - f_true)}")
        # print(q)
        numb_of_iter += 1
        grad = grad_func(x_k, y_k)


def compare_dich_golden_section():
    eps_set = [1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9]
    n_gold_sect, n_dich = [], []
    for eps in eps_set:
        n_gold_sect.append(method_of_steepest_descent(f.f, f.grad_f, eps, golden_section_method)[2])
        n_dich.append(method_of_steepest_descent(f.f, f.grad_f, eps, dichotomy_method)[2])

    plt.figure()
    # plt.semilogx(eps_set, n_gold_sect, linewidth=0.5, label="golden section method")
    plt.semilogx(eps_set, n_dich, linewidth=0.5, label="dichotomy method")
    plt.legend()
    plt.xlabel("eps")
    plt.ylabel("number of call function")
    plt.title("Зависимость числа вызовов функции от точности")
    plt.show()


def plot_and_solve():
    # x, y = np.meshgrid(np.linspace(-2, -0.5, 1000), np.linspace(-5, 0.5, 1000))
    x, y = np.meshgrid(np.linspace(-2, 2, 1000), np.linspace(-2, 2, 1000))
    result = method_of_steepest_descent(f.f, f.grad_f, 10e-5)

    print(f'first order: [{result[0]}, {result[1]}]')

    z = f.f(x, y)
    plt.figure()
    plt.contour(x, y, z, 20)
    # plt.scatter(-1, -1, color='red')
    plt.scatter(result[0], result[1], color='blue')
    # plt.scatter(-1, -f.SQRT_5, color='red')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.show()
    fig = plt.figure(figsize=(7, 4))
    ax_3d = fig.add_subplot(projection='3d')
    ax_3d.plot_wireframe(x, y, z)
    plt.show()


def norm_grad():
    result = method_of_steepest_descent(f.f, f.grad_f, 10e-5, golden_section_method)

    numb_of_iter = np.linspace(0, result[3], result[3] + 1)
    print(f'result: [{result[0]},{result[1]}]')
    print(f'alpha = {result[5]}')

    plt.figure()
    plt.semilogy(numb_of_iter, result[4])
    plt.xlabel("number of iteration")
    plt.ylabel("gradient norm")
    plt.title("Норма градиента от номера итерации")
    plt.show()
