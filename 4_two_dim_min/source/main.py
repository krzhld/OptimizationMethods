import math
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d as ax, Axes3D


def f_1(x, y):
    return np.power((x - 2), 2) + np.power((y + 3), 2)


def f_2(x, y):
    return x + y + 4 * np.sqrt(1 + 2 * np.power(x,2) + 3 * np.power(y,2))


def f(x, y):
    return np.power(x,2) + 5 * np.power(y,2) + np.sin(4 * x + 5 * y) + 3 * x + 2 * y


def grad_f(x, y):
    return 2 * x + 4 * math.cos(4 * x + 5 * y) + 3, 10 * y + 5 * math.cos(4 * x + 5 * y) + 2


def grad_f_1(x, y):
    return 2 * (x - 2), 2 * (y + 3)


def grad_f_2(x, y):
    return 1 + (8 * x) / math.sqrt(1 + 2 * x**2 + 3 * y**2) , 1 + (12 * y) / math.sqrt(1 + 2 * x**2 + 3 * y**2) 


alpha = (3 - 5 ** 0.5) / 2


#модификая метода дихотомии для поиска alpha_k в методе наискорейшего спуска
def golden_section_method(func, eps : float, grad : list, x : list):
    n = 0
    a = 0
    b = 1
    alpha_1 = -10e10
    alpha_2 = -10e10
    func_alpha_1 = -10e10
    func_alpha_2 = -10e10
    alpha_2_flag = 1 #флаг, что alpha_2 нужно пересчитать
    alpha_1_flag = 1 #флаг, что alpha_1 нужно пересчитать
    while True:
        if abs(b - a) < eps:
            return (a + b) / 2, n
        if (alpha_1_flag == 1):
            alpha_1 = a + alpha * (b - a)
            func_alpha_1 = func(x[0] - alpha_1 * grad[0], x[1] - alpha_1 * grad[1])
            n += 1

        if (alpha_2_flag == 1):
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



#модификая метода дихотомии для поиска alpha_k в методе наискорейшего спуска
def dichotomy_method(func, eps : float, grad : list, x : list):
    n = 0
    a = 0 
    b = 1
    while True:
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

# def method_of_steepest_descent_dich(func, grad_func, eps):
#     x_k = -5
#     y_k = -2
#     num_of_func_call = 0
#     grad = grad_func(x_k, y_k)
#     while(True):
#         if((grad[0]**2 + grad[1]**2)**0.5 < eps):
#             return x_k, y_k, num_of_func_call
#         result_find_alpha_k = dichotomy_method(func, eps / 1000, grad, [x_k, y_k])
#         alpha_k = result_find_alpha_k[0]
#         num_of_func_call += result_find_alpha_k[1]
#         x_k -= alpha_k * grad[0]
#         y_k -= alpha_k * grad[1]
#         grad = grad_func(x_k, y_k)

def method_of_steepest_descent_gold_sect(func, grad_func, eps):
    x_k = -5
    y_k = -2
    num_of_func_call = 0
    numb_of_iter = 0 
    norms_grad = []
    grad = grad_func(x_k, y_k)
    while(True):
        norm_grad = (grad[0]**2 + grad[1]**2)**0.5
        norms_grad.append(norm_grad)
        if(norm_grad < eps):
            return x_k, y_k, num_of_func_call, numb_of_iter, norms_grad
        result_find_alpha_k = golden_section_method(func, eps / 1000, grad, [x_k, y_k])
        alpha_k = result_find_alpha_k[0]
        num_of_func_call += result_find_alpha_k[1]
        x_k -= alpha_k * grad[0]
        y_k -= alpha_k * grad[1]
        grad = grad_func(x_k, y_k)
        numb_of_iter += 1




# eps_set = [10e-2, 10e-3, 10e-4, 10e-5, 10e-6, 10e-7, 10e-8]

# n_gold_sect, n_dich = [], []

# for eps in eps_set:
#     n_gold_sect.append(method_of_steepest_descent_gold_sect(f, grad_f, eps)[2])
#     n_dich.append(method_of_steepest_descent_dich(f, grad_f, eps)[2])

# plt.figure()
# plt.semilogx(eps_set, n_gold_sect,  linewidth=0.5, label = "golden section method")
# plt.semilogx(eps_set, n_dich, linewidth=0.5, label = "dichotomy method")
# plt.legend()
# plt.xlabel("eps")
# plt.ylabel("number of call function")
# plt.title("Зависимость числа вызовов функции от точности")
# plt.show()


x, y = np.meshgrid(np.linspace(-5, 1, 1000), np.linspace(-2, 1.5, 1000))

result = method_of_steepest_descent_gold_sect(f, grad_f, 10e-7)

print(f'[{result[0]}, {result[1]}]')

z = f(x, y)
plt.figure()
plt.contour(x, y, z, 20)
plt.scatter(result[0], result[1])
plt.show()

numb_of_iter = np.linspace(0, result[3], result[3] + 1)

plt.figure()
plt.semilogy(numb_of_iter, result[4])
plt.xlabel("number of iteration")
plt.ylabel("gradient norm")
plt.title("Норма градиента от номера итерации")
plt.show()

# fig = plt.figure(figsize=(7, 4))
# ax_3d = fig.add_subplot(projection='3d')
# ax_3d.plot_wireframe(x, y, z)
# plt.show()
