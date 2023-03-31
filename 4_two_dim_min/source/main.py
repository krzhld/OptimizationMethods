import math
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d as ax, Axes3D


def f_1(x, y):
    return np.power((x - 2), 2) + np.power((y + 3), 2)

def f_2(x,y):
    return x + y + 4 * np.sqrt(1 + 2 * np.power(x,2) + 3 * np.power(y,2))

def f(x, y):
    return np.power(x,2) + 5 * np.power(y,2) + np.sin(4 * x + 5 * y) + 3 * x + 2 * y

def grad_f(x,y):
    return 2 * x + 4 * math.cos(4 * x + 5 * y) + 3, 10 * y + 5 * math.cos(4 * x + 5 * y) + 2

def grad_f_1(x, y):
    return 2 * (x - 2), 2 * (y + 3)

def grad_f_2(x, y):
    return 1 + (8 * x) / math.sqrt(1 + 2 * x**2 + 3 * y**2) , 1 + (12 * y) / math.sqrt(1 + 2 * x**2 + 3 * y**2) 

alpha = (3 - 5 ** 0.5) / 2

def golden_section_method(func, a, b, eps):
    n_eps = (math.log(eps) - math.log(b - a)) / math.log(1 - alpha)
    n_eps = math.ceil(n_eps)
    n_eps += 1
    n = 0
    x_1 = -10e10
    x_2 = -10e10
    func_x_1 = -10e10
    func_x_2 = -10e10
    x_2_flag = 1 #флаг, что x_2 нужно пересчитать
    x_1_flag = 1 #флаг, что x_1 нужно пересчитать
    print(f'predictive n: {n_eps}')
    while True:
        if abs(b - a) < eps:
            return (a + b) / 2, n, a, b
        if (x_1_flag == 1):
            x_1 = a + alpha * (b - a)
            func_x_1 = func(x_1)
            n += 1

        if (x_2_flag == 1):
            x_2 = b - alpha * (b - a)
            func_x_2 = func(x_2)
            n += 1

        if func_x_1 > func_x_2:
            a = x_1
            x_1 = x_2
            func_x_1 = func_x_2
            x_1_flag = 0
            x_2_flag = 1

        else:
            b = x_2
            x_2 = x_1
            func_x_2 = func_x_1
            x_1_flag = 1
            x_2_flag = 0



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


def method_of_steepest_descent(func, grad_func, eps):
    x_k = 0
    y_k = 0.314
    grad = grad_func(x_k, y_k)
    while(True):
        if((grad[0]**2 + grad[1]**2)**0.5 < eps):
            return x_k, y_k
        alpha_k = dichotomy_method(func, eps, grad, [x_k, y_k])[0]
        x_k -= alpha_k * grad[0]
        y_k -= alpha_k * grad[1]
        grad = grad_func(x_k, y_k)









x, y = np.meshgrid(np.linspace(-5, 1, 1000), np.linspace(-2, 1.5, 1000))

result = method_of_steepest_descent(f, grad_f, 10e-7)

print(f'[{result[0]}, {result[1]}]')


z = f(x, y)
plt.figure()
plt.contour(x, y, z, 20)
plt.scatter(result[0], result[1])
plt.show()

# fig = plt.figure(figsize=(7, 4))
# ax_3d = fig.add_subplot(projection='3d')
# ax_3d.plot_wireframe(x, y, z)
# plt.show()
