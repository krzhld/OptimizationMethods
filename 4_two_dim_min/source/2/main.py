import math
import numpy as np
import matplotlib.pyplot as plt

def f_1(x, y):
    return np.power((x - 2), 2) + np.power((y + 3), 2)


def f_2(x, y):
    return x + y + 4 * np.sqrt(1 + 2 * np.power(x,2) + 3 * np.power(y,2))


def f(x, y):
    return np.power(x,2) + 5 * np.power(y,2) + np.sin(4 * x + 5 * y) + 3 * x + 2 * y


def grad_f(x, y):
    return 2 * x + 4 * math.cos(4 * x + 5 * y) + 3, 10 * y + 5 * math.cos(4 * x + 5 * y) + 2

def hess_f(x,y):
    return [[2 - 16 * math.sin(4 * x + 5 * y), -20 * math.sin(4 * x + 5 * y)], [-20 * math.sin(4 * x + 5 * y), 10 - 25 * math.sin(4 * x + 5 * y)]]

def grad_f_1(x, y):
    return 2 * (x - 2), 2 * (y + 3)


def grad_f_2(x, y):
    return 1 + (8 * x) / math.sqrt(1 + 2 * x**2 + 3 * y**2) , 1 + (12 * y) / math.sqrt(1 + 2 * x**2 + 3 * y**2) 


alpha = (3 - 5 ** 0.5) / 2

def norm(x, y):
    return (x**2 + y**2)**0.5 


def mult_matrix_and_vector(A, b) -> list:
    return A[0][0] * b[0] + A[0][1] * b[1], A[1][0] * b[0] + A[1][1] * b[1]

def mult_vectors(v1 : list, v2 : list) -> float:
    return v1[0] * v2[0] + v1[1] * v2[1]

def newton_method(func, grad_func, hess_func, eps):
    x_k = 15
    y_k = 2
    alpha_0 = 1
    delta = 0.5 
    lambd = 0.5 

    while(True):
        plt.scatter(x_k, y_k)
        grad = grad_func(x_k, y_k)
        norm_grad = (grad[0]**2 + grad[1]**2)**0.5
        if(norm_grad < eps):
            return x_k, y_k 
        
        # print(f'{norm_grad}')
        d_k = np.array(mult_matrix_and_vector(np.linalg.inv(hess_func(x_k,y_k)), grad)) 
        alpha_k = 0.01 
        alpha_k = alpha_0

        f = func(x_k - alpha_k * d_k[0], y_k - alpha_k * d_k[1])
        f_k = func(x_k, y_k)
        while(f - f_k > -delta * alpha_k * mult_vectors(grad, d_k)):
            alpha_k *= lambd

            f = func(x_k - alpha_k * d_k[0], y_k - alpha_k * d_k[1])
            

        # alpha_0 = alpha_k 
        x_k -= alpha_k * d_k[0] 
        y_k -= alpha_k * d_k[1]




x, y = np.meshgrid(np.linspace(-4, 1, 1000), np.linspace(-1.5, 1, 1000))




# print(f'[{result[0]}, {result[1]}]')

h = np.array(hess_f(1, 1))
hinv = np.array(np.linalg.inv(h))

print(f'{h}')
print(f'{hinv}')
print(f'{h * hinv}')

z = f(x, y)
plt.figure()
plt.contour(x, y, z, 50)
result = newton_method(f, grad_f, hess_f, 10e-6)
plt.scatter(result[0], result[1])
plt.show()
