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

def hess_f_1(x,y):
    return [[2, 0], [0, 2]]

def grad_f_2(x, y):
    return 1 + (8 * x) / math.sqrt(1 + 2 * x**2 + 3 * y**2) , 1 + (12 * y) / math.sqrt(1 + 2 * x**2 + 3 * y**2) 

def hess_f_2(x,y):
    return [[8 * (2 * x**2 + 3 * y**2 + 1)**(-0.5) - 16 * x**2 * (2 * x**2 + 3 * y**2 + 1)**(-1.5), -24 * x * y * (2 * x**2 + 3 * y**2 + 1)**(-1.5) ], [-24 * x * y * (2 * x**2 + 3 * y**2 + 1)**(-1.5), 12 * (2 * x**2 + 3 * y**2 + 1)**(-0.5) - 36 * y**2 * (2 * x**2 + 3 * y**2 + 1)**(-1.5)]]

alpha = (3 - 5 ** 0.5) / 2

def norm(x, y):
    return (x**2 + y**2)**0.5 


def mult_matrix_and_vector(A, b) -> list:
    return A[0][0] * b[0] + A[0][1] * b[1], A[1][0] * b[0] + A[1][1] * b[1]

def mult_vectors(v1 : list, v2 : list) -> float:
    return v1[0] * v2[0] + v1[1] * v2[1]

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

def method_of_steepest_descent_gold_sect(func, grad_func, eps):
    x_k = -1
    y_k = -1
    numb_of_iter = 0 
    grad = grad_func(x_k, y_k)
    while(True):
        norm_grad = (grad[0]**2 + grad[1]**2)**0.5
        if(norm_grad < eps):
            return x_k, y_k, numb_of_iter
        result_find_alpha_k = golden_section_method(func, eps, grad, [x_k, y_k])
        alpha_k = result_find_alpha_k[0]

        x_k_1 = x_k - alpha_k * grad[0]
        y_k_1 = y_k - alpha_k * grad[1] 
        
        x_k = x_k_1
        y_k = y_k_1
        grad = grad_func(x_k, y_k)
        numb_of_iter += 1



def newton_method(func, grad_func, hess_func, eps):
    x_k = -1
    y_k = -1
    alpha_0 = 1
    delta = 0.5 
    lambd = 0.5 
    numb_of_iter = 0

    while(True):
        grad = grad_func(x_k, y_k)
        norm_grad = (grad[0]**2 + grad[1]**2)**0.5
        if(norm_grad < eps):
            return x_k, y_k, numb_of_iter 
        
        print(f'{norm_grad}')
        d_k = np.array(mult_matrix_and_vector(np.linalg.inv(hess_func(x_k,y_k)), grad)) 
        alpha_k = alpha_0

        f = func(x_k - alpha_k * d_k[0], y_k - alpha_k * d_k[1])
        f_k = func(x_k, y_k)
        while(f - f_k > -delta * alpha_k * mult_vectors(grad, d_k)):
            alpha_k *= lambd

            f = func(x_k - alpha_k * d_k[0], y_k - alpha_k * d_k[1])
            
        x_k -= alpha_k * d_k[0] 
        y_k -= alpha_k * d_k[1]
        numb_of_iter += 1




x, y = np.meshgrid(np.linspace(-4, 1, 1000), np.linspace(-1.5, 1, 1000))

# x, y = np.meshgrid(np.linspace(-5, 5, 1000), np.linspace(-5, 5, 1000))


# print(f'[{result[0]}, {result[1]}]')



z = f(x, y)
plt.figure()
plt.contour(x, y, z, 50)
result = newton_method(f, grad_f, hess_f, 10e-4)
plt.scatter(result[0], result[1])
plt.show()

resultt = method_of_steepest_descent_gold_sect(f, grad_f, 10e-4)
print(f'{result}')


epsilons = [10e-3, 10e-4, 10e-5, 10e-6, 10e-7, 10e-8, 10e-9]
result_1, result_2 = [], []

for eps in epsilons:
    result_1.append(method_of_steepest_descent_gold_sect(f, grad_f, eps)[2])
    result_2.append(newton_method(f, grad_f, hess_f, eps)[2])


plt.figure()
plt.semilogx(epsilons, result_1, linewidth=0.5, label = "method of steepest descent")
plt.semilogx(epsilons, result_2, linewidth=0.5, label = "newton method")
plt.legend()
plt.xlabel("eps")
plt.ylabel("number of iterations")
plt.title("Зависимость числа итераций от точности")
plt.show()