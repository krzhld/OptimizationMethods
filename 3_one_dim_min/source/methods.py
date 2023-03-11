import matplotlib.pyplot as plt
import numpy as np
import math

def golden_section_method(func, a, b, eps):
    alpha = (3 - 5 ** 0.5) / 2
    n = 0
    # n_eps = (math.log(eps)-math.log(b-a)) / math.log(1-alpha)
    # print(n_eps)
    while (1):
        if(abs(b - a) < eps):
            return (a + b) / 2, n
        x_1 = a + alpha * (b - a)
        x_2 = b - alpha * (b - a)
        func_x_1 = func(x_1)
        func_x_2 = func(x_2)
        n += 2
        if(func_x_1 > func_x_2):
            a = x_1
        if(func_x_1 <= func_x_2):
            b = x_2


def dichotomy_method(func, a, b, eps):
    n = 0
    # n_eps = (math.log(eps) - math.log(b - a)) / math.log(0.501)
    # print(n_eps)
    while(1):
        if (abs(b - a) < eps):
            return (a + b) / 2, n
        delta = (b - a) * 0.001
        x_1 = (a + b) / 2 - delta
        x_2 = (a + b) / 2 + delta
        func_x_1 = func(x_1)
        func_x_2 = func(x_2)
        n += 2
        if (func_x_1 > func_x_2):
            a = x_1
        if (func_x_1 <= func_x_2):
            b = x_2