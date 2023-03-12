import math


def golden_section_method(func, a, b, eps):
    alpha = (3 - 5 ** 0.5) / 2
    n_eps = 2 * (math.log(eps) - math.log(b - a)) / math.log(1 - alpha)
    n_eps = math.ceil(n_eps)
    n_eps += (n_eps % 2)
    n = 0
    print(f'predictive n: {n_eps}')
    while True:
        if abs(b - a) < eps:
            return (a + b) / 2, n, a, b
        x_1 = a + alpha * (b - a)
        x_2 = b - alpha * (b - a)
        func_x_1 = func(x_1)
        func_x_2 = func(x_2)
        n += 2
        if func_x_1 > func_x_2:
            a = x_1
        if func_x_1 <= func_x_2:
            b = x_2


def dichotomy_method(func, a, b, eps):
    n_eps = 2 * (math.log(eps) - math.log(b - a)) / math.log(0.501)
    n_eps = math.ceil(n_eps)
    n_eps += (n_eps % 2)
    n = 0
    print(f'predictive n: {n_eps}')
    while True:
        if abs(b - a) < eps:
            return (a + b) / 2, n, a, b
        delta = (b - a) * 0.001
        x_1 = (a + b) / 2 - delta
        x_2 = (a + b) / 2 + delta
        func_x_1 = func(x_1)
        func_x_2 = func(x_2)
        n += 2
        if func_x_1 > func_x_2:
            a = x_1
        if func_x_1 <= func_x_2:
            b = x_2


def trial_points_method(func, a, b, eps):
    n_eps_min = 2 * (math.log(eps) - math.log(b - a)) / math.log(0.5)
    n_eps_min = math.ceil(n_eps_min)
    n_eps_min += (n_eps_min % 2)
    n_eps_max = 3 * (math.log(eps) - math.log(b - a)) / math.log(0.5)
    n_eps_max = math.ceil(n_eps_max)
    while n_eps_max % 3 != 0:
        n_eps_max += 1
    print(f'predictive n belongs [{n_eps_min}, {n_eps_max}]')
    n = 0
    while True:
        if abs(b - a) < eps:
            return (a + b) / 2, n, a, b
        quarter = (b - a) / 4
        x_1 = a + quarter
        x_2 = x_1 + quarter
        func_x_1 = func(x_1)
        func_x_2 = func(x_2)
        n += 2
        if func_x_1 < func_x_2:
            b = x_2
        else:
            x_3 = x_2 + quarter
            func_x_3 = func(x_3)
            n += 1
            if func_x_2 < func_x_3:
                a = x_1
                b = x_3
            else:
                a = x_2
