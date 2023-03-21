import math

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
    n_eps_min = (math.log(eps) - math.log(b - a)) / math.log(0.5)
    n_eps_min = math.ceil(n_eps_min)
    n_eps_min += 1
    n_eps_max = 2 * (math.log(eps) - math.log(b - a)) / math.log(0.5)
    n_eps_max = math.ceil(n_eps_max)
    n_eps_max += (n_eps_max % 2)
    while n_eps_max % 3 != 0:
        n_eps_max += 1
    print(f'predictive n belongs [{n_eps_min}, {n_eps_max}]')
    n = 0
    while True:
        if abs(b - a) < eps:
            return (a + b) / 2, n, a, b
        quarter = (b - a) / 4
        x_1 = a + quarter
        func_x_1 = func(x_1)
        n += 1
        if(n == 1):
            x_2 = x_1 + quarter
            func_x_2 = func(x_2)
            n += 1

        print(f'x_1={x_1}, x_2={x_2}, [{a}, {b}]')
        if func_x_1 < func_x_2:
            b = x_2
            x_2 = x_1
            func_x_2 = func_x_1
        else:
            x_3 = x_2 + quarter
            func_x_3 = func(x_3)
            n += 1
            if func_x_2 < func_x_3:
                a = x_1
                b = x_3
            else:
                a = x_2
                x_2 = x_3
                func_x_2 = func_x_3
