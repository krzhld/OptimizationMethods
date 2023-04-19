import numpy as np


def f_1(x: float, y: float):
    return np.power((x - 2), 2) + np.power((y + 3), 2)


def f_2(x: float, y: float):
    return x + y + 4 * np.sqrt(1 + 2 * np.power(x, 2) + 3 * np.power(y, 2))


def grad_f_1(x: float, y: float):
    return 2 * (x - 2), 2 * (y + 3)


def grad_f_2(x: float, y: float):
    return 1 + (8 * x) / np.sqrt(1 + 2 * x ** 2 + 3 * y ** 2), 1 + (12 * y) / np.sqrt(1 + 2 * x ** 2 + 3 * y ** 2)


def hess_f_1(x: float, y: float):
    return [[2, 0], [0, 2]]


def hess_f_2(x: float, y: float):
    return [[8 * (2 * x ** 2 + 3 * y ** 2 + 1) ** (-0.5) - 16 * x ** 2 * (2 * x ** 2 + 3 * y ** 2 + 1) ** (-1.5),
             -24 * x * y * (2 * x ** 2 + 3 * y ** 2 + 1) ** (-1.5)],
            [-24 * x * y * (2 * x ** 2 + 3 * y ** 2 + 1) ** (-1.5),
             12 * (2 * x ** 2 + 3 * y ** 2 + 1) ** (-0.5) - 36 * y ** 2 * (2 * x ** 2 + 3 * y ** 2 + 1) ** (-1.5)]]
