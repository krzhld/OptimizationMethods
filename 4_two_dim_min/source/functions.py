import numpy as np

SQRT_5 = 2.23606797749978


def f(x: float, y: float):
    return np.power(x, 2) + 5 * np.power(y, 2) + np.sin(4 * x + 5 * y) + 3 * x + 2 * y


def f_zoom(x: float, y: float):
    return np.power(x, 2) + np.power(y, 2) + np.sin(4 * x + SQRT_5 * y) + 3 * x + 2 * y / SQRT_5


def grad_f(x: float, y: float):
    return 2 * x + 4 * np.cos(4 * x + 5 * y) + 3, 10 * y + 5 * np.cos(4 * x + 5 * y) + 2


def grad_f_zoom(x: float, y: float):
    return 2 * x + 4 * np.cos(4 * x + SQRT_5 * y) + 3, 2 * y + SQRT_5 * np.cos(4 * x + SQRT_5 * y) + 2 / SQRT_5


def hess_f(x: float, y: float):
    return [[2 - 16 * np.sin(4 * x + 5 * y), -20 * np.sin(4 * x + 5 * y)],
            [-20 * np.sin(4 * x + 5 * y), 10 - 25 * np.sin(4 * x + 5 * y)]]


def hess_f_zoom(x: float, y: float):
    return [[2 - 16 * np.sin(4 * x + SQRT_5 * y), -4 * SQRT_5 * np.sin(4 * x + SQRT_5 * y)],
            [-4 * SQRT_5 * np.sin(4 * x + SQRT_5 * y), 2 - 5 * np.sin(4 * x + SQRT_5 * y)]]


def norm(x: float, y: float):
    return (x ** 2 + y ** 2) ** 0.5
