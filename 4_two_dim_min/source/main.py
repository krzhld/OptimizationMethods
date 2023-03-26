import numpy as np
import matplotlib.pyplot as plt


a = 5
b = 4
c = 5


# def func(x):
#     return x[0] ** 2 + a * x[1] ** 2 + m.sin(b * x[0] + c * x[1]) + 3 * x[0] + 2 * x[1]
def func(x1, x2):
    return x1 ** 2 + a * x2 ** 2 + np.sin(b * x1 + c * x2) + 3 * x1 + 2 * x2


def grad_func(x):
    return [2 * x[0] + b * np.cos(b * x[0] + c * x[1]) + 3, 2 * a * x[1] + c * np.cos(b * x[0] + c * x[1]) + 2]


x = np.arange(-4, 2, 0.05)
y = np.arange(-2, 2, 0.05)
x_grid, y_grid = np.meshgrid(x, y)

z_grid = func(x_grid, y_grid)

plt.contour(x_grid, y_grid, z_grid, levels=30)
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(x_grid, y_grid, z_grid, cmap='inferno')

plt.show()
