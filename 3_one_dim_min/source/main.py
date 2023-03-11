import matplotlib.pyplot as plt
import numpy as np


def f(x):
    result = (x * (x - 1) * (x - 2) * (x - 3)) / 8
    return result


plt.close('all')
a = 2.1
b = 3.2
x = np.linspace(a, b, 100)
y = f(x)
plt.plot(x, y, 'b-')
plt.show()
