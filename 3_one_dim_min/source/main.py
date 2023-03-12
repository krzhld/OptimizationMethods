import matplotlib.pyplot as plt
import numpy as np
import methods


def f(x):
    result = (x * (x - 1) * (x - 2) * (x - 3)) / 8
    return result


a = 2.1
b = 3.2
eps = 0.0001

answer_golden_section, n_golden_section = methods.golden_section_method(f, a, b, eps)
print('golden section method')
print(f'answer: {answer_golden_section}, n: {n_golden_section}')
print()

answer_dichotomy, n_dichotomy = methods.dichotomy_method(f, a, b, eps)
print('dichotomy method')
print(f'answer: {answer_dichotomy}, n: {n_dichotomy}')
print()

answer_trial_points, n_trial_points = methods.trial_points_method(f, a, b, eps)
print('trial points method')
print(f'answer: {answer_trial_points}, n: {n_trial_points}')
print()

x = np.linspace(a, b, 100)
y = f(x)

plt.close('all')
plt.figure()
plt.plot(x, y, 'b-')
plt.plot(answer_golden_section, f(answer_golden_section), 'ro--', linewidth=2, markersize=10)
plt.title(f'golden section method, eps = {eps}, n = {n_golden_section}')
plt.xlabel('x-axis')
plt.ylabel('y-axis')
plt.grid()

plt.figure()
plt.plot(x, y, 'b-')
plt.plot(answer_dichotomy, f(answer_dichotomy), 'ro--', linewidth=2, markersize=10)
plt.title(f'dichotomy method, eps = {eps}, n = {n_dichotomy}')
plt.xlabel('x-axis')
plt.ylabel('y-axis')
plt.grid()

plt.figure()
plt.plot(x, y, 'b-')
plt.plot(answer_trial_points, f(answer_trial_points), 'ro--', linewidth=2, markersize=10)
plt.title(f'trial points method, eps = {eps}, n = {n_trial_points}')
plt.xlabel('x-axis')
plt.ylabel('y-axis')
plt.grid()

# plt.show()
