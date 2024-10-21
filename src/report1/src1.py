import numpy as np
from matplotlib import pyplot as plt


N = 1001

func = np.cos
x = np.linspace(0, 2 * np.pi, N)
h = x[1] - x[0]

y_forward = np.zeros(N)
y_backward = np.zeros(N)
y_central = np.zeros(N)

for i in range(1, N - 1):
    y_forward[i] = (func(x[i + 1]) - func(x[i])) / h
    y_backward[i] = (func(x[i]) - func(x[i - 1])) / h
    y_central[i] = (func(x[i + 1]) - func(x[i - 1])) / (2 * h)

y_analytical = -np.sin(x)

plt.title(r"$N = 1001, h = 6.28 \times 10^{-3}$")
plt.plot(x, y_forward, label="Forward")
plt.plot(x, y_backward, label="Backward")
plt.plot(x, y_central, label="Central")
plt.plot(x, y_analytical, label="Analytical")
plt.xlabel("$x$")
plt.ylabel("$y$")
plt.legend()
plt.savefig("figs/prob1-1.png")


error_forward = np.abs(y_forward - y_analytical)
error_backward = np.abs(y_backward - y_analytical)
error_central = np.abs(y_central - y_analytical)

plt.figure()
plt.title(r"$N = 1001, h = 6.28 \times 10^{-3}$")
plt.plot(x, error_forward, label="Forward")
plt.plot(x, error_backward, label="Backward")
plt.plot(x, error_central, label="Central")
plt.xlabel("$x$")
plt.ylabel("Error")
plt.legend()

plt.savefig("figs/prob1-2-1.png")
