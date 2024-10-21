import numpy as np
import scipy.optimize as opt
from matplotlib import pyplot as plt


def avg_residual(N: int):
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

    error_forward = np.abs(y_forward - y_analytical)
    error_backward = np.abs(y_backward - y_analytical)
    error_central = np.abs(y_central - y_analytical)

    return h, np.mean(error_forward), np.mean(error_backward), np.mean(error_central)


forward_residual = []
backward_residual = []
central_residual = []

hs = []

for N in [11, 101, 1001, 10001]:
    h, forward, backward, central = avg_residual(N)
    hs.append(h)
    forward_residual.append(forward)
    backward_residual.append(backward)
    central_residual.append(central)

plt.plot(hs, forward_residual, label="Forward", marker="o")
plt.plot(hs, backward_residual, label="Backward", marker="*")
plt.plot(hs, central_residual, label="Central", marker="x")
plt.xscale("log")
plt.yscale("log")
plt.xlabel(r"$\ln h$")
plt.ylabel(r"$\ln (\text{average residual})$")
plt.legend()

plt.savefig("figs/prob1-2-2.png")


def predict_residual_log(hs, residua):
    linear = np.polyfit(np.log10(hs), np.log10(residua), 1)

    def predict(logh):
        return linear[0] * logh + linear[1]

    return predict


pred_forward_log = predict_residual_log(hs, forward_residual)
pred_backward_log = predict_residual_log(hs, backward_residual)
pred_central_log = predict_residual_log(hs, central_residual)

# 残差が 10^{-6} になるような h を求める
logh_forward = opt.root(
    lambda logh: pred_forward_log(logh) - np.log10(1e-6), np.log10(hs[0])
).x[0]
logh_backward = opt.root(
    lambda logh: pred_backward_log(logh) - np.log10(1e-6), np.log10(hs[0])
).x[0]
logh_central = opt.root(
    lambda logh: pred_central_log(logh) - np.log10(1e-6), np.log10(hs[0])
).x[0]

print(f"Forward: h = {10 ** logh_forward}, logh = {logh_forward}")
print(f"Backward: h = {10 ** logh_backward}, logh = {logh_backward}")
print(f"Central: h = {10 ** logh_central}, logh = {logh_central}")
