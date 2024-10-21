import numpy as np
from matplotlib import pyplot as plt
from pathlib import Path

cwd = Path(__file__).parent

N = 1001

a = -np.pi
b = np.pi


def delta(i, j):
    return 1 if i == j else 0


ham = np.zeros((N, N))
x = np.linspace(a, b, N)
h = (b - a) / (N - 1)

for i in range(N):
    for j in range(N):
        ham[i, j] = -0.5 * (delta(i + 1, j) - 2 * delta(i, j) + delta(i - 1, j)) / h**2

eigenvalues, eigenvectors = np.linalg.eigh(ham)

print("numerical:", eigenvalues[:3])

# normalize
wf = eigenvectors / np.sqrt(h)


def analytical_wf(i):
    fnc = np.sin if i % 2 == 0 else np.cos
    return 1/np.sqrt(np.pi) * fnc(i * x / 2.0)

def analytical_energy(i):
    return i**2 / 8.0

print("analytical:", [analytical_energy(i) for i in range(1, 4)])


plt.plot(x, wf[:, 0], label="Ground state")
plt.plot(x, wf[:, 1], label="First excited state")
plt.plot(x, wf[:, 2], label="Second excited state")
plt.plot(x, analytical_wf(1), label="Analytical ground state", linestyle="--")
plt.plot(x, analytical_wf(2), label="Analytical first excited state", linestyle="--")
plt.plot(x, analytical_wf(3), label="Analytical second excited state", linestyle="--")
plt.xlabel("$x$")
plt.ylabel(r"$\psi$")
plt.legend()

plt.savefig(cwd / "figs" / "prob1-3.png")

# numerical
for i in range(3):
    for j in range(3):
        print(f"N: overlap({i}, {j}):", np.dot(wf[:, i], wf[:, j])*h)

# analytical
for i in range(3):
    for j in range(3):
        print(f"A: overlap({i}, {j}):", np.dot(analytical_wf(i + 1), analytical_wf(j + 1)) * h)

