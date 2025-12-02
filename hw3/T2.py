import numpy as np
from scipy.special import roots_legendre

def f(x):
    return np.exp(-x) / x

# 梯形法则
def trapezoidal_rule(func, a, b, n):
    x = np.linspace(a, b, n+1)
    y = func(x)
    h = (b - a) / n
    integral = (h / 2) * (y[0] + 2 * np.sum(y[1:n]) + y[n])
    return integral

# 辛普森法则
def simpsons_rule(func, a, b, n):
    x = np.linspace(a, b, 2*n+1)
    y = func(x)
    h = (b - a) / (2 * n)
    integral = (h / 3) * (y[0] + 4 * np.sum(y[1:(2*n):2]) + 2 * np.sum(y[2:(2*n-1):2]) + y[2*n])
    return integral

# Gauss-Legendre 方法
def gauss_legendre_rule(func, a, b, n):
    [nodes, weights] = roots_legendre(n+1)
    # 变换节点到 [a, b]
    nodes = 0.5 * (b - a) * nodes + 0.5 * (a + b)
    weights = 0.5 * (b - a) * weights
    integral = np.sum(weights * func(nodes))
    return integral


# 题目
if __name__ == "__main__":
    a = 1
    b = 100

    # 梯形法则
    for n in [10, 100, 1000]:
        result = trapezoidal_rule(f, a, b, n)
        print(f"Trapezoidal Rule with {n+1} points: {result:.7f}")

    # 辛普森法则
    for n in [10, 100, 1000]:
        result = simpsons_rule(f, a, b, n)
        print(f"Simpson's Rule with {2*n+1} points: {result:.7f}")

    # Gauss-Legendre 方法
    for n in [10, 100]:
        result = gauss_legendre_rule(f, a, b, n)
        print(f"Gauss-Legendre Rule with {n+1} points: {result:.7f}")
