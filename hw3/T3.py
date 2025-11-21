import numpy as np
import sympy as sp
import matplotlib.pyplot as plt

def f(x, k):
    return np.exp(x) * np.cos(k * x)

def exact_integral(k):
    x = sp.symbols('x')
    integral = sp.integrate(sp.exp(x) * sp.cos(k * x), (x, 0, 1))
    return float(integral)

def simpsons_integral(func, a, b, n, k):
    if n % 2 == 0:
        n += 1  # Simpson's rule requires an odd number of points
    x = np.linspace(a, b, n)
    y = func(x, k)
    h = (b - a) / (n - 1)
    integral = (h / 3) * (y[0] + 4 * np.sum(y[1:n-1:2]) + 2 * np.sum(y[2:n-2:2]) + y[n-1])
    return integral

def I0_I2_J1(k, h):
    """
    计算积分 I0, I2, J1
    """
    kh = k * h
    if k == 0:
        I0 = 2.0 * h
        I2 = 2.0 * h**3 / 3.0
        J1 = 0.0
    else:
        I0 = 2.0 * np.sin(kh) / k
        I2 = 2.0 * (h**2 * k**2 * np.sin(kh) + 2.0 * h * k * np.cos(kh) - 2.0 * np.sin(kh)) / (k**3)
        J1 = 2.0 * (-h * k * np.cos(kh) + np.sin(kh)) / (k**2)
    return I0, I2, J1

def composite_filon_exp_cos_vectorized(k, N):
    """
    使用 Filon 方法计算积分 ∫[0,1] exp(x) * cos(kx) dx
    """
    h = 1.0 / (2.0 * N)
    # 中心点 x_mid = (2j+1)*h, j=0..N-1
    j = np.arange(N, dtype=np.float64)
    x_mid = (2.0 * j + 1.0) * h

    # 三点的函数值
    y0 = np.exp(x_mid - h)   # left
    y1 = np.exp(x_mid)       # mid
    y2 = np.exp(x_mid + h)   # right

    # 解出 A,B,C
    A = y1
    B = (y2 - y0) / (2.0 * h)
    C = (y0 + y2 - 2.0 * y1) / (2.0 * h * h)

    # 计算积分 I0, I2, J1
    I0, I2, J1 = I0_I2_J1(k, h)

    # 计算最终结果
    cos_kxm = np.cos(k * x_mid)
    sin_kxm = np.sin(k * x_mid)
    block_vals = cos_kxm * (A * I0 + C * I2) - sin_kxm * (B * J1)
    total = np.sum(block_vals)
    return total


# 题目
if __name__ == "__main__":
    a = 0
    b = 1
    ks = [10**3, 10**4, 10**5, 10**6]
    alphas = [10**(-2), 10**(-1), 10**0, 10**1, 10**2]

    results = {}

    print("第一题")
    print("k\t\t\tN\t\t\tNumerical\t\t\tExact\t\t\tError")
    for k in ks:
        exact_value = exact_integral(k)
        results[k] = []
        for alpha in alphas:
            N = int(alpha * k)
            numerical_value = simpsons_integral(f, a, b, N, k)
            error = abs(numerical_value - exact_value)
            results[k].append((N, numerical_value, exact_value, error))
            print(f"{k:.0e}\t\t{N:.0e}\t\t{numerical_value:.10e}\t\t{exact_value:.10e}\t\t{error:.10e}")

    # 可视化误差变化趋势
    plt.figure(figsize=(10, 6))
    for k in ks:
        Ns = [res[0] for res in results[k]]
        errors = [res[3] for res in results[k]]
        relative_error = [err / abs(results[k][i][2]) for i, err in enumerate(errors)]
        plt.loglog(Ns, relative_error, marker='o', label=f'k={k}')
    plt.xlabel('Number of Subintervals N')
    plt.ylabel('Relative Error')
    plt.title('Error Trend in Simpson\'s Numerical Integration')
    plt.legend()
    plt.grid(True, which="major", ls="--")
    plt.show()



    print("第三题")
    print("k\t\t\tN\t\t\tNumerical\t\t\tExact\t\t\tError")
    results = []
    for k in ks:
        exact_value = exact_integral(k)
        for alpha in alphas:
            N = int(alpha * k)
            numerical_value = composite_filon_exp_cos_vectorized(k, N)
            error = abs(numerical_value - exact_value)
            results.append((k, N, numerical_value, exact_value, error))
            print(f"{k:.0e}\t\t{N:.0e}\t\t{numerical_value:.10e}\t\t{exact_value:.10e}\t\t{error:.10e}")

    # 可视化误差变化趋势
    plt.figure(figsize=(10, 6))
    for k in ks:
        Ns = [res[1] for res in results if res[0] == k]
        errors = [res[4] for res in results if res[0] == k]
        relative_error = [err / abs([res[3] for res in results if res[0] == k][i]) for i, err in enumerate(errors)]
        plt.loglog(Ns, relative_error, marker='o', label=f'k={k}')
    plt.xlabel('Number of Subintervals N')
    plt.ylabel('Relative Error')
    plt.title('Error Trend in Filon\'s Numerical Integration')
    plt.legend()
    plt.grid(True, which="major", ls="--")
    plt.show()
