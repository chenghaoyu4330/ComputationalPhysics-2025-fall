import numpy as np
import sympy as sp
import matplotlib.pyplot as plt

# 生成 Chebyshev 节点
def Chebyshev_nodes(a, b, n):
    """
    生成区间 [a, b] 上的 Chebyshev 节点
    :param a: 区间起点
    :param b: 区间终点
    :param n: 节点数量
    :return: Chebyshev 节点列表
    """
    nodes = []
    for k in range(n):
        x_k = 0.5 * (a + b) + 0.5 * (b - a) * np.cos((2 * k + 1) * np.pi / (2 * n))
        nodes.append(x_k)
    return nodes


# Chebyshev 多项式插值
def Chebyshev_poly_interpolation(f, a, b, n, x):
    """
    使用 Chebyshev 多项式进行函数近似计算
    :param f: 目标函数
    :param a: 区间起点
    :param b: 区间终点
    :param n: 节点数量
    :param x: 需要插值的点列表，保证在 [-1, 1] 范围内
    :return: 在指定点处的插值结果列表
    """
    # 生成 Chebyshev 节点
    nodes = Chebyshev_nodes(a, b, n)

    # 计算节点处的函数值
    y_values = [f(xi) for xi in nodes]

    # 使用Clenshaw算法计算Chebyshev展开：对于Chebyshev多项式Tk(x)，满足 Tk+1(x) = 2xTk(x) - Tk-1(x)，则构造新函数bk(x)，满足：bn+1(x) = 0,
    # bn+2(x) = 0, bk(x) = 2x * bk+1(x) - bk+2(x) + ck，最终f(x) = c0 + x * b1(x) - b2(x)，因此需要计算ck，bk(x) 计算系数ck
    c = []
    for k in range(n):
        sum_ck = 0
        for j in range(0, n):
            sum_ck += y_values[j] * np.cos(k * (2 * j + 1) * np.pi / (2 * n))
        ck = (2 / n) * sum_ck
        # 特殊处理c0
        if k == 0:
            ck = ck / 2
        c.append(ck)

    # 计算bk(x)并计算插值结果
    results = []
    for xi in x:
        b_kplus2 = 0
        b_kplus1 = 0
        for k in range(n - 1, 0, -1):
            b_k = 2 * xi * b_kplus1 - b_kplus2 + c[k]
            b_kplus2 = b_kplus1
            b_kplus1 = b_k
        fx = xi * b_kplus1 - b_kplus2 + c[0]
        results.append(fx)
    return results


# 题目
if __name__ == "__main__":
    # 定义目标函数
    f = lambda x: np.sin(2*np.pi*x) / (1+x*x)

    # 定义区间和节点数量
    a = 2
    b = 5
    n = [5, 10, 15]

    # 定义插值点
    x_vals = np.linspace(a, b, 400)
    x_mapped = (2 * x_vals - (a + b)) / (b - a)  # 映射到 [-1, 1]
    f_exact = f(x_vals)

    # 绘图，同时把插值的点也画出来
    for ni in n:
        ni = ni + 1  # n阶插值需要n+1个节点
        plt.figure(figsize=(12, 8))
        f_interp = Chebyshev_poly_interpolation(f, a, b, ni, x_mapped)
        plt.plot(x_vals, f_interp, label=f'Chebyshev Interpolation n={ni}')
        plt.plot(x_vals, f_exact, 'k--', label='Exact Function', linewidth=1)
        plt.plot(Chebyshev_nodes(a, b, ni), [f(x) for x in Chebyshev_nodes(a, b, ni)], 'ro', label='Chebyshev Nodes')
        plt.title('Chebyshev Polynomial Interpolation')
        plt.xlabel('x')
        plt.ylabel('f(x)')
        plt.legend()
        plt.grid()
        plt.show()