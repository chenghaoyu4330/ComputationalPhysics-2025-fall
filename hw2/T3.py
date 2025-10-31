# Runge效应
import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
from T2 import lagrange_interpolation


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
    使用 Chebyshev 近似进行函数计算：f(x) = -c0/2 + Σ (ck * Tk(x))，其中 Tk(x) 是第 k 个 Chebyshev 多项式，ck 是对应的系数
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
        for j in range(1, n):
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


# 三次样条函数插值
def cubic_spline_interpolation(x_values, y_values, x=None, need_piecewise_expression=False):
    """
    使用三次样条插值法计算在给定点处的插值
    :param x_values: 已知数据点的 x 坐标
    :param y_values: 已知数据点的 y 坐标
    :param x: 需要插值的点列表
    :param need_piecewise_expression: 是否返回分段函数表达式
    :return: 在指定点处的插值结果列表或分段函数表达式
    """
    n = len(x_values)
    if n != len(y_values):
        raise ValueError("x_values 和 y_values 必须具有相同的长度")
    if n <= 2:
        raise ValueError("至少需要三个数据点进行三次样条插值")

    # 定义区间长度hi，以及构造线性方程组求解三次样条函数的矩m（二阶导数），其中mi满足μi * m_(i-1) + 2 * mi + λi * m_(i+1) = di，i = 1, 2, ..., n-2，其中μi =
    # h_(i-1) / (hi + h_(i-1))，λi = hi / (hi + h_(i-1))，di = 6 * (y_(i+1) / h_(i-1)*(h_(i-1)+hi) - yi / hi*h_(i-1) +
    # y_(i-1) / h_(i-1)*(h_(i-1)+hi))
    h = [x_values[i + 1] - x_values[i] for i in range(n - 1)]
    A = np.zeros((n-2, n-2))
    b = np.zeros(n-2)
    # 取自然边界条件m_0 = m_(n-1) = 0
    m = np.zeros(n)
    for i in range(1, n - 1):
        if i != 1:
            A[i - 1][i - 2] = h[i - 1] / (h[i - 1] + h[i])  # μi
        A[i - 1][i - 1] = 2
        if i != n - 2:
            A[i - 1][i] = h[i] / (h[i - 1] + h[i])  # λi
        b[i - 1] = 6 * ((y_values[i + 1] / h[i]) - (y_values[i] * (h[i - 1] + h[i]) / (h[i - 1] * h[i])) + (y_values[i - 1] / h[i - 1])) / (h[i - 1] + h[i])

    # 求解线性方程组
    m[1:n-1] = np.linalg.solve(A, b)

    # 计算插值结果
    if x is not None:
        results = []
        for xi in x:
            # 找到xi所在的区间[xi, xi+1]
            for i in range(n - 1):
                if x_values[i] <= xi <= x_values[i + 1]:
                    break
            hi = h[i]
            Ai = -m[i] / (6 * hi)
            Bi = m[i + 1] / (6 * hi)
            Ci = (y_values[i + 1] - y_values[i]) / hi - hi * (m[i + 1] - m[i]) / 6
            Di = y_values[i] - m[i] * hi ** 2 / 6
            dxi = xi - x_values[i]
            dxi1 = xi - x_values[i + 1]
            fx = Ai * dxi1 ** 3 + Bi * dxi ** 3 + Ci * dxi + Di
            results.append(fx)

        return results

    if need_piecewise_expression:
        # 返回分段函数表达式
        piecewise_expressions = []
        sym = 't'
        x_sym = sp.symbols(sym)
        for i in range(n - 1):
            hi = h[i]
            Ai = -m[i] / (6 * hi)
            Bi = m[i + 1] / (6 * hi)
            Ci = (y_values[i + 1] - y_values[i]) / hi - hi * (m[i + 1] - m[i]) / 6
            Di = y_values[i] - m[i] * hi ** 2 / 6
            dxi = x_sym - x_values[i]
            dxi1 = x_sym - x_values[i + 1]
            fx = Ai * dxi1 ** 3 + Bi * dxi ** 3 + Ci * dxi + Di
            print(f"S_{i}({sym}) = {sp.simplify(fx).evalf(6)} for {sym} in [{x_values[i]}, {x_values[i + 1]}]")
            piecewise_expressions.append((sp.simplify(fx), x_values[i], x_values[i + 1]))
        return piecewise_expressions


# 题目测试
if __name__ == "__main__":
    # 第一题
    print("第一题：Runge效应，采用拉格朗日插值法")
    func = lambda x: 1 / (1 + 25 * x**2)
    x_values = np.linspace(-1, 1, 21)  # 21个等距节点
    y_values = func(x_values)

    # 需要插值的点，包括节点和每个小段的中点
    x_interp = []
    for i in range(len(x_values) - 1):
        x_interp.append(x_values[i])
        x_interp.append((x_values[i] + x_values[i + 1]) / 2)
    x_interp.append(x_values[-1])

    # 计算拉格朗日插值结果
    y_interp_lagrange = lagrange_interpolation(x_values, y_values, x_interp)

    # 计算真实函数值
    y_true = func(np.array(x_interp))

    # 输出结果
    print(f"{'x':>10} {'f(x)':>15} {'P20(x)':>15} {'|f(x)-P20(x)|':>20}")
    for xi, fi, pi in zip(x_interp, y_true, y_interp_lagrange):
        print(f"{xi:10.6f} {fi:15.6f} {pi:15.6f} {abs(fi - pi):20.6e}")

    # 绘图
    plt.figure(figsize=(10, 6))
    # x_plot = np.linspace(-1, 1, 400)
    x_plot = np.array(x_interp)
    y_plot_true = func(x_plot)
    y_plot_lagrange = lagrange_interpolation(x_values, y_values, x_plot)
    plt.plot(x_plot, y_plot_true, label='True Function', linestyle='-')
    plt.plot(x_plot, y_plot_lagrange, label='Lagrange Interpolation', linestyle='--')
    plt.scatter(x_interp, y_true, color='red', s=10, label='Interpolation Points')
    plt.title('Prob1\nRunge Phenomenon: Lagrange Interpolation vs True Function')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend()
    plt.grid()
    plt.show()

    # 分析
    print('可以看到，使用均匀分布节点的拉格朗日插值法在区间端点附近出现了较大的误差，即Runge现象。')


    # 第二题
    print("\n第二题：Chebyshev 多项式插值法")

    n_chebyshev = 21
    a, b = -1, 1

    # 需要插值的点，包括节点和每个小段的中点
    x_interp_chebyshev = []
    chebyshev_nodes = Chebyshev_nodes(a, b, n_chebyshev)
    for i in range(len(chebyshev_nodes) - 1):
        x_interp_chebyshev.append(chebyshev_nodes[i])
        x_interp_chebyshev.append((chebyshev_nodes[i] + chebyshev_nodes[i + 1]) / 2)
    x_interp_chebyshev.append(chebyshev_nodes[-1])

    # 计算 Chebyshev 多项式插值结果
    y_interp_chebyshev = Chebyshev_poly_interpolation(func, a, b, n_chebyshev, x_interp_chebyshev)

    # 计算真实函数值
    y_true_chebyshev = func(np.array(x_interp_chebyshev))

    # 输出结果
    print(f"{'x':>10} {'f(x)':>15} {'Chebyshev P(x)':>20} {'|f(x)-P(x)|':>20}")
    for xi, fi, pi in zip(x_interp_chebyshev, y_true_chebyshev, y_interp_chebyshev):
        print(f"{xi:10.6f} {fi:15.6f} {pi:20.6f} {abs(fi - pi):20.6e}")

    # 绘图
    plt.figure(figsize=(10, 6))
    x_plot_chebyshev = np.array(x_interp_chebyshev)
    y_plot_true_chebyshev = func(x_plot_chebyshev)
    y_plot_chebyshev = Chebyshev_poly_interpolation(func, a, b, n_chebyshev, x_plot_chebyshev)
    plt.plot(x_plot_chebyshev, y_plot_true_chebyshev, label='True Function', linestyle='-')
    plt.plot(x_plot_chebyshev, y_plot_chebyshev, label='Chebyshev Interpolation', linestyle='--')
    plt.scatter(x_interp_chebyshev, y_true_chebyshev, color='red', s=10, label='Interpolation Points')
    plt.title('Prob2\nChebyshev Interpolation vs True Function')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend()
    plt.grid()

    # 分析
    print('可以看到，Chebyshev 多项式插值法在整个区间上的误差显著小于均匀分布节点的拉格朗日插值法，避免了Runge现象。')


    # 第三题
    print("\n第三题：f(x) = |x| 的插值比较")
    func_abs = lambda x: abs(x)

    # 均匀分布节点
    x_values_abs = np.linspace(-1, 1, 21)
    y_values_abs = func_abs(x_values_abs)

    # 需要插值的点，包括节点和每个小段的中点
    x_interp_abs = []
    for i in range(len(x_values_abs) - 1):
        x_interp_abs.append(x_values_abs[i])
        x_interp_abs.append((x_values_abs[i] + x_values_abs[i + 1]) / 2)
    x_interp_abs.append(x_values_abs[-1])

    # 拉格朗日插值
    y_interp_lagrange_abs = lagrange_interpolation(x_values_abs, y_values_abs, x_interp_abs)

    # Chebyshev插值
    y_interp_chebyshev_abs = Chebyshev_poly_interpolation(func_abs, a, b, n_chebyshev, x_interp_abs)

    # 真实值
    y_true_abs = func_abs(np.array(x_interp_abs))

    # 输出结果
    print(f"{'x':>10} {'f(x)':>15} {'Lagrange P(x)':>20} {'Chebyshev P(x)':>20} {'|f(x)-Lagrange P(x)|':>20} {'|f(x)-Chebyshev P(x)|':>20}")
    for xi, fi, pi_lag, pi_cheb in zip(x_interp_abs, y_true_abs, y_interp_lagrange_abs, y_interp_chebyshev_abs):
        print(f"{xi:10.6f} {fi:15.6f} {pi_lag:20.6f} {pi_cheb:20.6f} {abs(fi - pi_lag):20.6e} {abs(fi - pi_cheb):20.6e}")

    # 绘图，分为左右两部分，左边为[-1, 1]上的整体图，右边为x接近0的局部放大图
    plt.figure(figsize=(12, 6))
    # 整体图
    plt.subplot(1, 2, 1)
    x_plot_abs = np.array(x_interp_abs)
    y_plot_true_abs = func_abs(x_plot_abs)
    y_plot_lagrange_abs = np.array(lagrange_interpolation(x_values_abs, y_values_abs, x_plot_abs))
    y_plot_chebyshev_abs = np.array(Chebyshev_poly_interpolation(func_abs, a, b, n_chebyshev, x_plot_abs))
    plt.plot(x_plot_abs, y_plot_true_abs, label='True Function', linestyle='-')
    plt.plot(x_plot_abs, y_plot_lagrange_abs, label='Lagrange Interpolation', linestyle='--')
    plt.plot(x_plot_abs, y_plot_chebyshev_abs, label='Chebyshev Interpolation', linestyle=':')
    plt.title('Prob3\nf(x) = |x| Interpolation Comparison (Overall)')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend()
    plt.grid()
    # 局部放大图
    plt.subplot(1, 2, 2)
    zoom_indices = np.where((x_plot_abs >= -0.2) & (x_plot_abs <= 0.2))
    plt.plot(x_plot_abs[zoom_indices], y_plot_true_abs[zoom_indices], label='True Function', linestyle='-')
    plt.plot(x_plot_abs[zoom_indices], y_plot_lagrange_abs[zoom_indices], label='Lagrange Interpolation', linestyle='--')
    plt.plot(x_plot_abs[zoom_indices], y_plot_chebyshev_abs[zoom_indices], label='Chebyshev Interpolation', linestyle=':')
    plt.title('Prob3\nf(x) = |x| Interpolation Comparison (Zoomed in at x≈0)')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.show()

    # 分析
    print('对于 f(x) = |x|，在 x ≈ 0 附近，均匀分布节点的拉格朗日插值法表现出较大的误差，而 Chebyshev 插值法则表现较好，但仍然存在一定的误差。\n这是因为 f(x) 在 x=0 处不可导，导致插值多项式难以准确逼近该点的行为。')


    # 第四题
    print("\n第四题：三次样条插值法")

    # 需要插值的点，这次利用 61 个点的三次样条函数
    x_interp_spline = np.linspace(-1, 1, 61)

    # 三次样条插值，构造样条函数的节点仍然是均匀分布的21个节点
    y_interp_spline = cubic_spline_interpolation(x_values, y_values, x_interp_spline)

    # 真实值
    y_true_spline = func(np.array(x_interp_spline))

    # 输出结果
    print(f"{'x':>10} {'f(x)':>15} {'Spline P(x)':>20} {'|f(x)-P(x)|':>20}")
    for xi, fi, pi in zip(x_interp_spline, y_true_spline, y_interp_spline):
        print(f"{xi:10.6f} {fi:15.6f} {pi:20.6f} {abs(fi - pi):20.6e}")

    # 绘图
    plt.figure(figsize=(10, 6))
    x_plot_spline = np.array(x_interp_spline)
    y_plot_true_spline = func(x_plot_spline)
    y_plot_spline = cubic_spline_interpolation(x_values, y_values, x_plot_spline)
    plt.plot(x_plot_spline, y_plot_true_spline, label='True Function', linestyle='-')
    plt.plot(x_plot_spline, y_plot_spline, label='Cubic Spline Interpolation', linestyle='--')
    plt.scatter(x_interp_spline, y_true_spline, color='red', s=10, label='Interpolation Points')
    plt.title('Prob4\nCubic Spline Interpolation vs True Function')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend()
    plt.grid()
    plt.show()

    # 分析
    print('三次样条插值法在整个区间上表现出良好的逼近效果，误差较小，且没有出现Runge现象。')


    # 第五题
    print("\n第五题：Lebesgue 数的计算与分析")

    def lebesgue_function(x_values, a, b, dots=1000):
        """
        计算 Lebesgue 函数值
        :param x_values: 插值节点
        :param a: 区间起点
        :param b: 区间终点
        :param dots: 计算点的数量
        :return: x 点列表和对应的 Lebesgue 函数值列表
        """
        n = len(x_values)
        x = np.linspace(a, b, dots)
        L_values = []
        for xi in x:
            L_xi = 0
            for i in range(n):
                term = 1
                for j in range(n):
                    if j != i:
                        term *= abs((xi - x_values[j]) / (x_values[i] - x_values[j]))
                L_xi += term
            L_values.append(L_xi)
        return x, L_values

    # 计算均匀分布节点的 Lebesgue 函数
    x_lebesgue_uniform, L_uniform = lebesgue_function(x_values, -1, 1)
    # 计算 Chebyshev 节点的 Lebesgue 函数
    chebyshev_nodes_21 = Chebyshev_nodes(-1, 1, 21)
    x_lebesgue_chebyshev, L_chebyshev = lebesgue_function(chebyshev_nodes_21, -1, 1)

    # 分别绘图，左图为均匀分布节点，右图为Chebyshev节点
    plt.figure(figsize=(12, 6))
    # 均匀分布节点
    plt.subplot(1, 2, 1)
    plt.plot(x_lebesgue_uniform, L_uniform, label='Uniform Nodes', color='blue')
    plt.title('Prob5\nLebesgue Function - Uniform Nodes')
    plt.xlabel('x')
    plt.ylabel('Lebesgue Function L(x)')
    plt.grid()
    # Chebyshev节点
    plt.subplot(1, 2, 2)
    plt.plot(x_lebesgue_chebyshev, L_chebyshev, label='Chebyshev Nodes', color='orange')
    plt.title('Prob5\nLebesgue Function - Chebyshev Nodes')
    plt.xlabel('x')
    plt.ylabel('Lebesgue Function L(x)')
    plt.grid()
    plt.tight_layout()
    plt.show()

    # 分析
    print('Lebesgue 数衡量了插值过程中的误差放大效应。对于均匀分布节点，Lebesgue 数在区间端点附近显著增大，表明误差可能被放大，导致Runge现象的出现。\n而Chebyshev节点的Lebesgue数相对较小且分布更均匀，说明其插值过程中的误差放大效应较小，更稳定。因此，选择合适的插值节点对于控制插值误差至关重要。')
    print('Lebesgue 数确实存在下限，且对于任意节点分布，其Lebesgue 数的下限均大于等于1。\n这是因为在插值过程中，至少会有一个基函数在某些点处取值为1。\n例如对于基函数Li(x)，当x等于xi时，Li(xi)=1，则Lebesgue函数L(x)在该点处至少为1。因此Lebesgue数有下限，且至少为1。')