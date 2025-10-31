# 样条函数在计算机绘图中的运用
import numpy as np
import matplotlib.pyplot as plt
from T3 import cubic_spline_interpolation

# 我们将尝试用三次样条函数平滑地连接若干个二维空间中已知的点。考虑二维空间的一系列点 (xi, yi), i = 0, 1, ···, n。我们现在希望按照顺序（由 0 到 n）将它们平滑地连接起来。一个方便的办法是引入一个连续参数 t ∈ [0, n]，取节点为ti = 0, 1, ···, n，然后分别建立两个样条函数：S∆(X; t) 和 S∆(Y ; t)，它们分别满足S∆(X; t) = xi; S∆(Y, t) = yi。这两个样条函数可以看作是 (x(t), y(t)) 的内插近似。因此绘制参数曲线 (x(t), y(t)) 的问题就化为求出两个样条函数并将它们画出的问题。

def plot_spline_curve(points, t_points, num_points=1000):
    """
    使用三次样条插值绘制二维空间中的平滑曲线，返回插值点的坐标。
    :param points: 已知数据点的列表，格式为 [(x0, y0), (x1, y1), ..., (xn, yn)]
    :param t_points: 参数点列表，通常为 [0, 1, 2, ..., n]
    :param num_points: 用于绘制曲线的插值点数量
    :return: 插值点的坐标列表 (x_fine, y_fine)，分别为 x 和 y 坐标
    """
    points = np.array(points)
    x_values = points[:, 0]
    y_values = points[:, 1]
    n = len(points)

    # 构造参数 t
    t_values = np.array(t_points)

    # 生成用于绘制的细分 t 值
    t_fine = np.linspace(0, n - 1, num_points)

    # 计算样条插值
    x_fine = cubic_spline_interpolation(t_values, x_values, t_fine)
    y_fine = cubic_spline_interpolation(t_values, y_values, t_fine)

    return x_fine, y_fine


# 题目测试
if __name__ == "__main__":
    # 第二题，选取 ϕ = tπ/4, t = 0, 1, 2, 3, 4, 5, 6, 7, 8 这 9 个点，给出 xt = r(ϕ) cos ϕ 和 yt = r(ϕ) sin ϕ 的数值。将这些数值作为精确的数值列在一个表里。我们考虑的函数是著名的心形线（cardioid）。它的极坐标方程是r(ϕ) = 2a(1 − cos ϕ) = 1 − cos ϕ，为了方便起见我们取了 2a = 1。给出过这 9 个点的两个三次样条函数 S∆(X; t) 和 S∆(Y ; t)。
    print('第二题\n输出中的小数保留6位有效数字。')
    a = 0.5
    t_points = np.arange(0, 9)
    phi_values = t_points * (np.pi / 4)
    r_values = 2 * a * (1 - np.cos(phi_values))
    x_values = r_values * np.cos(phi_values)
    y_values = r_values * np.sin(phi_values)

    # 将这些x和y的数值作为精确的数值列在一个表里
    # print("t\tphi\t\trx\t\try")
    # for t, phi, x, y in zip(t_points, phi_values, x_values, y_values):
    #     print(f"{t}\t{phi:.6f}\t{x:.6f}\t{y:.6f}")

    # 给出过这 9 个点的两个三次样条函数 S∆(X; t) 和 S∆(Y ; t)
    print("S(X; t):")
    piecewise_expressions_x = cubic_spline_interpolation(t_points, x_values, need_piecewise_expression=True)
    print("\nS(Y; t):")
    piecewise_expressions_y = cubic_spline_interpolation(t_points, y_values, need_piecewise_expression=True)


    # 第三题
    print('\n第三题\n如图所示。')

    # 计算样条插值曲线
    points = list(zip(x_values, y_values))
    x_fine, y_fine = plot_spline_curve(points, t_points)

    # 原曲线
    t_original = np.linspace(0, 2 * np.pi, 1000)
    r_original = 2 * a * (1 - np.cos(t_original))
    x_original = r_original * np.cos(t_original)
    y_original = r_original * np.sin(t_original)

    # 绘图
    plt.figure(figsize=(8, 8))
    plt.plot(x_original, y_original, label='Original Curve', color='blue', linewidth=2)
    plt.plot(x_fine, y_fine, label='Spline Interpolated Curve', color='red', linestyle='--', linewidth=2)
    plt.scatter(x_values, y_values, color='green', marker='o', s=100, label='Data Points')
    for i, (x, y) in enumerate(zip(x_values, y_values)):
        plt.text(x, y, f'P{i}', fontsize=12, ha='right')
    plt.title('Cubic Spline Interpolation of Cardioid')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.axis('equal')
    plt.legend()
    plt.grid()
    plt.show()


    # 第四题，简要说明为什么这个算法可以平滑地连接所有的点 (这实际上是很多画图软件中 spline 曲线所采用的算法)。
    print('\n第四题')

    print("三次样条插值通过在每个相邻数据点之间构造三次多项式，并确保这些多项式在数据点处连续且具有连续的一阶和二阶导数，从而实现了平滑连接。\n这种方法不仅保证了曲线通过所有给定点，还确保了曲线的光滑性，避免了不连续或尖锐的转折。\n因此，三次样条插值是绘图软件中常用的平滑曲线生成方法。")
