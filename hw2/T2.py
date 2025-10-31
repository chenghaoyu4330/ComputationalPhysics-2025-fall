# 插值,包括有理分式内插法和拉格朗日插值法
import numpy as np
import matplotlib.pyplot as plt

def rational_function_interpolation_phi(x_values, y_values):
    """
    计算使用有理分式内插法时用于构造连分式的函数phi的值。对y = f(x) 进行有理分式插值，给定插值点xi以及yi = f(xi),输出用于构造连分式的phi(x0, x1), phi(x0, x1, x2), ...的值
    :param x_values: 已知数据点的 x 坐标
    :param y_values: 已知数据点的 y 坐标
    :return: 用于构造连分式的phi的值列表
    """
    n = len(x_values)
    if n != len(y_values):
        raise ValueError("x_values 和 y_values 必须具有相同的长度")

    # 初始化一个二维数组来存储中间计算结果
    phi = np.zeros((n, n))

    # 填充 phi 的第一列y_values
    for i in range(n):
        phi[i][0] = y_values[i]

    # 递推计算 phi 的其他列
    for i in range(1, n):  # phi 的行
        for j in range(1, i+1):  # phi 的列
            numerator = x_values[j-1] - x_values[i]
            denominator = phi[j-1][j-1] - phi[i][j-1]
            phi[i][j] = numerator / denominator

    # 提取对角线上的值作为结果
    result = [phi[i][i] for i in range(1, n)]
    return result

def rational_function_interpolation(x_values, y_values, x):
    """
    使用有理分式内插法计算在给定点处的插值
    :param x_values: 已知数据点的 x 坐标
    :param y_values: 已知数据点的 y 坐标
    :param x: 需要插值的点列表
    :return: 在指定点处的插值结果列表
    """
    n = len(x_values)
    if n != len(y_values):
        raise ValueError("x_values 和 y_values 必须具有相同的长度")

    # 计算 phi 的值
    phi_values = rational_function_interpolation_phi(x_values, y_values)

    # 将y_values的第一个值插入到phi_values的开头
    phi_values.insert(0, y_values[0])

    # 计算插值结果
    results = []
    for xi in x:
        result = phi_values[-1]  # phi(x0, x1, ..., xn-1)
        for i in range(n-1, 0, -1):
            result = phi_values[i-1] + (xi - x_values[i-1]) / result
        results.append(result)

    return results


# 拉格朗日插值法
def lagrange_interpolation(x_values, y_values, x):
    """
    使用拉格朗日插值法计算在点 x 处的插值
    :param x_values: 已知数据点的 x 坐标
    :param y_values: 已知数据点的 y 坐标
    :param x: 需要插值的点列表
    :return: 插值结果
    """
    n = len(x_values)
    if n != len(y_values):
        raise ValueError("x_values 和 y_values 必须具有相同的长度")

    results = []
    for xi in x:
        result = 0
        for i in range(n):
            term = y_values[i]
            for j in range(n):
                if j != i:
                    term *= (xi - x_values[j]) / (x_values[i] - x_values[j])
            result += term
        results.append(result)
    return results


# 题目测试
if __name__ == "__main__":
    # 第一小问
    print("第一小问：有理分式内插法中，用于构造连分式的phi值")
    func = lambda x: np.tan(x)  # 需要插值的原函数
    x_values = np.array([0, 0.3, 0.6, 0.9, 1.2, 1.5])  # 构造插值函数所用的节点
    y_values = func(x_values)  # 原函数的真实值
    result = rational_function_interpolation_phi(x_values, y_values)
    for val in result:
        print(f"{val:.6f}")

    # 第二小问
    print("\n第二小问：有理分式内插法和拉格朗日插值法的比较")
    x = np.linspace(0, 1.57, 500)  # 需要插值的点
    rational_results = rational_function_interpolation(x_values, y_values, x)
    lagrange_results = lagrange_interpolation(x_values, y_values, x)

    # print(f'有理分式内插法在 x=1.57 处的插值结果: {rational_results[np.searchsorted(x, 1.57)]:.6f}')
    # print(f'拉格朗日插值法在 x=1.57 处的插值结果: {lagrange_results[np.searchsorted(x, 1.57)]:.6f}')
    # print(f'真实值 tan(1.57) = {func(1.57):.6f}')

    # 绘制插值结果的对比图,取前500个点
    dots = 500
    plt.figure(figsize=(10, 6))
    plt.plot(x[:dots], rational_results[:dots], label='rational', linestyle='--')
    plt.plot(x[:dots], lagrange_results[:dots], label='lagrange', linestyle=':')
    plt.plot(x[:dots], func(x[:dots]), label='true value', linestyle='-')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('compare')
    plt.legend()
    plt.grid()
    plt.show()

    # 分别计算两种方法的结果在各插值点处的误差的平方和
    func_values = func(x)
    rational_error = np.sum((np.array(rational_results) - func_values) ** 2)
    lagrange_error = np.sum((np.array(lagrange_results) - func_values) ** 2)
    print(f"有理分式内插法的误差平方和: {rational_error:.6e}")
    print(f"拉格朗日插值法的误差平方和: {lagrange_error:.6e}")

    # 分析
    print("\n分析：")
    print("从误差平方和可以看出，有理分式内插法的误差明显小于拉格朗日插值法，说明有理分式内插法在本例中表现更好，能够更准确地逼近函数tan(x)的值。")