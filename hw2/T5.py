import numpy as np
import matplotlib.pyplot as plt
from T1 import gaussian_elimination


def cholesky(M, eps=1e-10):
    """
    Cholesky 分解, 将对称正定矩阵M分解为L·L^T，L为下三角矩阵
    :param M: 对称正定矩阵
    :param eps: 用于判断正定性的阈值
    :return: L矩阵
    """
    n = M.shape[0]
    L = np.zeros((n, n), dtype=float)

    for i in range(n):
        for j in range(i + 1):
            # 计算对角线元素（i == j）
            if i == j:
                # L[i,i] = sqrt(M[i,i] - sum(L[i,k]^2 for k=0..i-1))
                sum_sq = np.sum(L[i, :i] ** 2)
                diag_val = M[i, i] - sum_sq
                if diag_val < eps:  # 若被开方数非正，矩阵非正定
                    raise ValueError("矩阵非正定，无法进行 Cholesky 分解")
                L[i, i] = np.sqrt(diag_val)
            # 计算非对角线元素（i > j）
            else:
                # L[i,j] = (M[i,j] - sum(L[i,k]·L[j,k] for k=0..j-1)) / L[j,j]
                sum_prod = np.sum(L[i, :j] * L[j, :j])
                L[i, j] = (M[i, j] - sum_prod) / L[j, j]

    return L


def forward_substitution(L, b):
    """
    求解下三角方程组 L·x = b
    :param L: 下三角矩阵
    :param b: 右侧向量
    :return: x: 解向量
    """
    n = L.shape[0]
    x = np.zeros(n, dtype=float)
    for i in range(n):
        # x[i] = (b[i] - sum(L[i,k]·x[k] for k=0..i-1)) / L[i,i]
        sum_val = np.sum(L[i, :i] * x[:i])
        x[i] = (b[i] - sum_val) / L[i, i]
    return x


def backward_substitution(LT, b):
    """
    求解上三角方程组 LT·x = b
    :param LT: 上三角矩阵
    :param b: 右侧向量
    :return: x: 解向量
    """
    n = LT.shape[0]
    x = np.zeros(n, dtype=float)
    for i in range(n - 1, -1, -1):  # 从最后一行向前遍历
        # x[i] = (b[i] - sum(LT[i,k]·x[k] for k=i+1..n-1)) / LT[i,i]
        sum_val = np.sum(LT[i, i + 1:] * x[i + 1:])
        x[i] = (b[i] - sum_val) / LT[i, i]
    return x


def matrix_inverse(M):
    """
    对称正定矩阵求逆（基于Cholesky分解）
    :param M: 对称正定矩阵
    :return: M的逆矩阵
    """
    n = M.shape[0]
    # 1. Cholesky分解
    L = cholesky(M)
    LT = L.T  # L的转置

    # 2. 求解L·Y = I，得到Y（每列对应单位矩阵的一列）
    Y = np.zeros((n, n))
    for j in range(n):
        e_j = np.zeros(n)
        e_j[j] = 1.0  # 单位向量（第j个元素为1）
        Y[:, j] = forward_substitution(L, e_j)  # 解第j列方程

    # 3. 求解LT·X = Y，得到X = M^{-1}
    M_inv = np.zeros((n, n))
    for j in range(n):
        M_inv[:, j] = backward_substitution(LT, Y[:, j])  # 解第j列方程

    return M_inv


def weighted_poly_fit(x, y, sigma, degree):
    """
    加权最小二乘法多项式拟合
    :param x: 自变量数据点
    :param y: 因变量数据点
    :param sigma: 每个数据点的标准差（误差）
    :param degree: 多项式阶数
    :return:
        coeffs: 拟合系数 [a0, a1, ..., am]
        coeffs_err: 系数误差 [σ_a0, σ_a1, ..., σ_am]
    """
    n = len(x)
    m = degree  # 多项式阶数

    # 构建矩阵M和向量b
    M = np.zeros((m + 1, m + 1))
    b = np.zeros(m + 1)

    for j in range(m + 1):
        for k in range(m + 1):
            # 计算M[j][k] = Σ(x_i^(j+k) / σ_i²)
            M[j, k] = np.sum(x ** (j + k) / sigma ** 2)

        # 计算b[j] = Σ(y_i * x_i^j / σ_i²)
        b[j] = np.sum(y * x ** j / sigma ** 2)

    # 求解线性方程组 M·a = b
    # coeffs = np.linalg.solve(M, b)
    coeffs = gaussian_elimination(M, b)

    # 计算协方差矩阵 (M的逆)
    # cov_matrix = np.linalg.inv(M)
    cov_matrix = matrix_inverse(M)

    # 系数误差是协方差矩阵对角线元素的平方根
    coeffs_err = np.sqrt(np.diag(cov_matrix))

    return coeffs, coeffs_err


def print_results(degree, coeffs, coeffs_err):
    """
    打印拟合结果
    :param degree: 多项式阶数
    :param coeffs: 拟合系数
    :param coeffs_err: 系数误差
    """
    print(f"{degree}次多项式拟合结果:")
    for i in range(degree + 1):
        print(f"a_{i} = {coeffs[i]:.6f} ± {coeffs_err[i]:.6f}")
    print()


def plot_fit(x, y, sigma, coeffs_list, degrees):
    """
    绘制数据点和拟合曲线
    :param x: 自变量数据点
    :param y: 因变量数据点
    :param sigma: 每个数据点的标准差（误差）
    :param coeffs_list: 拟合系数列表，每个元素对应一个多项式阶数
    :param degrees: 多项式阶数列表
    """
    plt.figure(figsize=(10, 6))
    plt.errorbar(x, y, yerr=sigma, fmt='o', label='data', ecolor='r', capsize=5)

    # 生成平滑曲线用于绘图
    x_smooth = np.linspace(min(x), max(x), 1000)

    for coeffs, degree in zip(coeffs_list, degrees):
        # 计算拟合曲线
        y_fit = np.polyval(np.flip(coeffs), x_smooth)
        plt.plot(x_smooth, y_fit, label=f'degree={degree}')

    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend()
    plt.title('Fit')
    plt.grid(True)
    plt.show()


# 题目测试
if __name__ == "__main__":
    print("第五题：最小二乘法（带权多项式拟合）")

    # 表格数据
    x = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9], dtype=float)
    y = np.array([1.2, 2.8, 4.5, 7.1, 9.8, 13.5, 17.2, 21.9, 27.1, 33.2], dtype=float)
    sigma = np.array([0.12, 0.24, 0.14, 0.33, 0.90, 0.11, 0.45, 0.43, 0.59, 0.10], dtype=float)

    degrees = [1, 2, 3]
    coeffs_list = []

    # 分别进行1次、2次、3次多项式拟合
    for degree in degrees:
        coeffs, coeffs_err = weighted_poly_fit(x, y, sigma, degree)
        coeffs_list.append(coeffs)
        print_results(degree, coeffs, coeffs_err)

    # 绘制拟合结果
    plot_fit(x, y, sigma, coeffs_list, degrees)
