# 求解线性方程组，包括高斯消元法和雅可比迭代法
import numpy as np

def gaussian_elimination(A, b):
    """
    高斯消元法求解线性方程组 Ax = b
    :param A: 系数矩阵
    :param b: 常数向量
    :return: 解向量 x
    """
    n = len(b)
    # 构造增广矩阵
    Ab = np.hstack([A, b.reshape(-1, 1)])

    # 前向消元
    for i in range(n):
        # 寻找主元
        max_row = np.argmax(np.abs(Ab[i:, i])) + i
        Ab[[i, max_row]] = Ab[[max_row, i]]  # 交换行

        # 消元
        for j in range(i + 1, n):
            factor = -Ab[j, i] / Ab[i, i]
            Ab[j] += factor * Ab[i]

    # 回代求解
    x = np.zeros(n)
    for i in range(n - 1, -1, -1):
        x[i] = (Ab[i, -1] - np.dot(Ab[i, i + 1:n], x[i + 1:n])) / Ab[i, i]

    return x


def jacobi_iteration(A, b, x0=None, tol=1e-10, max_iterations=1000):
    """
    雅可比迭代法求解线性方程组 Ax = b，核心思想：x^(k+1) = D^(-1)(b - (L + U)x^(k))，其中 D 是 A 的对角矩阵，L 是 A 的严格下三角矩阵，U 是 A 的严格上三角矩阵。
    :param A: 系数矩阵
    :param b: 常数向量
    :param x0: 初始猜测解
    :param tol: 收敛容差
    :param max_iterations: 最大迭代次数
    :return: 解向量 x
    """
    n = len(b)
    if x0 is None:
        x0 = np.zeros(n)

    x = np.copy(x0)

    for iteration in range(max_iterations):
        x_new = np.zeros(n)
        for i in range(n):
            sum1 = np.dot(A[i, :i], x[:i])  # L部分，即Σ_(j=0)^(i-1) a_ij * x_j
            sum2 = np.dot(A[i, i + 1:], x[i + 1:])  # U部分，即Σ_(j=i+1)^(n-1) a_ij * x_j
            x_new[i] = (b[i] - sum1 - sum2) / A[i, i]  # 迭代公式：x_i^(k+1)=(b_i - Σ_(j≠i) a_ij * x_j^(k)) / a_ii

        # 检查收敛性
        if np.linalg.norm(x_new - x, ord=np.inf) < tol:
            return x_new

        x = x_new
        # print(x)

    raise ValueError(f"雅可比迭代未能在指定的迭代次数{max_iterations}次内收敛")

# 题目测试
if __name__ == "__main__":
    A = np.array([[6, -2, 1],
                  [1, 5, 2],
                  [-1, 1, 4]], dtype=float)
    b = np.array([5, 17, 13], dtype=float)

    print("使用高斯消元法求解:")
    x_gaussian = gaussian_elimination(A, b)
    for i, xi in enumerate(x_gaussian):
        print(f"x_{i+1} = {xi:.6f}")

    print("使用雅可比迭代法求解:")
    x_jacobi = jacobi_iteration(A, b)
    for i, xi in enumerate(x_jacobi):
        print(f"x_{i+1} = {xi:.6f}")