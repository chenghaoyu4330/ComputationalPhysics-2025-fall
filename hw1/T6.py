# 第6题
import numpy as np
import math


def psi_series(x, epsilon=1e-6):
    """
    计算级数 ψ(x)的近似值，保证截断误差小于 epsilon

    参数:
        x: 级数参数
        epsilon: 允许的最大截断误差

    返回:
        (近似值, 实际使用的项数)
    """
    # 当 x=0 时，单独考虑其情况，剩余项 R_n= ∑[k=n+1→∞] 1/k^2 < ∫[n,∞] 1/t^2 dt = 1/n < epsilon
    # 解不等式得到 n > 1/epsilon
    if x == 0:
        # 当 x=0 时，实际上级数为 ∑ 1/k^2 = π^2/6
        n_min = int(1 / epsilon) + 1
        partial_sum = 0.0
        for k in range(1, n_min + 1):
            partial_sum += 1 / (k * (k + x))
        return partial_sum, n_min

    # 当x≠0时，剩余项 R_n= ∑[k=n+1→∞] 1/(k(k+x)) < ∫[n,∞] 1/(t(t+x)) dt = (1/x) * ln(1 + x/n) < epsilon
    # 解不等式得到 n > x / (exp(x * epsilon) - 1)

    n_min = int(x / (math.exp(x * epsilon) - 1)) + 1

    # 对于大x，需要的项数很少，但为了安全我们设置一个下限
    n_min = max(n_min, 1000)  # 至少计算1000项

    # 计算部分和
    partial_sum = 0.0
    for k in range(1, n_min + 1):
        partial_sum += 1 / (k * (k + x))

    return partial_sum, n_min


def main():
    # 给定的x值
    x_values = [0.0, 0.5, 1.0, math.sqrt(2), 10.0, 100.0, 300.0]
    # 给定的截断误差要求
    epsilon = 1e-6

    print(f"级数求和结果 (截断误差 < {epsilon}):")
    print("=" * 50)

    results = []
    for i, x in enumerate(x_values):
        psi_approx, n_used = psi_series(x, epsilon)
        results.append((x, psi_approx, n_used))
        print(f"{i + 1}. x = {x:.6f}, psi(x) = {psi_approx:.15f}, 使用的项数: {n_used}")

    return results


# 主函数
if __name__ == "__main__":
    results = main()