# 第5题
import numpy as np
import matplotlib.pyplot as plt


# 定义函数
def f(x):
    return np.exp(5 * x)

# 精确导数值
x0 = -1
exact_derivative = 5 * np.exp(5 * x0)
print(f"精确导数值: f'(-1) = {exact_derivative:.13f}")

# 生成不同步长
h_values = np.logspace(-14, -2, 10000)  # 从10^-14到10^-2，共10000个点

# 计算数值微分和误差
numerical_derivatives = []
errors = []
relative_errors = []

for h in h_values:
    # 差分法求微分
    numerical_deriv = (f(x0 + h) - f(x0)) / h
    numerical_derivatives.append(numerical_deriv)

    # 绝对误差
    error = abs(numerical_deriv - exact_derivative)
    errors.append(error)

    # 相对误差
    relative_error = error / abs(exact_derivative)
    relative_errors.append(relative_error)

# 转换为numpy数组以便处理
numerical_derivatives = np.array(numerical_derivatives)
errors = np.array(errors)
relative_errors = np.array(relative_errors)

# 找到最小误差对应的步长
min_error_idx = np.argmin(errors)
optimal_h = h_values[min_error_idx]
min_error = errors[min_error_idx]

print(f"\n最优步长: h = {optimal_h:.2e}")
print(f"最小绝对误差: {min_error:.2e}")
print(f"此时数值微分值: {numerical_derivatives[min_error_idx]:.13f}")


# 创建图形
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))

# 左图：数值微分值随步长变化
ax1.loglog(h_values, numerical_derivatives, 'b-', linewidth=2, label='Numerical Derivative Value')
ax1.axhline(y=exact_derivative, color='r', linestyle='--', linewidth=2, label='Exact Value')
ax1.set_xlabel('Step size h')
ax1.set_ylabel('Derivative Value')
ax1.set_title('Numerical Derivative Value vs Step Size')
ax1.grid(True, alpha=0.3)
ax1.legend()

# 右图：误差随步长变化
ax2.loglog(h_values, errors, 'g-', linewidth=2, label='Absolute Error')
ax2.axvline(x=optimal_h, color='r', linestyle='--', linewidth=2,
           label=f'Optimal h={optimal_h:.2e}')
ax2.set_xlabel('Step size h')
ax2.set_ylabel('Absolute Error')
ax2.set_title('Absolute Error vs Step Size')
ax2.grid(True, alpha=0.3)
ax2.legend()

plt.tight_layout()
plt.show()

# 相对误差图
plt.figure(figsize=(10, 6))
plt.loglog(h_values, relative_errors, 'purple', linewidth=2, label='Relative Error')
plt.axvline(x=optimal_h, color='r', linestyle='--', linewidth=2,
           label=f'Optimal h={optimal_h:.2e}')
plt.xlabel('Step size h')
plt.ylabel('Relative Error')
plt.title('Relative Error vs Step Size')
plt.grid(True, alpha=0.3)
plt.legend()
plt.show()

# 误差分析
print("\n=== 误差分析 ===")

# 截断误差主导区域 (h较大)
large_h_indices = np.where(h_values > 1e-4)[0][-5:]
print("\n大步长区域 (截断误差主导):")
for i in large_h_indices:
    print(f"h = {h_values[i]:.2e}, 误差 = {errors[i]:.2e}")

# 舍入误差主导区域 (h很小)
small_h_indices = np.where(h_values < 1e-12)[0][:5]
print("\n很小步长区域 (舍入误差主导):")
for i in small_h_indices:
    print(f"h = {h_values[i]:.2e}, 误差 = {errors[i]:.2e}")

# 最小值附近
optimal_region = np.where((h_values > optimal_h*0.95) & (h_values < optimal_h*1.05))[0]
print(f"\n最小值附近 (h ≈ {optimal_h:.2e}):")
for i in optimal_region[::7]:  # 每隔7个点取样
    print(f"h = {h_values[i]:.2e}, 误差 = {errors[i]:.2e}")

# 误差分析
print("\n=== 误差分析 ===")
print("差分公式的截断误差: O(h)，具体为M·h/2，M为|f''(x)|在区间[x0, x0+h]上的上界")
print("舍入误差: O(ε/h)，具体为2ε·|f(x0)|/h，ε为机器精度，含有2是因为计算了两次f(x)")
print("总误差 = |截断误差| + |舍入误差| ≈ M·h/2 + 2ε·|f(x0)|/h")
print("可以看到，总误差是h的函数，且在某个h值处取得最小值；当h较大时，截断误差主导，总误差随h增大而增大；当h较小时，舍入误差主导，总误差随h减小而增大")
print("通过对总误差表达式求导并令导数为零，可以找到最优步长 h_opt ≈ √(4ε·|f(x0)|/M)")

# 机器精度
machine_epsilon = np.finfo(float).eps
print(f"\n机器精度 ε ≈ {machine_epsilon:.2e}")

estimated_optimal_h = np.sqrt(4 * machine_epsilon * abs(f(x0)) / (25 * abs(f(x0))))
print(f"代入公式可得理论估计的最优步长 h_opt ≈ {estimated_optimal_h:.2e}")