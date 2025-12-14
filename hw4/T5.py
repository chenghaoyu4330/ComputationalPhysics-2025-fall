import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm, chi2

def run_simulation(N_samples=1000, N_experiments=10000, n_bins=10):
    # 理论上每个区间的期望频数
    expected_count = N_samples / n_bins
    
    print(f"正在运行 {N_experiments} 次实验，每次 {N_samples} 个样本...")
    
    # 生成所有随机数
    all_randoms = np.random.rand(N_experiments, N_samples)
    
    # 区间划分 [0,0.1), [0.1,0.2), ..., [0.9,1.0], bin 索引 0 到 9
    bin_indices = np.floor(all_randoms * n_bins).astype(int)

    # 处理 rand 正好为 1.0 的边缘情况
    bin_indices[bin_indices == n_bins] = n_bins - 1
    
    # 统计每个实验在各区间的频数
    bin_counts = np.zeros((N_experiments, n_bins))
    for b in range(n_bins):
        # 统计落在bin索引为 b 的样本在每行中出现的次数
        bin_counts[:, b] = np.sum(bin_indices == b, axis=1)
        
    # (a) 部分数据，区间 [0.6, 0.7] 对应 bin 索引 6
    counts_06_07 = bin_counts[:, 6]
    proportions_06_07 = counts_06_07 / N_samples
    
    # (b) 部分数据，每个实验的卡方统计量: sum (O - E)^2 / E
    chi_squared = np.sum((bin_counts - expected_count)**2 / expected_count, axis=1)
    
    return proportions_06_07, chi_squared

def plot_part_a(proportions, N_samples):
    plt.figure(figsize=(10, 6))
    
    # 比例的直方图
    count, bins, _ = plt.hist(proportions, bins=50, density=True, alpha=0.6, color='b', label='Simulation')
    
    # 理论正态分布
    # 均值 = p = 0.1
    # 方差 = p(1-p)/N = 0.1 * 0.9 / 1000 = 0.00009
    # 标准差 = sqrt(0.00009) 约 0.009487
    mu = 0.1
    sigma = np.sqrt(0.1 * 0.9 / N_samples)
    
    x = np.linspace(min(bins), max(bins), 100)
    plt.plot(x, norm.pdf(x, mu, sigma), linewidth=2, color='r', label='Normal theory')
    
    # 图像中避免中文
    plt.title('Proportion in [0.6, 0.7] (Part a)')
    plt.xlabel('Proportion')
    plt.ylabel('Density')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.savefig('hw4/assets/random_dist_a.png')
    print("已保存 (a) 部分图表至 hw4/assets/random_dist_a.png")

def plot_part_b(chi_sq_values, n_bins):
    plt.figure(figsize=(10, 6))
    
    # 卡方值的直方图
    count, bins, _ = plt.hist(chi_sq_values, bins=50, density=True, alpha=0.6, color='g', label='Simulation')
    
    # 理论卡方分布
    # 自由度 = bin数量 - 1 = 9
    df = n_bins - 1
    x = np.linspace(min(bins), max(bins), 100)
    plt.plot(x, chi2.pdf(x, df), linewidth=2, color='r', label=r'$\chi^2$ theory (df=9)')
    
    # 图像中避免中文
    plt.title(r'$\chi^2$ statistic (Part b)')
    plt.xlabel(r'$\chi^2$ value')
    plt.ylabel('Density')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.savefig('hw4/assets/random_chi2_b.png')
    print("已保存 (b) 部分图表至 hw4/assets/random_chi2_b.png")

    # 均值与方差检验
    print(f"卡方统计量均值: {np.mean(chi_sq_values):.4f} (理论值: {df})")
    print(f"卡方统计量方差: {np.var(chi_sq_values):.4f} (理论值: {2*df})")

class LCG16807:
    def __init__(self, seed=1):
        self.state = seed
        self.a = 16807
        self.m = 2147483647  # 2^31 - 1
        
    def next(self):
        self.state = (self.a * self.state) % self.m
        return self.state
    
    def next_float(self):
        return self.next() / self.m

def verify_lcg():
    print("\n--- (c) 16807 随机数生成器验证 ---")
    generator = LCG16807(seed=1)
    
    target_n = 10000
    final_val = 0
    
    # 需要 x(10000)，x(0)=seed，因此迭代 10000 次
    for i in range(target_n):
        final_val = generator.next()
        
    print(f"x({target_n}) = {final_val}")
    expected = 1043618065
    
    if final_val == expected:
        print("验证成功，结果与题目一致。")
    else:
        print(f"验证失败，期望值为 {expected}")

if __name__ == "__main__":
    N_samples = 1000
    N_experiments = 10000
    n_bins = 10

    # 计算 (a) 与 (b)
    props, chi_sqs = run_simulation(N_samples, N_experiments, n_bins)
    
    # 绘制 (a)
    plot_part_a(props, N_samples)
    
    # 绘制 (b)
    plot_part_b(chi_sqs, n_bins)
    
    # 验证 (c)
    verify_lcg()
