import numpy as np

def construct_matrix(N):
    """
    构建一维原子链周期性边界条件下的矩阵 A。
    -A_ij = delta_{i, j-1} + delta_{i, j+1} - 2 delta_{i, j}
    因此 A_ij = 2 delta_{i, j} - delta_{i, j-1} - delta_{i, j+1}
    由于取周期性边界条件，索引可取模 N。
    """
    A = np.zeros((N, N))
    for i in range(N):
        A[i, i] = 2.0
        A[i, (i + 1) % N] = -1.0
        A[i, (i - 1) % N] = -1.0
    return A

def power_method(A, max_iter=2000, tol=1e-10):
    """
    使用幂次法求解矩阵的最大特征值及其对应的特征向量。
    """
    n = A.shape[0]
    # 从随机向量开始，使用固定种子以保证结果可复现
    np.random.seed(42)
    q = np.random.rand(n)
    q = q / np.linalg.norm(q)
    
    lambda_k = 0.0
    
    print(f"{'迭代次数':<10} {'特征值':<20} {'相邻两次迭代的误差':<20}")
    print("-" * 50)
    
    for k in range(1, max_iter + 1):
        z = A @ q
        
        # 计算瑞利商作为特征值估计
        lambda_next = np.dot(q, z)
        
        # 归一化得到新的q
        norm_z = np.linalg.norm(z)
        q_next = z / norm_z
        
        diff = np.abs(lambda_next - lambda_k)
        
        if k % 10 == 0 or k == 1:
            print(f"{k:<10} {lambda_next:<20.8f} {diff:<20.2e}")
            
        if diff < tol:
            print(f"\n在第 {k} 次迭代收敛")
            return lambda_next, q_next
            
        q = q_next
        lambda_k = lambda_next
        
    print("\n达到最大迭代次数。")
    return lambda_k, q

if __name__ == "__main__":
    N = 10
    print(f"正在构建 N = {N} 的矩阵...")
    A = construct_matrix(N)
    print("矩阵 A:")
    print(A)
    print("-" * 30)
    
    print("正在运行幂次法...")
    eigenvalue, eigenvector = power_method(A)
    
    print("-" * 30)
    print(f"最大特征值 (omega_max^2): {eigenvalue:.8f}")
    print("对应的特征向量:")
    print(eigenvector)
        
    # 验证
    print("-" * 30)
    print("验证:")
    print(f"A @ v = \n{A @ eigenvector}")
    print(f"lambda * v = \n{eigenvalue * eigenvector}")
    
    # N=10 时的理论最大特征值
    # omega_k^2 = 4 sin^2(pi k / N)
    # 最大值在 k = N/2 时取得 (如果 N 是偶数) -> 4 sin^2(pi/2) = 4
    print(f"理论最大特征值: 4.0")
