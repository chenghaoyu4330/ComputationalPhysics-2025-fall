import numpy as np

def construct_matrix(N):
    """
    构建一维原子链周期性边界条件下的矩阵 A。
    -A_ij = delta_{i, j-1} + delta_{i, j+1} - 2 delta_{i, j}
    因此 A_ij = 2 delta_{i, j} - delta_{i, j-1} - delta_{i, j+1}
    对于周期性边界条件：索引取模 N。
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
    # 从随机向量开始（为了通用性使用复数，尽管 A 是实矩阵）
    # 使用固定种子以保证结果可复现
    np.random.seed(42)
    q = np.random.rand(n) + 0.0j
    q = q / np.linalg.norm(q)
    
    lambda_k = 0.0
    
    print(f"{'迭代次数':<10} {'特征值':<20} {'误差':<20}")
    print("-" * 50)
    
    for k in range(1, max_iter + 1):
        z = A @ q
        norm_z = np.linalg.norm(z)
        q_next = z / norm_z
        
        # 瑞利商: v^(k) = <q^(k)>^dagger A q^(k)
        # 由于 q_next 是新的 q，我们基于它计算 lambda
        # 实际上公式是 v^(k) = (q^(k))^H A q^(k)
        # q^(k) 是第 k 步归一化后的向量
        
        lambda_next = np.vdot(q_next, A @ q_next).real 
        # 对于厄米矩阵，特征值是实数
        
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
    # 如果虚部可以忽略，则打印实部
    if np.allclose(eigenvector.imag, 0):
        print(eigenvector.real)
    else:
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
