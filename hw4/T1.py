import numpy as np

def get_matrix():
    return np.array([
        [1, -1, 0, 0],
        [-1, 2, -1, 0],
        [0, -1, 3, -1],
        [0, 0, -1, 4]
    ], dtype=float)

def householder_qr(A: np.ndarray):
    """
    使用Householder变换进行QR分解
    :param A: 输入矩阵(m*n, m≥n, 列满秩)
    :return: 正交矩阵Q(m*m), 上三角矩阵R(m*n)
    """
    # 复制矩阵并转换为浮点型，避免原矩阵被修改
    m, n = A.shape
    R = A.copy().astype(float)
    # 初始化Q为m阶单位矩阵
    Q = np.eye(m, dtype=float)

    # 遍历每一列
    for k in range(n):
        # 提取当前列的子向量（从第k行到最后一行）
        x = R[k:, k].copy()
        # 计算x的2-范数
        norm_x = np.linalg.norm(x, ord=2)
        # 如果范数为0，说明该列已经是0，跳过
        if norm_x == 0:
            continue

        # 构造单位向量e1（长度和x一致）
        e1 = np.zeros_like(x)
        e1[0] = 1.0
        # 选择符号，避免数值抵消（和x[0]同号）
        sign = np.sign(x[0]) if x[0] != 0 else 1.0
        # 构造Householder向量v并归一化
        v = x + sign * norm_x * e1
        v /= np.linalg.norm(v, ord=2)

        # 构造m-k阶的Householder矩阵H
        h_size = len(v)
        H = np.eye(h_size, dtype=float) - 2.0 * np.outer(v, v)

        # 步骤1：将Householder变换应用到R的子矩阵（左乘H）
        R[k:, k:] = H @ R[k:, k:]

        # 步骤2：构造m阶的完整Householder矩阵（嵌入到单位矩阵的右下角）
        H_full = np.eye(m, dtype=float)
        H_full[k:, k:] = H

        # 步骤3：将完整的Householder矩阵应用到Q（右乘H_full）
        Q = Q @ H_full

    # 对于n×n矩阵，R的下三角部分（除了上三角）可以置0，避免数值误差导致的小值
    for i in range(m):
        for j in range(n):
            if i > j:
                R[i, j] = 0.0

    return Q, R

def qr_algorithm(A, max_iter=50, tol=1e-6):
    n = A.shape[0]
    Ak = A.copy()
    print(f"初始矩阵:\n{Ak}")
    
    for k in range(1, max_iter + 1):
        Q, R = householder_qr(Ak)
        Ak = R @ Q
        
        if k % 5 == 0:
            print(f"\n第 {k} 次迭代:")
            print(Ak)
        
        # 检查收敛性，判断非对角元素的绝对值和是否小于容差
        off_diagonal = np.sum(np.abs(Ak)) - np.sum(np.abs(np.diag(Ak)))
        if off_diagonal < tol:
            print(f"\n{k} 次迭代后收敛。")
            break
            
    print("\nQR算法特征值:", np.diag(Ak))
    return np.diag(Ak)

def jacobi_algorithm(A, max_iter=1000, tol=1e-8):
    n = A.shape[0]
    Ak = A.copy()
    
    for k in range(max_iter):
        # 找到最大的非对角元素
        max_val = 0.0
        p, q = -1, -1
        for i in range(n):
            for j in range(i + 1, n):
                if abs(Ak[i, j]) > max_val:
                    max_val = abs(Ak[i, j])
                    p, q = i, j
        
        if max_val < tol:
            print(f"\n{k}次迭代后收敛。")
            break
            
        # 计算旋转参数
        if Ak[p, q] == 0:
            c = 1.0
            s = 0.0
        else:
            eta = (Ak[q, q] - Ak[p, p]) / (2 * Ak[p, q])
            if eta >= 0:
                t = 1.0 / (eta + np.sqrt(1 + eta**2))
            else:
                t = -1.0 / (-eta + np.sqrt(1 + eta**2))
            
            c = 1.0 / np.sqrt(1 + t**2)
            s = t * c
        
        # 构造Jacobi旋转矩阵
        J = np.eye(n)
        J[p, p] = c
        J[q, q] = c
        J[p, q] = s
        J[q, p] = -s
        
        Ak = J.T @ Ak @ J
        
    print("\nJacobi算法特征值:", np.diag(Ak))
    return np.diag(Ak)

def sturm_sequence_count(d, e, lam):
    """
    计算Sturm序列在λ处的符号变化次数
    :param d: 矩阵的对角元素
    :param e: 矩阵的副对角元素
    """
    n = len(d)
    
    # 初始化Sturm序列，P0(λ) = 1.0
    seq = [1.0]

    # 计算P1(λ) = d[0] - λ
    val = d[0] - lam

    # 用一个非常小的数替换零值，方便符号变化的计算
    if val == 0:
        val = 1e-15
    seq.append(val)
    
    # 递推计算P_k(λ)： P_k(λ) = (d[k-1] - λ) * P_{k-1}(λ) - (e[k-2])^2 * P_{k-2}(λ)
    for k in range(2, n + 1):
        val = (d[k-1] - lam) * seq[-1] - (e[k-2]**2) * seq[-2]
        # 替换零值
        if val == 0:
            val = 1e-15
        seq.append(val)
        
    # 计算符号变化次数
    changes = 0
    for i in range(len(seq) - 1):
        if (seq[i] > 0 and seq[i+1] < 0) or (seq[i] < 0 and seq[i+1] > 0):
            changes += 1
    return changes

def sturm_bisection(A, tol=1e-6):
    # 提取对角线和副对角线元素
    d = np.diag(A)
    e = np.diag(A, k=1) 
    n = len(d)
    
    # 利用Gerschgorin圆盘定理计算特征值的初始搜索区间
    max_val = -np.inf
    min_val = np.inf
    for i in range(n):
        row_sum = np.sum(np.abs(A[i, :])) - abs(A[i, i])
        if (row_sum + A[i, i]) > max_val:
            max_val = row_sum + A[i, i]
        if (A[i, i] - row_sum) < min_val:
            min_val = A[i, i] - row_sum
    
    # 避免边界问题，稍微扩大区间
    low = min_val - 1.0
    high = max_val + 1.0
    
    print(f"\n搜索区间: [{low}, {high}]")
    
    found_evals = []
    # 求解第1到第n个特征值（对应k=1到n）
    for k in range(1, n + 1):
        a, b = low, high
        
        for _ in range(100): # 最多迭代100次
            mid = (a + b) / 2
            count = sturm_sequence_count(d, e, mid)
            if count < k:
                # 符号变化次数小于k，说明第k个特征值在mid右侧，更新左边界
                a = mid
            else:
                # 符号变化次数≥k，说明第k个特征值在mid左侧或等于mid，更新右边界
                b = mid
            
            if abs(b - a) < tol:
                break
        
        found_evals.append((a + b) / 2)
        
    print("\nSturm序列+二分法特征值:", np.array(found_evals))
    return np.array(found_evals)

if __name__ == "__main__":
    T = get_matrix()
    print("初始矩阵T:")
    print(T)
    print("-" * 30)
    
    print("1. QR算法")
    qr_algorithm(T)
    print("-" * 30)
    
    print("2. Jacobi算法")
    jacobi_algorithm(T)
    print("-" * 30)
    
    print("3. Sturm序列+二分法")
    sturm_bisection(T)
