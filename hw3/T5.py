import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import line_search

def f(p):
    x, y = p
    return (x**2 + y - 11)**2 + (x + y**2 - 7)**2

def grad_f(p):
    x, y = p
    df_dx = 4 * x * (x**2 + y - 11) + 2 * (x + y**2 - 7)
    df_dy = 2 * (x**2 + y - 11) + 4 * y * (x + y**2 - 7)
    return np.array([df_dx, df_dy])

# 最速下降法
def steepest_descent(func, grad, x0, tol=1e-5, max_iter=1000):
    x = np.array(x0, dtype=float)  # 初始点
    path = [x.copy()]  # 记录路径
    
    for i in range(max_iter):
        g = grad(x)  # 计算梯度
        # 检查收敛，如果梯度范数小于容差则停止
        if np.linalg.norm(g) < tol:
            break
            
        p = -g  # 最速下降方向，负梯度方向
        
        # 线搜索确定步长alpha
        res = line_search(func, grad, x, p)
        alpha = res[0]
        
        if alpha is None:
            # 如果线搜索失败，使用默认步长
            alpha = 1e-3
            
        x = x + alpha * p  # 更新位置
        path.append(x.copy())
        
    return x, np.array(path), i+1

# 共轭梯度法
def conjugate_gradient(func, grad, x0, tol=1e-5, max_iter=1000):
    x = np.array(x0, dtype=float)
    path = [x.copy()]
    
    g = grad(x)
    d = -g  # 初始搜索方向为负梯度方向
    
    for i in range(max_iter):
        # 检查收敛，如果梯度范数小于容差则停止
        if np.linalg.norm(g) < tol:
            break
            
        # 线搜索确定步长alpha
        res = line_search(func, grad, x, d)
        alpha = res[0]
        
        if alpha is None:
            alpha = 1e-3  # 如果线搜索失败，使用默认步长
            
        x_new = x + alpha * d  # 更新位置
        path.append(x_new.copy())
        
        g_new = grad(x_new)  # 计算新的梯度
        
        # 计算beta，使用Fletcher-Reeves公式
        beta = np.dot(g_new, g_new) / np.dot(g, g)
        
        d = -g_new + beta * d  # 更新搜索方向
        
        g = g_new  # 更新梯度
        x = x_new  # 更新位置
        
        # 每两次迭代重置搜索方向为负梯度方向
        if (i + 1) % 2 == 0:
             d = -g
             
    return x, np.array(path), i+1


# 题目
if __name__ == "__main__":
    start_points = [
        [0, 0],
        [-1, -1],
        [4, 4],
        [-4, 4]
    ]
    
    print(f"{'Method':<20} | {'Start':<15} | {'End Point':<25} | {'Value':<10} | {'Iter':<5}")
    print("-" * 85)
    
    for p0 in start_points:
        # 最速下降法
        sd_x, sd_path, sd_iter = steepest_descent(f, grad_f, p0)
        sd_val = f(sd_x)
        print(f"{'Steepest Descent':<20} | {str(p0):<15} | {str(np.round(sd_x, 4)):<25} | {sd_val:.2e}   | {sd_iter:<5}")
        
        # 共轭梯度法
        cg_x, cg_path, cg_iter = conjugate_gradient(f, grad_f, p0)
        cg_val = f(cg_x)
        print(f"{'Conjugate Gradient':<20} | {str(p0):<15} | {str(np.round(cg_x, 4)):<25} | {cg_val:.2e}   | {cg_iter:<5}")
        print("-" * 85)

    # 可视化路径
    x = np.linspace(-6, 6, 400)
    y = np.linspace(-6, 6, 400)
    X, Y = np.meshgrid(x, y)
    Z = f([X, Y])
    fig, axs = plt.subplots(2, 2, figsize=(12, 10))
    axs = axs.flatten()
    for idx, p0 in enumerate(start_points):
        ax = axs[idx]
        ax.contour(X, Y, Z, levels=np.logspace(0, 5, 35), cmap='jet')
        
        # 最速下降法路径
        sd_x, sd_path, _ = steepest_descent(f, grad_f, p0)
        ax.plot(sd_path[:, 0], sd_path[:, 1], 'ro-', label='Steepest Descent', markersize=4)
        
        # 共轭梯度法路径
        cg_x, cg_path, _ = conjugate_gradient(f, grad_f, p0)
        ax.plot(cg_path[:, 0], cg_path[:, 1], 'go-', label='Conjugate Gradient', markersize=4)
        
        # 标出起点和终点
        ax.plot(p0[0], p0[1], 's', color='orange', label='Start', markersize=8)
        ax.plot(sd_x[0], sd_x[1], 'r*', label='SD End', markersize=10)
        ax.plot(cg_x[0], cg_x[1], 'g*', label='CG End', markersize=10)
        
        ax.set_title(f'Start Point: {p0}')
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.legend()
        ax.grid()
    plt.tight_layout()
    plt.show()
    