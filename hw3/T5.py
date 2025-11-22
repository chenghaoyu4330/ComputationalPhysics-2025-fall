import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import line_search

def himmelblau(p):
    x, y = p
    return (x**2 + y - 11)**2 + (x + y**2 - 7)**2

def grad_himmelblau(p):
    x, y = p
    df_dx = 4 * x * (x**2 + y - 11) + 2 * (x + y**2 - 7)
    df_dy = 2 * (x**2 + y - 11) + 4 * y * (x + y**2 - 7)
    return np.array([df_dx, df_dy])

def steepest_descent(func, grad, x0, tol=1e-5, max_iter=1000):
    x = np.array(x0, dtype=float)
    path = [x.copy()]
    
    for i in range(max_iter):
        g = grad(x)
        if np.linalg.norm(g) < tol:
            break
            
        p = -g # Direction is negative gradient
        
        # Line search to find alpha
        # func(x + alpha * p)
        res = line_search(func, grad, x, p)
        alpha = res[0]
        
        if alpha is None:
            # Fallback if line search fails
            alpha = 1e-3
            
        x = x + alpha * p
        path.append(x.copy())
        
    return x, np.array(path), i+1

def conjugate_gradient(func, grad, x0, tol=1e-5, max_iter=1000):
    x = np.array(x0, dtype=float)
    path = [x.copy()]
    
    g = grad(x)
    d = -g
    
    for i in range(max_iter):
        if np.linalg.norm(g) < tol:
            break
            
        # Line search
        res = line_search(func, grad, x, d)
        alpha = res[0]
        
        if alpha is None:
            alpha = 1e-3
            
        x_new = x + alpha * d
        path.append(x_new.copy())
        
        g_new = grad(x_new)
        
        # Fletcher-Reeves formula
        beta = np.dot(g_new, g_new) / np.dot(g, g)
        
        d = -g_new + beta * d
        
        g = g_new
        x = x_new
        
        # Reset direction every N iterations (optional but good practice for non-quadratic)
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
        # Steepest Descent
        sd_x, sd_path, sd_iter = steepest_descent(himmelblau, grad_himmelblau, p0)
        sd_val = himmelblau(sd_x)
        print(f"{'Steepest Descent':<20} | {str(p0):<15} | {str(np.round(sd_x, 4)):<25} | {sd_val:.2e}   | {sd_iter:<5}")
        
        # Conjugate Gradient
        cg_x, cg_path, cg_iter = conjugate_gradient(himmelblau, grad_himmelblau, p0)
        cg_val = himmelblau(cg_x)
        print(f"{'Conjugate Gradient':<20} | {str(p0):<15} | {str(np.round(cg_x, 4)):<25} | {cg_val:.2e}   | {cg_iter:<5}")
        print("-" * 85)
