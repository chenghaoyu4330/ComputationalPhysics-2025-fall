import numpy as np
import pandas as pd

# Set precision for display
pd.set_option('display.precision', 6)

def f(x):
    return x - 2 * np.sin(x)

def df(x):
    return 1 - 2 * np.cos(x)

def g(x):
    return (x - 2 * np.sin(x))**2

def dg(x):
    return 2 * (x - 2 * np.sin(x)) * (1 - 2 * np.cos(x))

def bisection(func, a, b, tol=1e-5, max_iter=100):
    results = []
    results.append({'iter': 0, 'a': a, 'b': b, 'x': np.nan, 'err': np.nan})
    
    if func(a) * func(b) > 0:
        return results, False

    for i in range(1, max_iter + 1):
        c = (a + b) / 2
        err = (b - a) / 2
        
        # Store the state BEFORE update (the interval used for this step)
        results.append({'iter': i, 'a': a, 'b': b, 'x': c, 'err': err})
        
        # Update interval
        if func(c) == 0:
            a = c
            b = c
        elif func(a) * func(c) < 0:
            b = c
        else:
            a = c
        
        if err < tol:
            break
            
    return results, True

def newton(func, dfunc, x0, tol=1e-5, max_iter=100):
    results = []
    x = x0
    results.append({'iter': 0, 'x': x, 'err': np.nan})
    
    for i in range(1, max_iter + 1):
        fx = func(x)
        dfx = dfunc(x)
        if dfx == 0:
            break
        
        x_new = x - fx / dfx
        err = abs(x_new - x)
        x = x_new
        results.append({'iter': i, 'x': x, 'err': err})
        
        if err < tol:
            break
            
    return results

def secant(func, x0, x1, tol=1e-5, max_iter=100):
    results = []

    results.append({'iter': 0, 'x': x0, 'err': np.nan})
    results.append({'iter': 1, 'x': x1, 'err': abs(x1 - x0)}) 
    
    x_prev = x0
    x_curr = x1
    
    for i in range(2, max_iter + 1):
        fx_curr = func(x_curr)
        fx_prev = func(x_prev)
        
        if fx_curr == fx_prev:
            break
            
        x_new = x_curr - fx_curr * (x_curr - x_prev) / (fx_curr - fx_prev)
        err = abs(x_new - x_curr)
        
        x_prev = x_curr
        x_curr = x_new
        
        results.append({'iter': i, 'x': x_curr, 'err': err})
        
        if err < tol:
            break
            
    return results


# 题目
if __name__ == "__main__":
    print("第一题")
    tol = 1e-5
    
    # Bisection
    bis_res, bis_valid = bisection(f, 1.5, 2.0, tol)
    
    # Newton
    newton_res = newton(f, df, 1.5, tol)
    
    # Secant
    x0_sec = 1.5
    if len(newton_res) > 1:
        x1_sec = newton_res[1]['x']
    else:
        x1_sec = 1.5
        
    secant_res = secant(f, x0_sec, x1_sec, tol)
    
    # 输出表格
    max_len = max(len(bis_res), len(newton_res), len(secant_res))
    
    print(f"{'Iter':<5} | {'Bisection [a, b]':<25} | {'Newton x':<12} | {'Secant x':<12}")
    print("-" * 65)
    
    for i in range(max_len):
        # Bisection
        if i < len(bis_res):
            r = bis_res[i]
            bis_str = f"[{r['a']:.5f}, {r['b']:.5f}]"
        else:
            bis_str = ""
            
        # Newton
        if i < len(newton_res):
            newt_str = f"{newton_res[i]['x']:.6f}"
        else:
            newt_str = ""
            
        # Secant
        if i < len(secant_res):
            sec_str = f"{secant_res[i]['x']:.6f}"
        else:
            sec_str = ""
            
        print(f"{i:<5} | {bis_str:<25} | {newt_str:<12} | {sec_str:<12}")

    # 收敛性
    # Bisection
    for r in bis_res:
        if r['iter'] > 0 and r['err'] < tol:
            print(f"二分法在第 {r['iter']} 次迭代时收敛 (误差 < {tol})")
            break
            
    # Newton
    for r in newton_res:
        if r['iter'] > 0 and r['err'] < tol:
            print(f"牛顿法在第 {r['iter']} 次迭代时收敛 (误差 < {tol})")
            break
            
    # Secant
    for r in secant_res:
        if r['iter'] > 1 and r['err'] < tol: # iter 0 and 1 are setup
            print(f"割线法在第 {r['iter']} 次迭代时收敛 (误差 < {tol})")
            break



    print("\n第二题")
    tol = 1e-5
    
    # 检查二分法
    print("检查二分法对于 g(x) = (x - 2sin(x))^2 在区间 [1.5, 2] 上的适用性:")
    ga = g(1.5)
    gb = g(2.0)
    print(f"g(1.5) = {ga:.6f}, g(2.0) = {gb:.6f}")
    if ga * gb > 0:
        print("二分法不适用 (g(a) 和 g(b) 同号)。")
    else:
        print("二分法适用。")
        
    # Newton
    print("\n牛顿法求解 g(x):")
    newton_res = newton(g, dg, 1.5, tol)
    
    # Secant
    print("割线法求解 g(x):")
    x0_sec = 1.5
    if len(newton_res) > 1:
        x1_sec = newton_res[1]['x']
    else:
        x1_sec = 1.5
    secant_res = secant(g, x0_sec, x1_sec, tol)
    
    # 输出表格
    max_len = max(len(newton_res), len(secant_res))
    print(f"{'Iter':<5} | {'Newton x':<15} | {'Secant x':<15}")
    print("-" * 45)
    
    for i in range(max_len):
        if i < len(newton_res):
            newt_str = f"{newton_res[i]['x']:.6f}"
        else:
            newt_str = ""
            
        if i < len(secant_res):
            sec_str = f"{secant_res[i]['x']:.6f}"
        else:
            sec_str = ""
            
        print(f"{i:<5} | {newt_str:<15} | {sec_str:<15}")

    # 收敛性
    for r in newton_res:
        if r['iter'] > 0 and r['err'] < tol:
            print(f"牛顿法在第 {r['iter']} 次迭代时收敛")
            break
    for r in secant_res:
        if r['iter'] > 1 and r['err'] < tol:
            print(f"割线法在第 {r['iter']} 次迭代时收敛")
            break
