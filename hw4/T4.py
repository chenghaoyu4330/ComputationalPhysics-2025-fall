import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

def analytical_solution(x, y, t):
    """
    解析解: u(x,y,t) = sin(pi*x)*sin(2*pi*y)*cos(sqrt(5)*pi*t)
    """
    return np.sin(np.pi * x) * np.sin(2 * np.pi * y) * np.cos(np.sqrt(5) * np.pi * t)

def solve_wave_equation(Nx, Ny, T_max, lambda_param=1.0):
    """
    使用显式有限差分法求解二维波动方程。
    Nx, Ny: x和y方向的网格划分数
    T_max: 最大模拟时间
    lambda_param: 稳定性参数, 满足 dt <= 1/sqrt(lambda) * (1/dx^2 + 1/dy^2)^(-1/2)
    """
    L = 1.0
    dx = L / Nx
    dy = L / Ny

    # 根据题目公式计算稳定性上界
    dt_limit = 1.0 / np.sqrt(1.0/dx**2 + 1.0/dy**2)
    dt = (1.0 / np.sqrt(lambda_param)) * dt_limit

    Nt = int(np.ceil(T_max / dt))

    x = np.linspace(0, L, Nx+1)
    y = np.linspace(0, L, Ny+1)
    X, Y = np.meshgrid(x, y, indexing='ij')

    # 存储 u^{n-1}, u^n 的数组（Leapfrog 形式）
    u_prev = np.zeros((Nx+1, Ny+1))

    # 初始条件: u(x,y,0)
    u_prev[:, :] = np.sin(np.pi * X) * np.sin(2 * np.pi * Y)

    # 第一步 (n=0 到 n=1) 使用 u_t(0) = 0
    rx2 = (dt / dx)**2
    ry2 = (dt / dy)**2

    u_nm1 = u_prev.copy() # u^{n-1}

    # 计算 u^1
    u_nm1_xx = u_nm1[2:, 1:-1] - 2*u_nm1[1:-1, 1:-1] + u_nm1[:-2, 1:-1]
    u_nm1_yy = u_nm1[1:-1, 2:] - 2*u_nm1[1:-1, 1:-1] + u_nm1[1:-1, :-2]

    u_n = u_nm1.copy()
    u_n[1:-1, 1:-1] = u_nm1[1:-1, 1:-1] + 0.5 * (rx2 * u_nm1_xx + ry2 * u_nm1_yy)

    # 记录需要保存的时刻 (使用最接近的数值时刻)
    targets = {t: {'diff': float('inf'), 'u': None} for t in [0.0, 1.0, 2.0]}
    targets[0.0] = {'diff': 0.0, 'u': u_nm1.copy()}

    current_time = dt
    frames = [u_nm1.copy()]

    for step in range(1, Nt):
        u_n_xx = u_n[2:, 1:-1] - 2*u_n[1:-1, 1:-1] + u_n[:-2, 1:-1]
        u_n_yy = u_n[1:-1, 2:] - 2*u_n[1:-1, 1:-1] + u_n[1:-1, :-2]

        u_np1 = np.zeros_like(u_n)
        u_np1[1:-1, 1:-1] = 2*u_n[1:-1, 1:-1] - u_nm1[1:-1, 1:-1] + \
                            rx2 * u_n_xx + ry2 * u_n_yy

        # 边界条件保持为0
        u_np1[0, :] = 0.0
        u_np1[-1, :] = 0.0
        u_np1[:, 0] = 0.0
        u_np1[:, -1] = 0.0

        # 检查是否保存历史记录（选择距离目标时刻最近的一帧）
        for target_t in targets:
            diff = abs(current_time - target_t)
            if diff < targets[target_t]['diff']:
                targets[target_t] = {'diff': diff, 'u': u_n.copy(), 't': current_time}

        if step % 5 == 0:
            frames.append(u_n.copy())

        u_nm1 = u_n
        u_n = u_np1
        current_time += dt

        # 检查不稳定性
        if np.max(np.abs(u_n)) > 1e5: # 提前捕捉发散
            history = [{'t': targets[t].get('t', t), 'u': targets[t]['u']} for t in targets if targets[t]['u'] is not None]
            return X, Y, history, frames, False, current_time

    history = [{'t': targets[t].get('t', t), 'u': targets[t]['u']} for t in targets if targets[t]['u'] is not None]
    history = sorted(history, key=lambda item: item['t'])
    return X, Y, history, frames, True, current_time

def plot_comparison(X, Y, history):
    """
    绘制 t=0,1,2 时刻数值解与解析解的对比。
    """
    if not history:
        print("无可用的历史数据，无法绘图。")
        return

    fig = plt.figure(figsize=(15, 10))

    for i, h in enumerate(history):
        t = h['t']
        u_num = h['u']
        u_ana = analytical_solution(X, Y, t)

        ax1 = fig.add_subplot(3, 3, i*3 + 1, projection='3d')
        ax1.plot_surface(X, Y, u_num, cmap='viridis')
        ax1.set_title(f'Numerical t={t:.2f}')
        ax1.set_zlim(-1, 1)

        ax2 = fig.add_subplot(3, 3, i*3 + 2, projection='3d')
        ax2.plot_surface(X, Y, u_ana, cmap='viridis')
        ax2.set_title(f'Analytical t={t:.2f}')
        ax2.set_zlim(-1, 1)

        ax3 = fig.add_subplot(3, 3, i*3 + 3)
        im = ax3.pcolormesh(X, Y, u_num - u_ana, cmap='coolwarm', shading='auto')
        plt.colorbar(im, ax=ax3)
        ax3.set_title(f'Error t={t:.2f}')

    plt.tight_layout()
    plt.savefig('hw4/assets/wave_comparison.png')
    print("对比图已保存到 hw4/assets/wave_comparison.png")

def run_stability_test():
    print("\n--- 稳定性测试 ---")
    Nx, Ny = 50, 50

    # 情形1: lambda = 2.0 (更严格的稳定步长)
    print("lambda = 2.0: 预期稳定，正在计算...")
    _, _, _, _, success, t_final = solve_wave_equation(Nx, Ny, T_max=2.1, lambda_param=2.0)
    if success:
        print(f"-> 计算稳定，结束时刻 t={t_final:.2f}")
    else:
        print(f"-> 计算发散，失稳时刻 t={t_final:.2f}")

    # 情形2: lambda = 0.8 (步长放宽，预期可能失稳)
    print("lambda = 0.8: 预期可能失稳，正在计算...")
    _, _, _, _, success, t_final = solve_wave_equation(Nx, Ny, T_max=2.1, lambda_param=0.8)
    if success:
        print(f"-> 出乎意料地稳定，结束时刻 t={t_final:.2f}")
    else:
        print(f"-> 如预期失稳，发散时刻 t={t_final:.2f}")

def create_animation(X, Y, frames):
    print("\n正在生成动画...")
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection='3d')

    # Initial plot
    surf = ax.plot_surface(X, Y, frames[0], cmap='viridis', vmin=-1, vmax=1)
    ax.set_zlim(-1, 1)
    ax.set_title('2D Wave Equation')

    def update(frame):
        ax.clear()
        ax.set_zlim(-1, 1)
        ax.set_title('2D Wave Equation')
        surf = ax.plot_surface(X, Y, frame, cmap='viridis', vmin=-1, vmax=1)
        return surf,

    ani = animation.FuncAnimation(fig, update, frames=frames[::2], interval=50)

    try:
        ani.save('hw4/assets/wave_2d.gif', writer='pillow', fps=20)
        print("动画已保存到 hw4/assets/wave_2d.gif")
    except Exception as e:
        print(f"保存动画时出错: {e}")

if __name__ == "__main__":
    # (a) 分离变量解析解
    print("(a) 解析解: u(x,y,t) = sin(pi*x)*sin(2*pi*y)*cos(sqrt(5)*pi*t)")

    # (b) 数值解与解析解对比 (采用 CFL 极限 lambda=1.0)
    print("\n(b) 数值求解并与解析解对比 (lambda=1.0)...")
    Nx, Ny = 50, 50
    X, Y, history, frames, success, _ = solve_wave_equation(Nx, Ny, T_max=2.1, lambda_param=1.0)

    if success:
        plot_comparison(X, Y, history)

        # (c) Stability Test
        run_stability_test()

        # (d) Animation
        create_animation(X, Y, frames)
    else:
        print("主模拟失败，可能由于数值发散或步长设置不当。")
