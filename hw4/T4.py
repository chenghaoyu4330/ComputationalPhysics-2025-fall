import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

def analytical_solution(x, y, t):
    """
    解析解: u(x,y,t) = sin(pi*x)*sin(2*pi*y)*cos(sqrt(5)*pi*t)
    """
    return np.sin(np.pi * x) * np.sin(2 * np.pi * y) * np.cos(np.sqrt(5) * np.pi * t)

def solve_wave_equation(Nx, Ny, T_max, dt, flag_use_lambda=False, lambda_param=1.0):
    """
    使用显式有限差分法求解二维波动方程。
    Nx, Ny: x和y方向的网格划分数
    T_max: 最大模拟时间
    lambda_param: 稳定性参数, 满足 dt <= 1/sqrt(lambda) * (1/dx^2 + 1/dy^2)^(-1/2)
    """
    L = 1.0
    dx = L / Nx
    dy = L / Ny

    x = np.linspace(0, L, Nx+1)
    y = np.linspace(0, L, Ny+1)
    X, Y = np.meshgrid(x, y, indexing='ij')

    if flag_use_lambda:
        dt = 1.0 / np.sqrt(lambda_param) * 1.0 / np.sqrt((1/dx**2) + (1/dy**2))
        print(f"\n取lambda={lambda_param}计算得到的时间步长 dt = {dt:.5e}")
    Nt = int(np.ceil(T_max / dt))

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
        u_np1[1:-1, 1:-1] = 2*u_n[1:-1, 1:-1] - u_nm1[1:-1, 1:-1] + rx2 * u_n_xx + ry2 * u_n_yy

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

def plot_comparison(X, Y, history, flag_stable_test=False, lambda_param=None):
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
    if flag_stable_test:
        plt.savefig(f'hw4/assets/wave_stability_test_lambda_{lambda_param}.png')
        print(f"稳定性测试对比图已保存到 hw4/assets/wave_stability_test_lambda_{lambda_param}.png")
    else:
        plt.savefig('hw4/assets/wave_comparison.png')
        print("对比图已保存到 hw4/assets/wave_comparison.png")

def run_stability_test(Nx, Ny, T_max, dt, lambda_values=None):
    """
    运行稳定性测试
    :param lambda_values: 要测试的 lambda 值列表。如果为 None，则默认为 [2.0, 0.8]
    """
    if lambda_values is None:
        lambda_values = [2.0, 0.8]

    print("\n--- 稳定性测试 ---")

    for lam in lambda_values:
        X, Y, history, frames, success, t_final = solve_wave_equation(Nx, Ny, T_max, dt, True, lam)
        if success:
            print(f"lambda = {lam}: 计算稳定，结束时刻 t={t_final:.2f}")
            plot_comparison(X, Y, history, True, lam)
        else:
            print(f"lambda = {lam}: 计算发散，失稳时刻 t={t_final:.2f}")
            plot_comparison(X, Y, history, True, lam)

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
    # (b) 数值解与解析解对比 (采用 CFL 极限 lambda=1.0)
    Nx, Ny = 100, 100
    T_max = 5.0
    dt = 0.001
    lambda_param = 1.0
    lambda_values = [3.0, 2.0, 1.01, 0.8]
    print(f"\n(b) 数值求解并与解析解对比 (lambda={lambda_param})...")
    X, Y, history, frames, success, _ = solve_wave_equation(Nx, Ny, T_max, dt, False, lambda_param)

    if success:
        plot_comparison(X, Y, history)

        # (c) Stability Test
        run_stability_test(Nx, Ny, T_max, dt, lambda_values)

        # (d) Animation
        create_animation(X, Y, frames)
    else:
        print("主模拟失败，可能由于数值发散或步长设置不当。")
