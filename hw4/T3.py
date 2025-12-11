import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# 参数设置
L = 2.0
N = 200
dx = L / N
delta = 0.022
delta_sq = delta**2
t_max = 6.0 # 模拟足够长的时间以观察孤立子
dt = 1e-4   # 稳定性要求
n_steps = int(t_max / dt)

# 网格设置
x = np.linspace(0, L, N, endpoint=False)

# 初始条件
u0 = np.cos(np.pi * x)

# 存储解的数组，Leapfrog格式只需要前一步、当前步和下一步
u_prev = np.copy(u0)
u_curr = np.copy(u0)
u_next = np.zeros_like(u0)

# 第一步使用前向欧拉法（或RK2）启动Leapfrog，方程: u_t = - (u u_x + delta^2 u_xxx)
# 使用空间导数的有限差分
def compute_rhs(u):
    # 使用 np.roll 处理周期性边界条件
    # u_x 近似为 (u_{j+1} - u_{j-1}) / 2dx
    # u u_x 采用 Zabusky-Kruskal 守恒格式: 1/3 (u_{j+1} + u_j + u_{j-1}) * (u_{j+1} - u_{j-1}) / 2dx
    
    u_jp1 = np.roll(u, -1)
    u_jm1 = np.roll(u, 1)
    u_jp2 = np.roll(u, -2)
    u_jm2 = np.roll(u, 2)
    
    # 非线性项: 1/3 * (u_{j+1} + u_j + u_{j-1}) * (u_{j+1} - u_{j-1}) / (2*dx)
    nonlinear = (u_jp1 + u + u_jm1) / 3.0 * (u_jp1 - u_jm1) / (2 * dx)
    
    # 色散项: delta^2 * (u_{j+2} - 2u_{j+1} + 2u_{j-1} - u_{j-2}) / (2*dx^3)
    dispersion = delta_sq * (u_jp2 - 2*u_jp1 + 2*u_jm1 - u_jm2) / (2 * dx**3)
    
    return -(nonlinear + dispersion)

# t=1 的欧拉步
rhs = compute_rhs(u_curr)
u_curr = u_prev + dt * rhs

# 存储动画帧
frames = []
save_interval = 100 # 每100步保存一次
frames.append(u_prev.copy())

print(f"开始模拟: N={N}, dt={dt}, 总步数={n_steps}")

for step in range(1, n_steps):
    # Leapfrog 步进：u^{n+1} = u^{n-1} + 2*dt * RHS(u^n)
    rhs = compute_rhs(u_curr)
    u_next = u_prev + 2 * dt * rhs
    
    # 更新数组
    u_prev[:] = u_curr[:]
    u_curr[:] = u_next[:]
    
    if step % save_interval == 0:
        frames.append(u_curr.copy())

print("模拟结束。")

# 动画生成
fig, ax = plt.subplots(figsize=(8, 5))
line, = ax.plot(x, frames[0], lw=2)
ax.set_xlim(0, L)
ax.set_ylim(-1.5, 3.0) # 孤立子会比初始的1.0更高
ax.set_xlabel('x')
ax.set_ylabel('u(x,t)')
ax.set_title('KdV Equation Solitons Evolution')
time_text = ax.text(0.02, 0.95, '', transform=ax.transAxes)

def update(frame_idx):
    u = frames[frame_idx]
    line.set_ydata(u)
    time_text.set_text(f't = {frame_idx * dt * save_interval:.2f}')
    return line, time_text

ani = animation.FuncAnimation(fig, update, frames=len(frames), interval=30, blit=True)

# 保存动画
output_file = 'hw4/assets/kdv_solitons.gif'
print(f"正在保存动画至 {output_file}...")
try:
    ani.save(output_file, writer='pillow', fps=30)
    print("动画保存成功。")
except Exception as e:
    print(f"保存动画出错: {e}")
    # 尝试保存为mp4
    try:
        output_file_mp4 = 'hw4/assets/kdv_solitons.mp4'
        ani.save(output_file_mp4, writer='ffmpeg', fps=30)
        print(f"动画已保存为 {output_file_mp4}")
    except Exception as e2:
        print(f"无法保存为mp4: {e2}")

plt.close()
