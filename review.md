# 计算物理学 2025Fall 复习笔记

2025.12.29 by Ritmo

授课教师：李强

## 0. 考后更新

12.30 考后更新：今年期末考试黑化了。。。七道题考了三道往年题（QR分解，二阶RK方法，常微分方程差分法及其稳定性分析），剩下四道都是新题：一道题考了高斯积分法，在这份资料的3.2节，考的内容基本上就是3.2.2的内容；一道考数值积分精度和误差，没有覆盖到；一道考利用插值法进行证明：给定 $n+1$ 个点 $x_0, x_1, \ldots, x_n$ 以及函数 $f(x)=cx^n+x$，证明下面的式子成立：

1. 
   $$
   \sum_{k=0}^n \frac{x_k^n}{\prod_{j=0,j\ne k}^n (x_k - x_j)} = 1,
   $$

2. 
   $$
   f[x_0, x_1, \ldots, x_{n-1}] = c\sum_{k=0}^{n-1} x_k,
   $$
   其中， $f[x_0, x_1, \ldots, x_{i}]$ 是对 $f(x)$ 在点 $x_0, x_1, \ldots, x_n$ 应用 Newton 插值法得到的 $i$ 阶差商。

最后一题考抽样，n维超球表面进行均匀抽样，从未见过，难。。

## 1. 线性方程组

### 1.1 高斯消元法

$$
\begin{aligned}
(U, \vec{c}) = G^{(n-1)} P^{(n-1)} \cdots G^{(1)} P^{(1)} (A, \vec{b})
\end{aligned}
$$

其中，$P^{(k)}$ 是第 $k$ 步的行交换矩阵，$G^{(k)}$ 是第 $k$ 步的消元矩阵。得到上三角矩阵 $U$ 和变换后的向量 $\vec{c}$ 后，通过回代求解方程组 $U \vec{x} = \vec{c}$。

> e.g. $A = \begin{pmatrix} 1 & 3 & 1 \\ 3 & 4 & 2 \\ -1 & -5 & 4 \end{pmatrix}, \vec{b} = \begin{pmatrix} 2 \\ 9 \\ 10 \end{pmatrix}$，求解 $A \vec{x} = \vec{b}$。
>
> 步骤如下：
> $$
> \begin{aligned}
> P^{(1)} &= \begin{pmatrix} 0 & 1 & 0 \\ 1 & 0 & 0 \\ 0 & 0 & 1 \end{pmatrix}, \quad P^{(1)}(A,\vec{b}) = \begin{pmatrix} 3 & 4 & 2 & 9 \\ 1 & 3 & 1 & 2 \\ -1 & -5 & 4 & 10 \end{pmatrix} \\
> G^{(1)} &= \begin{pmatrix} 1 & 0 & 0 \\ -\frac{1}{3} & 1 & 0 \\ \frac{1}{3} & 0 & 1 \end{pmatrix} , \quad (A_1,\vec{b_1})=G^{(1)} P^{(1)}(A,\vec{b}) = \begin{pmatrix} 3 & 4 & 2 & 9 \\ 0 & \frac{5}{3} & \frac{1}{3} & -1 \\ 0 & -\frac{11}{3} & \frac{14}{3} & 13 \end{pmatrix} \\
> P^{(2)} &= \begin{pmatrix} 1 & 0 & 0 \\ 0 & 0 & 1 \\ 0 & 1 & 0 \end{pmatrix}, \quad P^{(2)}(A_1,\vec{b_1}) = \begin{pmatrix} 3 & 4 & 2 & 9 \\ 0 & -\frac{11}{3} & \frac{14}{3} & 13 \\ 0 & \frac{5}{3} & \frac{1}{3} & -1 \end{pmatrix} \\
> G^{(2)} &= \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & \frac{5}{11} & 1 \end{pmatrix}, \quad (U,\vec{c})=G^{(2)} P^{(2)}(A_1,\vec{b_1}) = \begin{pmatrix} 3 & 4 & 2 & 9 \\ 0 & -\frac{11}{3} & \frac{14}{3} & 13 \\ 0 & 0 & \frac{27}{11} & \frac{54}{11} \end{pmatrix}
> \end{aligned}
> $$
> 回代求解可得 $\vec{x} = \begin{pmatrix} 3 \\ -1 \\ 2 \end{pmatrix}$。

### 1.2 LU 分解

在高斯消元法的基础上，令$P = P^{(n-1)} P^{(n-2)} \cdots P^{(1)}$，$M = G^{(n-1)} P^{(n-1)} \cdots G^{(1)} P^{(1)}$，则有

$$
\begin{aligned}
PA = LU, \quad \text{其中} L = PM^{-1}， U = M A
\end{aligned}
$$

> e.g. 继续上例，求解 $PA = LU$ 分解。
> 步骤如下：
> $$
> \begin{aligned}
> P &= P^{(2)} P^{(1)} = \begin{pmatrix} 0 & 1 & 0 \\ 0 & 0 & 1 \\ 1 & 0 & 0 \end{pmatrix} \\
> L &= PM^{-1} = \begin{pmatrix} 1 & 0 & 0 \\ -\frac{1}{3} & 1 & 0 \\ \frac{1}{3} & -\frac{5}{11} & 1 \end{pmatrix} \\
> U &= MA = \begin{pmatrix} 3 & 4 & 2 \\ 0 & -\frac{11}{3} & \frac{14}{3} \\ 0 & 0 & \frac{27}{11} \end{pmatrix}
> \end{aligned}
> $$
> 验证 $PA = LU$：
> $$
> \begin{aligned}
> PA &= \begin{pmatrix} 0 & 1 & 0 \\ 0 & 0 & 1 \\ 1 & 0 & 0 \end{pmatrix} \begin{pmatrix} 1 & 3 & 1 \\ 3 & 4 & 2 \\ -1 & -5 & 4 \end{pmatrix} = \begin{pmatrix} 3 & 4 & 2 \\ -1 & -5 & 4 \\ 1 & 3 & 1 \end{pmatrix} \\
> LU &= \begin{pmatrix} 1 & 0 & 0 \\ -\frac{1}{3} & 1 & 0 \\ \frac{1}{3} & -\frac{5}{11} & 1 \end{pmatrix} \begin{pmatrix} 3 & 4 & 2 \\ 0 & -\frac{11}{3} & \frac{14}{3} \\ 0 & 0 & \frac{27}{11} \end{pmatrix} = \begin{pmatrix} 3 & 4 & 2 \\ -1 & -5 & 4 \\ 1 & 3 & 1 \end{pmatrix}
> \end{aligned}
> $$

### 1.3 Cholesky 分解

对于正定厄米矩阵 $A$，存在唯一的下三角矩阵 $L$，使得

$$
\begin{aligned}
A = LL^{\dagger}
\end{aligned}
$$

记 $A_n = \begin{pmatrix} a_{11} & a_{12} & \cdots & a_{1n} \\ a_{21} & a_{22} & \cdots & a_{2n} \\ \vdots & \vdots & \ddots & \vdots \\ a_{n1} & a_{n2} & \cdots & a_{nn} \end{pmatrix} = \begin{pmatrix} A_{n-1} & \vec{v_n} \\ \vec{v_n}^{\dagger} & a_{nn} \end{pmatrix} = L_n L_n^{\dagger} = \begin{pmatrix} L_{n-1} & 0 \\ \vec{l_n}^{\dagger} & l_{nn} \end{pmatrix} \begin{pmatrix} L_{n-1}^{\dagger} & \vec{l_n} \\ 0 & l_{nn} \end{pmatrix}$，则有

$$
\begin{aligned}
\vec{v_n} &= L_{n-1} \vec{l_n} \\
a_{nn} &= \vec{l_n}^{\dagger} \vec{l_n} + l_{nn}^2
\end{aligned}
$$

由此可递推地计算出 $L$ 的各个元素：

$$
\begin{aligned}
l_{ij} &= \frac{1}{l_{jj}} \left( a_{ij} - \sum_{k=1}^{j-1} l_{ik} l_{jk} \right), \quad j = 1, 2, \ldots, i-1 \\
l_{ii} &= \sqrt{a_{ii} - \sum_{k=1}^{i-1} l_{ik}^2}, \quad i = 1, 2, \ldots, n
\end{aligned}
$$

> e.g. $A = \begin{pmatrix} 1 & 2 & 3 \\ 2 & 5 & 8 \\ 3 & 8 & 14 \end{pmatrix}$，求$A$的 Cholesky 分解。
>
> 步骤如下：
>
> 1. $i=1$:
>    $$
>    \begin{aligned}
>    l_{11} &= \sqrt{a_{11}} = \sqrt{1} = 1
>    \end{aligned}
>    $$
> 2. $i=2$:
>    $$
>    \begin{aligned}
>    l_{21} &= \frac{1}{l_{11}} (a_{21}) = \frac{1}{1} (2) = 2 \\
>    l_{22} &= \sqrt{a_{22} - l_{21}^2} = \sqrt{5 - 2^2} = 1
>    \end{aligned}
>    $$
> 3. $i=3$:
>    $$
>    \begin{aligned}
>    l_{31} &= \frac{1}{l_{11}} (a_{31}) = \frac{1}{1} (3) = 3 \\
>    l_{32} &= \frac{1}{l_{22}} (a_{32} - l_{31} l_{21}) = \frac{1}{1} (8 - 3 \cdot 2) = 2 \\
>    l_{33} &= \sqrt{a_{33} - l_{31}^2 - l_{32}^2} = \sqrt{14 - 3^2 - 2^2} = 1
>    \end{aligned}
>    $$
>
> 最终得到
> $$
> L = \begin{pmatrix} 1 & 0 & 0 \\ 2 & 1 & 0 \\ 3 & 2 & 1 \end{pmatrix}
> $$
>
> 验证 $A = LL^{\dagger}$：
> $$
> \begin{aligned}
> LL^{\dagger} &= \begin{pmatrix} 1 & 0 & 0 \\ 2 & 1 & 0 \\ 3 & 2 & 1 \end{pmatrix} \begin{pmatrix} 1 & 2 & 3 \\ 0 & 1 & 2 \\ 0 & 0 & 1 \end{pmatrix} = \begin{pmatrix} 1 & 2 & 3 \\ 2 & 5 & 8 \\ 3 & 8 & 14 \end{pmatrix} = A
> \end{aligned}
> $$

### 1.4 Thomas 算法

用于求解三对角矩阵方程组 $A \vec{x} = \vec{d}$，其中
$$
A = \begin{pmatrix} a_1 & c_1 & 0 & \cdots & 0 \\ b_2 & a_2 & c_2 & \cdots & 0 \\ 0 & b_3 & a_3 & \cdots & 0 \\ \vdots & \vdots & \vdots & \ddots & \vdots \\ 0 & 0 & 0 & b_n & a_n \end{pmatrix}, \quad \vec{d} = \begin{pmatrix} d_1 \\ d_2 \\ d_3 \\ \vdots \\ d_n \end{pmatrix}
$$

解法：将 $A$ 分解为 $LU$：
$$
L = \begin{pmatrix} 1 & 0 & 0 & \cdots & 0 \\ l_2 & 1 & 0 & \cdots & 0 \\ 0 & l_3 & 1 & \cdots & 0 \\ \vdots & \vdots & \vdots & \ddots & \vdots \\ 0 & 0 & 0 & l_n & 1 \end{pmatrix}, \quad U = \begin{pmatrix} u_1 & c_1 & 0 & \cdots & 0 \\ 0 & u_2 & c_2 & \cdots & 0 \\ 0 & 0 & u_3 & \cdots & 0 \\ \vdots & \vdots & \vdots & \ddots & \vdots \\ 0 & 0 & 0 & 0 & u_n \end{pmatrix}
$$
其中，$l_i$ 和 $u_i$ 通过以下递推关系计算得到：
$$
\begin{aligned}
u_1 &= a_1, \\
l_i &= \frac{b_i}{u_{i-1}}, \quad i = 2, 3, \ldots, n \\
u_i &= a_i - l_i c_{i-1}, \quad i = 2, 3, \ldots, n
\end{aligned}
$$

此时有 $A\vec{x}=LU\vec{x}=L(U\vec{x})=L\vec{y}=\vec{d}$，然后通过前向替换求解 $L \vec{y} = \vec{d}$，再通过回代求解 $U \vec{x} = \vec{y}$：
$$
\begin{aligned}
y_1 &= d_1, \\
y_i &= d_i - l_i y_{i-1}, \quad i = 2, 3, \ldots, n \\
x_n &= \frac{y_n}{u_n}, \\
x_i &= \frac{y_i - c_i x_{i+1}}{u_i}, \quad i = n-1, n-2, \ldots, 1
\end{aligned}
$$

变式：矩阵的左下角和右上角非零：
$$
A = \begin{pmatrix} a_1 & c_1 & 0 & \cdots & b_1 \\ b_2 & a_2 & c_2 & \cdots & 0 \\ 0 & b_3 & a_3 & \cdots & 0 \\ \vdots & \vdots & \vdots & \ddots & \vdots \\ c_n & 0 & 0 & b_n & a_n \end{pmatrix},\quad \vec{f} = \begin{pmatrix} f_1 \\ f_2 \\ f_3 \\ \vdots \\ f_n \end{pmatrix}
$$

解法：令矩阵除去第一行和第一列后的子矩阵为 $A_1$，对应除去第一行的向量为 $\vec{f}_1$，并设向量 $\vec{u}=(u_2, u_3, \ldots, u_n)^{\mathrm{T}}$ 满足 $A_1 \vec{u} = \vec{f}_1$，则有
$$
\begin{aligned}
\sum_{j=1}^{n} a_{ij} x_j &= f_i,\quad i = 1, 2, \ldots, n \\
\sum_{j=2}^{n} a_{ij} u_j &= f_i, \quad i = 2, 3, \ldots, n
\end{aligned}
$$

两式相减得
$$
\sum_{j=2}^{n} a_{ij} (x_j - u_j) = -a_{i1} x_1, \quad i = 2, 3, \ldots, n
$$

定义 $\vec{v}=(v_2, v_3, \ldots, v_n)^{\mathrm{T}}$ 满足 $v_j = \frac{x_j - u_j}{x_1}$，则 $x_j = u_j + v_j x_1$。此时 $\vec{u}, \vec{v}$ 满足
$$
\begin{aligned}
A_1 \vec{u} &= \vec{f}_1, \\
A_1 \vec{v} &= (-b_2, 0, \ldots, 0, -c_n)^{\mathrm{T}},
\end{aligned}
$$

通过 Thomas 算法分别求解 $\vec{u}, \vec{v}$。再由 $A \vec{x} = \vec{f}$ 的第一行可得
$$
a_{11} x_1 + c_1 x_2 + b_1 x_n = f_1
$$

将 $x_2, x_n$ 用 $\vec{u}, \vec{v}$ 表示后解出 $x_1=\frac{f_1 - c_1 u_2 - b_1 u_n}{a_1 + c_1 v_2 + b_1 v_n}$。

最后，代入 $x_j = u_j + v_j x_1$ 求出其余 $x_j$。

> e.g. 设
> $$
> A = \begin{pmatrix} 40 & 8 & 0 & 2 \\ 1 & 2 & 3 & 0 \\ 0 & 4 & 8 & 3 \\ 1 & 0 & 4 & 8 \end{pmatrix}, \quad \vec{f} = \begin{pmatrix} 0 \\ 12 \\ 31 \\ 16 \end{pmatrix}
> $$
>
> 解方程组 $A \vec{x} = \vec{f}$。
>
> 解：
>
> 1. 去掉第一行和第一列，得到子矩阵
>    $$
>    A_1 = \begin{pmatrix} 2 & 3 & 0 \\ 4 & 8 & 3 \\ 0 & 4 & 8 \end{pmatrix}, \quad \vec{f}_1 = \begin{pmatrix} 12 \\ 31 \\ 16 \end{pmatrix}
>    $$
> 2. 求解 $A_1 \vec{u} = \vec{f}_1,\; A_1 \vec{v} = \begin{pmatrix} -1 \\ 0 \\ -1 \end{pmatrix}$：
>    - 通过 Thomas 算法，得到 $\vec{u} = \begin{pmatrix} 3 \\ 2 \\ 1 \end{pmatrix}, \; \vec{v} = \begin{pmatrix} -\frac{61}{8} \\ \frac{19}{4} \\ -\frac{5}{2} \end{pmatrix}$
> 3. 由第一行方程求解 $x_1$：
>    $$
>    \begin{aligned}
>    40 x_1 + 8 x_2 + 2 x_4 &= 0
>    \end{aligned}
>    $$
>
>    代入 $x_2 = u_2 + v_2 x_1, \; x_4 = u_4 + v_4 x_1$，解得
>    $$
>    x_1 = \frac{0 - 8 \cdot 3 - 2 \cdot 1}{40 + 8 \cdot \frac{-61}{8} + 2 \cdot \left(-\frac{5}{2}\right)} = 1
>    $$
> 4. 计算其余 $x_j$得
>    $$
>    \begin{aligned}
>    \vec{x} = \begin{pmatrix} 1 \\ -\frac{37}{8} \\ \frac{27}{4} \\ -\frac{3}{2} \end{pmatrix}
>    \end{aligned}
>    $$

## 2. 内插

### 2.1 多项式内插

#### 2.1.1 Lagrange 内插

给定 $n+1$ 个数据点 $(x_0, y_0), (x_1, y_1), \ldots, (x_n, y_n)$，Lagrange 内插多项式为
$$
P_n(x) = \sum_{i=0}^{n} y_i L_i(x)
$$
其中，$L_i(x)$ 为 Lagrange 基函数，定义为
$$
L_i(x) = \prod_{\substack{0 \leq j \leq n \\ j \neq i}} \frac{x - x_j}{x_i - x_j}
$$

#### 2.1.2 Newton 内插

给定 $n+1$ 个数据点 $(x_0, y_0), (x_1, y_1), \ldots, (x_n, y_n)$，Newton 内插多项式为
$$
P_n(x) = \sum_{i=0}^{n} a_i n_i(x)
$$
其中，$n_i(x)$ 为 Newton 基函数，定义为
$$
n_i(x) = \prod_{j=0}^{i-1} (x - x_j)
$$
系数 $a_i$ 通过令 $P_n(x_i) = y_i$ 计算得到，具体为
$$
\begin{pmatrix}
1 \\
1 & n_1(x_1) \\
1 & n_1(x_2) & n_2(x_2) \\
\vdots & \vdots & \vdots & \ddots \\
1 & n_1(x_n) & n_2(x_n) & \cdots & n_n(x_n)
\end{pmatrix} \begin{pmatrix} a_0 \\ a_1 \\ a_2 \\ \vdots \\ a_n \end{pmatrix} = \begin{pmatrix} y_0 \\ y_1 \\ y_2 \\ \vdots \\ y_n \end{pmatrix}
$$

### 2.1.3 有理分式插值

给定 $2n+1$ 个数据点 $(x_0, y_0), (x_1, y_1), \ldots, (x_2n, y_2n)$，可以构造有理分式插值多项式
$$
\begin{aligned}
\phi(x_i) &= y_i, \quad i = 0, 1, \ldots, 2n \\
\phi(x_i, x_j) &= \frac{x_i-x_j}{\phi(x_i) - \phi(x_j)}, \\
\phi(x_i, x_j, x_k) &= \frac{x_j - x_k}{\phi(x_i, x_j) - \phi(x_i, x_k)}, \\
\vdots \\
\phi(x_{i_1}, x_{i_2}, \ldots, x_{i_m}) &= \frac{x_{i_{m-1}} - x_{i_m}}{\phi(x_{i_1}, x_{i_2}, \ldots, x_{i_{m-2}}, x_{i_{m-1}}) - \phi(x_{i_1}, x_{i_2}, \ldots, x_{i_{m-2}}, x_{i_m})}
\end{aligned}
$$

最终的有理分式插值多项式为
$$
\begin{aligned}
R(x) = \phi(x_0) + \frac{x-x_0}{\phi(x_0, x_1) + \frac{x-x_1}{\phi(x_0, x_1, x_2) + \frac{x-x_2}{\ddots +\frac{x - x_{2n-1}}{\phi(x_0, x_1, \ldots, x_{2n})}}}}
\end{aligned}
$$

## 3. 数值积分

### 3.1 Newton-Cotes 公式

Newton-Cotes 公式通过在积分区间上使用等距节点的插值多项式来近似定积分。设积分区间为 $[a, b]$，将其划分为 $n$ 个子区间，节点为 $x_i = a + i h$，其中 $h = \frac{b-a}{n}$，则 Newton-Cotes 公式为
$$
\int_a^b f(x) \, dx \approx h\sum_{i=0}^{n} \alpha_i f(x_i)
$$

其中，权重系数 $\alpha_i$ 由以下公式计算得到
$$
\alpha_i = \int_0^n L_i(t) \, dt
$$
这里，$L_i(t)$ 为变量替换 $x = a + t h$ 后 Lagrange 基函数，定义为
$$
L_i(t) = \prod_{\substack{0 \leq j \leq n \\ j \neq i}} \frac{t - j}{i - j}
$$

### 3.2 Gauss 积分公式

#### 3.2.1 一般形式的 Gauss 积分公式

Gauss 积分公式通过选择最佳的节点和权重来最大化积分的精度。对于 $n+1$ 个节点的 Gauss 积分公式，取节点 $x_i$ 为 $n+1$ 阶正交多项式 $P_{n+1}(x)$（如 Legendre 多项式）的 $n+1$ 个根，有下面的结论：

将区间 $[a, b]$ 上的函数 $f(x)$ 按正交多项式展开为
$$
f(x) \approx \sum_{i=0}^{n} \lambda_i P_i(x)
$$
$P_k(x)$ 是 $[a, b]$ 上的正交多项式，满足
$$
P_0(x) = 1, \quad  \text{deg}(P_k) = k, \quad k = 0, 1, \ldots, n, \\
(P_i(x), P_j(x))=\int_a^b P_i(x) P_j(x) w(x) \, dx = 0, \quad i \neq j
$$
其中，$w(x)$ 是权函数。则 Gauss 积分公式为
$$
\int_a^b f(x) w(x) \, dx \approx \sum_{i=0}^{n} \lambda_i \int_a^b P_i(x) w(x) \, dx = \sum_{i=0}^{n} \sum_{j=0}^{n} \lambda_i A_j P_i(x_j) = \sum_{j=0}^{n} A_j f(x_j)
$$
其中，权重 $A_j$ 满足
$$
\int_a^b P_i(x) w(x) \, dx = \sum_{j=0}^{n} A_j P_i(x_j), \quad  i = 0, 1, \ldots, n
$$
权重 $A_j$ 可通过解线性方程组计算得到：
$$
\begin{pmatrix}
P_0(x_0) & P_0(x_1) & \cdots & P_0(x_n) \\
P_1(x_0) & P_1(x_1) & \cdots & P_1(x_n) \\
\vdots & \vdots & \ddots & \vdots \\
P_n(x_0) & P_n(x_1) & \cdots & P_n(x_n)
\end{pmatrix}
\begin{pmatrix}
A_0 \\ A_1 \\ \vdots \\ A_n
\end{pmatrix} = \begin{pmatrix}
\int_a^b P_0(x) w(x) \, dx \\
\int_a^b P_1(x) w(x) \, dx \\
\vdots \\
\int_a^b P_n(x) w(x) \, dx
\end{pmatrix} = \begin{pmatrix}
\int_a^b w(x) \, dx \\
0 \\ \vdots \\ 0
\end{pmatrix}
$$
注意，由于正交多项式的性质 $(P_0(x), P_j(x))=0$（$j \neq 0$），上式右侧的向量中，除第一个元素外，其余元素均为零。

#### 3.2.2 Legendre-Gauss 积分公式

生成正交多项式的方法：Gram-Schmidt 正交化。
$$
\begin{aligned}
P_0(x) &= 1, \\
P_1(x) &= x - \frac{(x, P_0)}{(P_0, P_0)} P_0(x), \\
P_k(x) &= x^k - \sum_{i=0}^{k-1} \frac{(x^k, P_i)}{(P_i, P_i)} P_i(x), \quad k = 1, 2, \ldots, n
\end{aligned}
$$
更简单地，正交多项式满足以下递推关系：
$$
\begin{aligned}
P_0(x) &= 1, \\
P_1(x) &= (x - \alpha_0) P_0(x), \\
P_{k+1}(x) &= (x - \alpha_k) P_k(x) - \beta_k P_{k-1}(x), \quad k = 1, 2, \ldots, n-1
\end{aligned}
$$
其中，$\alpha_k$ 和 $\beta_k$ 由下式计算得到：
$$
\begin{aligned}
\alpha_k &= \frac{(x P_k, P_k)}{(P_k, P_k)},\quad k = 0, 1, \ldots, n-1, \\
\beta_k &= \frac{(P_k, P_k)}{(P_{k-1}, P_{k-1})},\quad k = 1, 2, \ldots, n-1
\end{aligned}
$$
上面得到的正交多项式 $P_n(x)$ 实际上就是 Legendre 多项式，其定义域为 $[-1, 1]$，权函数为 $w(x) = 1$。

由此可得 Legendre-Gauss 积分公式：
$$
\int_{-1}^{1} f(x) \, dx \approx \sum_{j=0}^{n} A_j f(x_j)
$$
其中，节点 $x_j$ 为 Legendre 多项式 $P_{n+1}(x)$ 的 $n+1$ 个根，权重 $A_j$ 通过解线性方程组计算得到：
$$
\begin{pmatrix}
1 & 1 & \cdots & 1 \\
P_1(x_0) & P_1(x_1) & \cdots & P_1(x_n) \\
\vdots & \vdots & \ddots & \vdots \\
P_n(x_0) & P_n(x_1) & \cdots & P_n(x_n)
\end{pmatrix}
\begin{pmatrix}
A_0 \\ A_1 \\ \vdots \\ A_n
\end{pmatrix} = \begin{pmatrix}
2 \\ 0 \\ \vdots \\ 0
\end{pmatrix}
$$

## 4. 方程求根

### 4.1 二分法

二分法适用于连续函数 $f(x)$ 在区间 $[a, b]$ 上有根的情况，即 $f(a)f(b) < 0$。

它是线性收敛的：$\epsilon_{i+1} \approx \frac{1}{2} \epsilon_i$，其中 $\epsilon_i = |x_i - x^*|$ 是第 $i$ 次迭代的误差，$x^*$ 是实际根。

### 4.2 不动点迭代法

不动点迭代法通过将方程 $f(x) = 0$ 转化为 $x = g(x)$ 的形式来求解根。选择适当的函数 $g(x)$，使得迭代过程收敛。

#### 4.2.1 Newton-Raphson 方法

Newton-Raphson 方法是一种常用的不动点迭代法，适用于求解非线性方程 $f(x) = 0$。迭代公式为
$$
x_{n+1} = x_n - \frac{f(x_n)}{f'(x_n)}
$$

它是二次收敛的：$\epsilon_{i+1} \approx C \epsilon_i^2$，其中 $\epsilon_i = |x_i - x^*|$ 是第 $i$ 次迭代的误差，$x^*$ 是实际根，$C$ 是与函数 $f(x)$ 及其导数有关的常数。

#### 4.2.2 割线法

割线法是一种不需要计算导数的迭代方法，适用于求解非线性方程 $f(x) = 0$。迭代公式为
$$
x_{n+1} = x_n - f(x_n) \frac{x_n - x_{n-1}}{f(x_n) - f(x_{n-1})}
$$

它的收敛速度介于线性和二次之间，通常称为超线性收敛：$\epsilon_{i+1} \approx C \epsilon_i^{\phi}$，其中 $\phi=\frac{1+\sqrt{5}}{2} \approx 1.618$（黄金比例），$\epsilon_i = |x_i - x^*|$ 是第 $i$ 次迭代的误差，$x^*$ 是实际根，$C$ 是与函数 $f(x)$ 有关的常数。

## 5. 矩阵本征值、本征矢量

### 5.1 QR 算法

QR 算法是一种用于计算矩阵本征值和本征矢量的迭代方法。其基本思想是将矩阵 $T_k$ 分解为 QR 形式，然后通过迭代更新矩阵，逐渐逼近上三角矩阵，从而得到本征值：
$$
\begin{aligned}
T_0 &= A, \\
T_k &= Q_k R_k, \\
R_k Q_k &= T_{k+1}, \quad k = 0, 1, 2, \ldots
\end{aligned}
$$
当 $k$ 足够大时，$T_k$ 将近似为上三角矩阵，其对角线元素即为矩阵 $A$ 的本征值。对应的本征矢量可以通过累积 QR 分解中的正交矩阵 $Q_k$ 来计算。

#### 5.1.1 Householder 变换

Householder 变换是一种用于将矩阵转换为上 Hessenberg 形式的正交变换。上 Hessenberg 矩阵是一种接近上三角矩阵的形式，便于后续的 QR 分解：
$$
\begin{pmatrix}*
& * & * & * & * \\
* & * & * & * & * \\
0 & * & * & * & * \\
0 & 0 & * & * & * \\
0 & 0 & 0 & * & * \end{pmatrix}
$$
Householder 变换通过构造一个反射矩阵，将矩阵的某一列下方的元素消去，从而实现矩阵的简化。

对于给定的向量 $\vec{x} \in \mathbb{R}^n$，Householder 变换定义为
$$
H = I - 2 \frac{\vec{u} \vec{u}^{\mathrm{T}}}{\vec{u}^{\mathrm{T}} \vec{u}}
$$
其中
$$
\vec{u} = \vec{x} - \alpha \vec{e_1}, \quad \alpha = \pm \|\vec{x}\|_2
$$
这里，$\vec{e_1}$ 是单位向量 $(1, 0, \ldots, 0)^{\mathrm{T}}$。通过适当选择 $\alpha$ 的符号，可以避免数值不稳定性。作用到矩阵 $A$ 上，得到上 Hessenberg 形式的矩阵 $A'$：
$$
A' = H_{(n-2)}^T \cdots H_{(1)}^T A H_{(1)} \cdots H_{(n-2)}
$$

> 另外，Householder 变换还可以直接对矩阵进行QR分解。对于矩阵 $A$，通过一系列 Householder 变换，可以将其转换为上三角矩阵 $R$，同时累积这些变换得到正交矩阵 $Q$，使得 $A = QR$：
> $$
> R = H_{(m-1)} \cdots H_{(1)} H_{(0)} A, \quad Q = H_{(0)}^T \cdots H_{(m-1)}^T
> $$
>
> 但这种方法在计算上不如先转换为 Hessenberg 形式再利用 Givens 变换进行 QR 分解高效，因为我们需要进行多次迭代，与 Hessenberg 形式下的 QR 分解相比，计算开销更大。

#### 5.1.2 Givens 变换

Givens 变换是一种用于对矩阵中特定元素置零的正交变换。它通过旋转二维子空间来实现对矩阵的简化，适用于将矩阵转换为上三角形式。

对于给定的向量 $\begin{pmatrix} a \\ b \end{pmatrix}$，Givens 变换定义为
$$
G = \begin{pmatrix} c & s \\ -s & c \end{pmatrix}
$$
其中
$$
c = \frac{a}{\sqrt{a^2 + b^2}}, \quad s = \frac{b}{\sqrt{a^2 + b^2}}
$$
通过选择适当的 $c$ 和 $s$，可以将向量的第二个分量置零，从而实现对矩阵的简化。

综合使用 Householder 变换和 Givens 变换，可以高效地实现 QR 分解：

1. 使用 Householder 变换对矩阵 $A$ 进行预处理，将其转换为上 Hessenberg 形式 $T_0$。
2. 在上 Hessenberg 矩阵 $T_0$ 上应用 Givens 变换，逐步将其转换为上三角矩阵 $R_0=G_{n-1} \cdots G_1 T_0$，得到正交矩阵 $Q_0 = G_1^T \cdots G_{n-1}^T$，即 $T_0 = Q_0 R_0$，进行一次 QR 分解。
3. 令 $T_1 = R_0 Q_0$，重复步骤 2，直到矩阵 $T_k$ 收敛为上三角矩阵。
4. 最终，矩阵 $T_k$ 的对角线元素即为矩阵 $A$ 的本征值，对应的本征矢量可以通过累积所有正交矩阵 $Q_k$ 来计算。

#### 5.1.3 Sturm 序列

Sturm 序列是一种用于计算实对称三对角矩阵本征值的工具。对于实对称三对角矩阵
$$
A=\begin{pmatrix}d_1 & b_1 \\
b_1 & d_2 & b_2 \\
& b_2 & d_3 & \ddots \\
& & \ddots & \ddots & b_{n-1} \\
& & & b_{n-1} & d_n \end{pmatrix}
$$
其 Sturm 序列定义为一系列多项式 $\{p_0(x), p_1(x), \ldots, p_n(x)\}$，其中
$$
\begin{aligned}
p_0(x) &= 1, \\
p_1(x) &= d_1 - x, \\
\ldots \\
p_k(x) &= (d_k - x) p_{k-1} - b_{k-1}^2 p_{k-2} = \det{T_i}, \quad k = 2, 3, \ldots, n
\end{aligned}
$$
这里，$T_i$ 是矩阵 $T =A - xI$ 的前 $i \times i$ 子矩阵。

Sturm 定理指出，对于任意实数 $\mu$，多项式序列 $\{p_0(\mu), p_1(\mu), \ldots, p_n(\mu)\}$ 中符号变化的次数 $V(\mu)$ 就是矩阵 $A$ 的本征值中严格小于 $\mu$ 的个数。

因此，可以结合二分法和 Sturm 序列来定位实对称矩阵的本征值。通过计算不同 $\mu$ 值下的符号变化次数，可以逐步缩小包含本征值的区间，最终确定所有本征值的位置:

1. 选择一个包含所有本征值的初始区间 $[\alpha, \beta]$。一般来说，可以通过 Gershgorin 圆盘定理来确定这个区间：
   $$
   \alpha = \min_i (d_i - R_i), \quad \beta = \max_i (d_i + R_i)
   $$
   其中，$R_i = \sum_{j \neq i} |a_{ij}|=|b_{i-1}| + |b_i|$ 是矩阵 $A$ 的第 $i$ 行的非对角元素之和。

2. 计算区间中点 $\mu = \frac{\alpha + \beta}{2}$，并计算 Sturm 序列在 $\mu$ 处的符号变化次数 $V(\mu)$。

3. 根据 $V(\mu)$ 的值，调整区间：
   - 如果 $V(\mu) < k$，则说明第 $k$ 个本征值在区间 $[\mu, \beta]$ 内，更新 $\alpha = \mu$。
   - 如果 $V(\mu) \geq k$，则说明第 $k$ 个本征值在区间 $[\alpha, \mu]$ 内，更新 $\beta = \mu$。

4. 重复步骤 2 和 3，直到区间 $[\alpha, \beta]$ 足够小，满足预设的精度要求。则区间内的中点即为第 $k$ 个本征值的近似值。

## 6. 常微分方程数值解法

### 6.1 Euler 方法

Euler 方法是一种用于求解初值问题的数值方法。对于常微分方程初值问题
$$
\begin{aligned}
\frac{du}{dt} &= f(t, u), \\
u(t_0) &= u_0
\end{aligned}
$$
Euler 方法通过在每个时间步长 $h$ 上使用差分近似来更新解。

1. 向前差分（显式 Euler 方法）：
$$
\begin{aligned}
&u'_{n} \approx \frac{u_{n+1} - u_n}{h} \\
\implies &u_{n+1} = u_n + h f(t_n, u_n)
\end{aligned}
$$

2. 向后差分（隐式 Euler 方法）：
$$
\begin{aligned}
&u'_{n+1} \approx \frac{u_{n+1} - u_n}{h} \\
\implies &u_{n+1} = u_n + h f(t_{n+1}, u_{n+1})
\end{aligned}
$$

3. 中心差分（隐式中点法；两步法）：
$$
\begin{aligned}
&u'_{n} \approx \frac{u_{n+1} - u_{n-1}}{2h} \\
\implies &u_{n+1} = u_{n-1} + 2h f\left(t_n, u_n \right)
\end{aligned}
$$

### 6.2 Runge-Kutta 方法

Runge-Kutta 方法是一类用于求解常微分方程初值问题的高阶数值方法。

对于 $n$ 阶 RK 法，设 $u_{k+1} = u_k + h \sum_{i=1}^{n} c_i k_i$，其中 $k_1 = f(t_k, u_k)$，$k_i = f\left(t_k + a_i h, u_k + h \sum_{j=1}^{i-1} b_{ij} k_j\right)$，$i = 2, 3, \ldots, n$，则系数 $a_i, b_{ij}, c_i$ 满足以下条件：
$$
\begin{aligned}
\sum_{i=1}^{n} c_i &= 1, \\
\sum_{i=2}^{n} c_i a_i &= \frac{1}{2}, \\
\sum_{j=1}^{i-1} b_{ij} &= a_i, \quad i = 2, 3, \ldots, n
\end{aligned}
$$
可以选取符合上述条件的系数来构造不同阶数的 RK 方法。

### 6.3 边值问题

边值问题通常涉及在区间 $[a, b]$ 上求解常微分方程，并满足特定的边界条件。设有二阶常微分方程
$$
\begin{aligned}
\frac{d^2 y}{dx^2} &= f(x, y, y'), \\
a_0 y(a) + b_0 y'(a) &= \alpha, \\
a_1 y(b) + b_1 y'(b) &= \beta
\end{aligned}
$$
可以通过有限差分法将其离散化。将区间 $[a, b]$ 划分为 $N$ 个子区间，节点为 $x_i = a + i h$，其中 $h = \frac{b-a}{N}$。使用中心差分近似二阶导数和一阶导数：
$$
\begin{aligned}
\frac{d^2 y}{dx^2} \bigg|_{x=x_i} &\approx \frac{y_{i+1} - 2y_i + y_{i-1}}{h^2} \\
\frac{dy}{dx} \bigg|_{x=x_i} &\approx \frac{y_{i+1} - y_{i-1}}{2h}
\end{aligned}
$$
最终可以得到一个三对角矩阵线性方程组 $A \vec{y} = \vec{d}$，形式为
$$
\begin{aligned}
\begin{pmatrix}
b_1 & c_1 \\
a_2 & b_2 & c_2 \\
& a_3 & b_3 & c_3 \\
& & \ddots & \ddots & \ddots \\
& & & a_{N-1} & b_{N-1} & c_{N-1} \\
& & & & a_N & b_N
\end{pmatrix} \begin{pmatrix}
y_1 \\ y_2 \\ y_3 \\ \vdots \\ y_{N-1} \\ y_N
\end{pmatrix} = \begin{pmatrix}
d_1 \\ d_2 \\ d_3 \\ \vdots \\ d_{N-1} \\ d_N
\end{pmatrix}
\end{aligned}
$$
其中，矩阵 $A$ 和向量 $\vec{d}$ 由离散化后的方程和边界条件构成。通过求解该线性方程组，可以得到在各个节点处的近似解 $y_i$

## 7. 偏微分方程数值解法

### 7.1 有限差分法

有限差分法通过在空间和时间上使用差分近似来求解偏微分方程。以一阶对流方程为例：
$$
\begin{aligned}
\frac{\partial u}{\partial t} + a \frac{\partial u}{\partial x} &= 0, \quad a \ne 0,
\end{aligned}
$$
将空间和时间离散化，节点为 $x_i = i h$，$t_n = n \Delta t$，其中 $h$ 是空间步长，$\Delta t$ 是时间步长。记 $u_i^n = u(x_i, t_n)$，则可以使用以下差分格式进行数值求解：

1. Upwind（迎风）差分：
   $$
   \frac{u_i^{n+1} - u_i^n}{\Delta t} + a \frac{u_i^n - u_{i-1}^n}{h} = 0, \quad a > 0 \\
   \frac{u_i^{n+1} - u_i^n}{\Delta t} + a \frac{u_{i+1}^n - u_i^n}{h} = 0, \quad a < 0
   $$
   von Neumann 稳定性分析：令 $u_i^n = \xi^n e^{j \omega i h}$，代入差分格式，得到增益因子
   $$
    \xi = 1 - \frac{a \Delta t}{h} (1 - e^{-j \omega h}), \quad a > 0 \\
    \xi = 1 - \frac{a \Delta t}{h} (e^{j \omega h} - 1), \quad a < 0
   $$
    计算 $|\xi|^2$，得到
   $$
    |\xi|^2 = 1 - 2 \left| \frac{a \Delta t}{h} \right| (1 - \cos(\omega h)) + \left( \frac{a \Delta t}{h} \right)^2 (1 - \cos(\omega h))^2
   $$
   当 $| \frac{a \Delta t}{h} | \leq 1$ 时，$|\xi|<1$，格式稳定，称为 CFL 条件。
2. 中心差分：
    $$
    \frac{u_i^{n+1} - u_i^n}{\Delta t} + a \frac{u_{i+1}^n - u_{i-1}^n}{2h} = 0
    $$
    计算表明 $|\xi|^2 = 1 + \left( \frac{a \Delta t}{h} \right)^2 \sin^2(\omega h)$，无论 $\Delta t$ 和 $h$ 取何值，均有 $|\xi|>1$，格式不稳定。
3. Lax 格式：
    $$
    \frac{u_i^{n+1} - \frac{1}{2}(u_{i+1}^n + u_{i-1}^n)}{\Delta t} + a \frac{u_{i+1}^n - u_{i-1}^n}{2h} = 0
    $$
    计算表明 $|\xi|^2 = 1 - (1 - \left( \frac{a \Delta t}{h} \right)^2) \sin^2(\omega h)$，当 $| \frac{a \Delta t}{h} | \leq 1$ 时，$|\xi| \leq 1$，格式稳定。
4. 蛙跳格式：
    $$
    \frac{u_i^{n+1} - u_i^{n-1}}{2 \Delta t} + a \frac{u_{i+1}^n - u_{i-1}^n}{2h} = 0
    $$
    计算表明当 $| \frac{a \Delta t}{h} | \leq 1$ 时，格式稳定。

### 7.2 Crank-Nicolson 方法

Crank-Nicolson 方法是一种用于求解时间依赖型偏微分方程的隐式有限差分方法，可以提高数值解的稳定性和精度。

在上面的中心差分基础上，Crank-Nicolson 方法通过在时间上采用中点近似来改进数值格式：
$$
\frac{u_i^{n+1} - u_i^n}{\Delta t} + a \frac{u_{i+1}^{n+1} - u_{i-1}^{n+1} + u_{i+1}^n - u_{i-1}^n}{4h} = 0
$$
通过整理，可以得到一个线性方程组，用于计算时间步长 $n+1$ 时刻的数值解 $u_i^{n+1}$。Crank-Nicolson 方法具有二阶精度，并且在适当的条件下是无条件稳定的：
$$
|\xi|^2 = \frac{1 + \frac{1}{2} \left( \frac{a \Delta t}{h} \right)^2 \sin^2(\omega h)}{1 - \frac{1}{2} \left( \frac{a \Delta t}{h} \right)^2 \sin^2(\omega h)} = 1
$$

## 8. 随机数与概率论

### 8.1 伪随机数生成

伪随机数生成器通过确定性算法生成近似随机的数值序列。常用的方法包括平方取中法、线性同余法和反馈移位寄存器法等等：

1. 平方取中法：选择一个初始种子数 $x_0$，通过平方该数并取中间几位作为下一个随机数 $x_{n+1}$。例如，若 $x_n = 1234$，则 $x_n^2 = 1522756$，取中间四位 $2275$ 作为 $x_{n+1}$。

2. 线性同余法：通过递推关系 $x_{n+1} = (a x_n + c) \mod m$ 生成随机数序列，其中 $a, c, m$ 是预先选择的常数，$x_0$ 是初始种子数。一般而言，选择 $m=2^{L}$，$a=4q+1$，$c=2r+1$ （$q, r$ 为整数）可以获得较好的周期性和均匀性。

3. 反馈移位寄存器法：利用线性反馈移位寄存器（LFSR）生成伪随机数序列。通过对寄存器中的位进行线性组合并反馈到寄存器的输入端，产生新的位序列。该方法适用于硬件实现，具有较高的速度和较长的周期。

### 8.2 常见概率分布

1. 均匀分布（Uniform Distribution）：在区间 $[a, b]$ 上均匀分布的随机变量 $X$ 的概率密度函数为
   $$
   f(x) = \begin{cases}
   \frac{1}{b-a}, & a \leq x \leq b \\
   0, & \text{otherwise}
   \end{cases}
   $$
   其期望值和方差分别为
   $$
   E(X) = \frac{a+b}{2}, \quad    V(X) = \frac{(b-a)^2}{12}
   $$

2. Poisson 分布（Poisson Distribution）：描述单位时间或单位空间内事件发生次数的概率分布。随机变量 $X$ 服从参数为 $\lambda$ 的 Poisson 分布，其概率质量函数为
   $$
   P(X=k) = \frac{\lambda^k e^{-\lambda}}{k!}, \quad k = 0, 1, 2, \ldots
   $$
   其期望值和方差均为
   $$
   E(X) = \lambda, \quad V(X) = \lambda
   $$

3. 正态分布（Normal Distribution, Gaussian Distribution）：描述连续随机变量的概率分布，具有钟形曲线的特征。随机变量 $X$ 服从均值为 $\mu$，方差为 $\sigma^2$ 的正态分布，其概率密度函数为
   $$
   f(x) = \frac{1}{\sigma \sqrt{2\pi}} e^{-\frac{(x-\mu)^2}{2\sigma^2}}, \quad -\infty < x < \infty
   $$
   其期望值和方差分别为
   $$
   E(X) = \mu, \quad V(X) = \sigma^2
   $$

### 8.3 随机抽样方法

1. 直接抽样法，又叫逆变换抽样法（Inverse Transform Sampling）。首先生成一个均匀分布的随机数 $u \sim U(0, 1)$，然后通过目标分布的累积分布函数（CDF）的反函数 $F^{-1}(u)$ 来获得样本值 $x = F^{-1}(u)$。
   > e.g. 对于指数分布 $f(x)=\lambda e^{-\lambda x}$，CDF 为 $F(x) = 1 - e^{-\lambda x}$，其反函数为 $F^{-1}(u) = -\frac{1}{\lambda} \ln(1-u)$，则可抽取均匀分布的 $u \sim U(0,1)$，得到样本值 $x=-\frac{1}{\lambda} \ln(1-u)$。

2. 变换抽样法，适用于从复杂分布中抽样。通过对简单分布的随机变量进行变换，得到目标分布的样本值：
   $$
   f(x)dx = f(g(y)) \left| \frac{dg(y)}{dy} \right| dy=\phi(y)dy,\quad x = g(y)
   $$
   设 $\eta,\delta$ 分别满足概率分布 $f(x),\phi(y)$，则可通过从 $\phi(y)$ 抽样得到 $\delta$，再通过 $\eta = g(\delta)$ 得到目标分布的样本值 $\eta$。

   > e.g. 二维标准正态分布:
   > $$
   > f(x, y)=\frac{1}{2\pi}e^{-\frac{x^2+y^2}{2}}
    > $$
    > 令 $x=\sqrt{-2\ln{u}}\cos{(2\pi v)},\;y=\sqrt{-2\ln{u}}\sin{(2\pi v)}$，其中 $u,v \sim U(0,1)$，则有
    > $$
    > f(x,y)dxdy = f(x(u,v), y(u,v))\left|\frac{\partial(x,y)}{\partial(u,v)}\right|dudv = dudv
    > $$
    > 则抽取 $u,v \sim U(0,1)$，即可得到二维标准正态分布的样本值 $(x,y)$。

3. 舍选抽样法
   - 第一类：适用于目标分布 $f(x)$ 有界且较为简单的情况。记 $L=\max_x f(x)$，步骤如下：
     1. 生成一个均匀随机数 $x^* \sim U(a, b)$，其中 $[a, b]$ 是 $f(x)$ 的定义域。
     2. 生成一个均匀随机数 $u \sim U(0, 1)$。
     3. 如果 $u \leq \frac{f(x^*)}{L}$，则接受样本 $x^*$；否则拒绝样本，返回步骤 1。

     抽样效率 $E=\frac{1}{L(b-a)}$。

   - 第二类：适用于目标分布 $f(x)$ 更复杂的情况。设有一个易于抽样的辅助分布 $h(x)$，满足 $f(x) \leq L h(x)$，其中 $L$ 是常数。定义 $g(x)=\frac{f(x)}{L h(x)}$，则有 $f(x)=L g(x) h(x)$。步骤如下：
     1. 从辅助分布 $h(x)$ 中生成样本 $x^*$。
     2. 生成一个均匀随机数 $u \sim U(0, 1)$。
     3. 如果 $u \leq g(x^*)$，则接受样本 $x^*$；否则拒绝样本，返回步骤 1。

     抽样效率 $E=\frac{1}{L}$。
     > e.g. 标准正态分布 $f(x)=\frac{1}{\sqrt{2\pi}}e^{-\frac{x^2}{2}}$，选择辅助分布 $h(x)=e^{-x}$，则有 $L=\sqrt{\frac{e}{2\pi}}$，$g(x)=\frac{f(x)}{L h(x)}=e^{-\frac{(x- 1)^2}{2}}$。
     >
     > 先对 $x \geq 0$ 抽样，生成 $x^*$ 服从指数分布 $h(x)$，然后生成均匀随机数 $u$，若 $u \leq g(x^*)$ 则接受样本 $x^*$，否则拒绝样本并重新抽样。对于 $x < 0$ 的情况，可以对称处理，即生成 $-x^*$ 作为样本值。

4. 加抽样，适用于目标分布 $f(x)$ 具有多个峰值的情况。设有 $m$ 个子分布 $h_i(x)$，满足 $f(x) = \sum_{i=1}^{m} p_i h_i(x)$，其中 $p_i \geq 0$ 且 $\sum_{i=1}^{m} p_i = 1$。步骤如下：
   1. 生成一个均匀随机数 $u \sim U(0, 1)$，根据 $u$ 的值选择子分布 $h_k(x)$： $\sum_{i=1}^{k-1} p_i < u \leq \sum_{i=1}^{k} p_i$。
   2. 从选定的子分布 $h_k(x)$ 中生成样本 $x^*$。

   该方法通过组合多个简单分布来近似复杂分布，从而提高抽样效率和准确性。
   > e.g. 均匀球壳抽样：在三维空间中，均匀地从内外半径分别为 $R_0$ 和 $R_1$ 的球壳中抽样。则概率密度函数为
   > $$
   > f(r)=\frac{3r^2}{R_1^3 - R_0^3}, \quad R_0 \leq r \leq R_1
   > $$
   > 设 $r=(R_1 - R_0)x + R_0$，则有
   > $$
   > f(r)dr=g(x)dx=\left[ \frac{(R_1 - R_0)^2}{\lambda}\cdot 3x^2 + \frac{3(R_1 - R_0)R_0}{\lambda} \cdot 2x + \frac{3R_0^2}{\lambda} \right] dx
   > $$
   > 其中，$\lambda = R_1^2+R_0^2+R_1R_0$。记 $p_1 = \frac{(R_1 - R_0)^2}{\lambda}$，$p_2 = \frac{3(R_1 - R_0)R_0}{\lambda}$，$p_3 = \frac{3R_0^2}{\lambda}$，则 $f(r) = p_1 h_1(x) + p_2 h_2(x) + p_3 h_3(x)$，其中 $h_1(x) = 3x^2$，$h_2(x) = 2x$，$h_3(x) = 1$。可以通过加抽样方法从 $h_1, h_2, h_3$ 中抽样来获得均匀球壳分布的样本值:
   > 1. 生成均匀随机数 $u \sim U(0, 1)$，根据 $u$ 的值选择子分布 $h_k(x)$：
   >    - 若 $0 < u \leq p_1$，选择 $h_1(x)$；
   >    - 若 $p_1 < u \leq p_1 + p_2$，选择 $h_2(x)$；
   >    - 若 $p_1 + p_2 < u \leq 1$，选择 $h_3(x)$。
   > 2. 从选定的子分布 $h_k(x)$ 中生成样本 $x^*$：
   >    - 若选择 $h_1(x)$，则 $x^* = u^{1/3}$；
   >    - 若选择 $h_2(x)$，则 $x^* = \sqrt{u}$；
   >    - 若选择 $h_3(x)$，则 $x^* = u$。
   > 3. 最终样本值为 $r = (R_1 - R_0)x^* + R_0$。
