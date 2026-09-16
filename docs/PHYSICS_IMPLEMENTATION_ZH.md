# 气泡溶解求解器物理与数值实现

本文用于论文“物理模型”和“数值方法”部分。公式对应当前代码实现，但论文仍须补充参数来源、实验条件、不确定度和标定过程。

## 1. 适用范围

模型描述孔隙网络中静止被困气泡向流动水相的单向溶解，包括：

- 液相压力驱动流动；
- 溶解气体的对流与扩散；
- 气液界面传质；
- 自由气体摩尔库存；
- 毛细压力闭合；
- 气泡消失active-set；
- 压力或总入口流量边界；
- 摩尔和液体体积守恒审计。

当前不包括气相网络流动、气泡迁移、聚并、破裂、再侵入、动态接触角、多组分分馏和热效应。因此应称为“静止被困气泡溶解孔隙网络模型”，而非通用动态两相流模型。

## 2. 几何与控制体

每个真实控制体具有体积 `Vi`、等效半径 `ri`、坐标、液相饱和度 `Sw,i` 和气相饱和度：

\[
S_{g,i}=1-S_{w,i}.
\]

输入长度、面积和体积转换为SI：

```text
length_to_m = 1e-6
area_to_m2 = 1e-12
volume_to_m3 = 1e-18
```

虚拟入口和出口节点只代表储层，不进入样品总体积、摩尔库存和饱和度统计。

零距离边通过union-find构成cluster。cluster成员共享液相压力和浓度；正体积含气成员仍保留各自的气相压力、饱和度、气体库存和界面面积。零体积massless junction没有积累项。

## 3. 状态变量与压力表示

无气泡控制体：

\[
\mathbf u_i=(P_{l,i},C_i).
\]

含气控制体：

\[
\mathbf u_i=(P_{l,i},C_i,P_{g,i},S_{w,i}).
\]

`Pl` 和 `Pg` 是相对压力。绝对压力为：

\[
P_l^{abs}=P_{ref}+P_l,
\qquad
P_g^{abs}=P_{ref}+P_g.
\]

EOS和Henry关系使用绝对压力；线性流量和毛细压差使用相对压力即可，因为固定背景压力相消。

## 4. 液体流动

边 `e=(i,j)` 的Poiseuille基准导流系数：

\[
G_{0,e}=\frac{\pi r_e^4}{8\mu_l L_e}.
\]

边饱和度取两端平均：

\[
\bar S_{w,e}=\frac{S_{w,i}+S_{w,j}}{2}.
\]

有效液相导流系数和体积流量：

\[
G_e=G_{0,e}\bar S_{w,e}^{m},
\qquad
Q_{ij}=G_e(P_{l,i}-P_{l,j}).
\]

`m=kr_exponent_m` 是待确认闭合参数。

## 5. 溶解组分输运

对流采用一阶上风：

\[
F_{adv,ij}=Q_{ij}C_{upwind}.
\]

湿润扩散面积：

\[
A_{w,e}=A_e\bar S_{w,e}^{\alpha}.
\]

扩散通量：

\[
F_{diff,ij}=\frac{D_m}{\tau}
\frac{A_{w,e}}{L_e}(C_i-C_j).
\]

`alpha=diffusion_area_exponent_alpha` 和 `tau=tortuosity` 需要物理标定。

## 6. 气液界面传质

局部气体体积：

\[
V_{g,i}=V_i(1-S_{w,i}).
\]

当前采用等体积球形界面面积：

\[
A_{gl,i}=(36\pi)^{1/3}V_{g,i}^{2/3}.
\]

界面摩尔传质率：

\[
J_{gl,i}=M_{tr}k_LA_{gl,i}
\left(H_{cp}P_{g,i}^{abs}-C_i\right).
\]

其中 `Mtr=mass_transfer_multiplier`。科研物理算例必须使用 `Mtr=1`；`100` 或 `1000` 仅用于加速预览，不能用于解释真实溶解时间。

同一个 `Jgl` 以相反符号进入溶解组分残量和自由气体残量，因此相间传质在总组分守恒中内部抵消。

## 7. 自由气体状态方程

\[
n_{g,i}=\frac{P_{g,i}^{abs}V_i(1-S_{w,i})}{ZRT}.
\]

当前采用常数压缩因子 `Z`。气体体积变化不等于摩尔数变化，因为压力和毛细压力也会改变体积。

## 8. 毛细压力

当前Qin形式：

\[
P_{c,i}=\frac{2\sigma\cos\theta}
{r_i[1-\exp(-bS_{w,i})]}.
\]

闭合残量：

\[
R_{c,i}=P_{g,i}-P_{l,i}-P_{c,i}=0.
\]

小 `bSw` 时使用级数稳定计算 `1-exp(-bSw)`。接触角、`b` 和将输入半径作为毛细有效半径的处理尚未实验标定。

## 9. Backward-Euler控制方程

液体体积残量：

\[
R_{l,i}=\frac{V_i(S_{w,i}^{n+1}-S_{w,i}^{n})}{\Delta t}
+\sum_j Q_{ij}^{n+1}=0.
\]

massless junction省略积累项。

溶解组分残量：

\[
R_{d,i}=\frac{V_i(S_{w,i}^{n+1}C_i^{n+1}-S_{w,i}^{n}C_i^n)}{\Delta t}
+\sum_j(F_{adv,ij}^{n+1}+F_{diff,ij}^{n+1})
-J_{gl,i}^{n+1}=0.
\]

自由气体残量：

\[
R_{g,i}=\frac{n_{g,i}^{n+1}-n_{g,i}^{n}}{\Delta t}
+J_{gl,i}^{n+1}=0.
\]

代码使用代数等价增量形式减轻相近库存相减的浮点消减误差。

## 10. 边界条件

### 10.1 压力Dirichlet

\[
P_l=P_{l,in}\quad\text{或}\quad P_l=P_{l,out}.
\]

脱气水入口使用 `C_in=0`。出口正常流出时使用内部上游浓度，出口扩散采用零梯度；只有出口回流才使用 `C_backflow`。

### 10.2 总入口流量约束

流量入口增加一个共享入口压力未知量 `Pin` 和一个全局方程：

\[
R_Q=\sum_{e\in\Gamma_{in}}Q_e-Q_{in,target}=0.
\]

各入口喉流量仍为：

\[
Q_e=G_e(P_{in}-P_i),
\]

因此各喉流量由方程自然分配，不是平均指定。Jacobian包含该全局残量对入口压力、内部压力和饱和度的导数。出口继续固定压力。配置中的入口压力在流量模式下仅是Newton初始猜测。

## 11. 局部反流

总入口流量为正不保证每条入口边正向；出口也可能局部回流。程序分别记录：

```text
inlet_reverse_flow_edges
inlet_reverse_flow_m3_s
outlet_backflow_edges
outlet_backflow_m3_s
```

上风处理为：

- 入口正向流入：使用 `C_in`；
- 入口反流：使用内部浓度；
- 出口正常流出：使用内部浓度；
- 出口回流：使用 `C_backflow`。

`C_backflow=0` 代表出口外部储层也是脱气水，是否合理取决于实验装置。

## 12. 气泡消失active-set

当：

\[
S_{g,i}<S_{g,off}
\quad\text{或}\quad
n_{g,i}<n_{g,off}
\]

求解器定位事件步。只有事件步通过Newton、状态边界和守恒后才执行：

```text
active_gas = false
permanently_inactive = true
Sw = 1
```

阈值处剩余自由气体转入cluster液相；饱和度置1产生的体积调整进入active-set体积账本。消失后禁止重新激活。

微小气泡不会永久控制全局时间步：事件邻域泡体和摩尔占比过小泡体不参与全局气体时间尺度；其余泡体使用时间尺度分位数而非最小值。

### 守恒型微小气泡批量退休

可选配置 `active_set.batch_microbubble_retirement` 默认关闭。开启后，只有单泡自由气体占比和同一接受步候选总占比同时低于硬上限的阈值穿越气泡，才允许在完整Backward-Euler步末尾批量退出active-set。显著气泡、批量总量超限和负自由气体预测仍使用精确事件定位。

批量路径不截断库存：每个退休气泡在隐式步末的正自由气体摩尔数逐个加入所属cluster的溶解库存，释放的气体体积逐个进入active-set液体体积账本，随后执行与普通步完全相同的摩尔和体积守恒硬门槛。`gas_disappearance_events.tsv`对每个节点单独记录估算穿越时间、接受时间、退休前后摩尔数、气体体积、区段和批量占比。

## 13. Newton和自动微分

`Rl/Rd/Rg/Rc` 由同一套FADBAD++反向自动微分生成稀疏Jacobian。每步执行：

1. 更新上风方向和传质active状态；
2. 组装残量/Jacobian；
3. 缩放线性系统；
4. Eigen、AMGX或PARDISO求解增量；
5. fraction-to-boundary和Armijo线搜索；
6. 检查残量和更新量双门槛；
7. 检查状态、非负气体和守恒；
8. 接受状态并执行active-set。

接受条件：

\[
\|R_{scaled}\|_\infty\le\varepsilon_R,
\]

\[
\|\Delta u_{scaled}\|_\infty
\le\varepsilon_u(1+\|u_{scaled}\|_\infty).
\]

## 14. 自适应时间步

```text
Newton迭代少：增长dt
Newton迭代正常：保持dt
Newton迭代多：缩小dt
Newton/线性/状态/守恒失败：拒绝并重算
```

气体时间尺度使用有效泡体的配置分位数，并排除事件邻域和摩尔占比过小的泡体。事件后从事件前请求步长的安全比例恢复，避免从微小事件步长长期爬升。

近零局部浓度不用于 `C/|dC/dt|` 全局限制。浓度变化限制默认关闭。

## 15. 摩尔守恒

\[
N_0+N_{in}-N_{out}+N_{ledger}
=N_{free}+N_{dissolved}+E_N.
\]

入口和出口累计量均按正值记录。`Nledger` 表示active-set数值账本。

## 16. 液体体积守恒

\[
V_{l,0}+V_{l,in}-V_{l,out}+V_{ledger}
=V_l(t)+E_V.
\]

气泡缩小时液体占据空间增加，因此通常：

\[
Q_{in}-Q_{out}\approx\frac{dV_l}{dt}.
\]

这解释了溶解阶段入口和出口流量不完全相等。

## 17. 科学统计

真实样品总体气体饱和度：

\[
\bar S_g=\frac{\sum_iV_i(1-S_{w,i})}{\sum_iV_i}.
\]

主贯通域指标只统计入口—出口主连通分量，排除虚拟节点和零体积节点。

同时记录：

\[
N_{free}(0)-N_{free}(t)
\]

和：

\[
\int_0^t\sum_iJ_{gl,i}\,dt.
\]

前者是自由相库存损失，后者是累计相间传质，二者不能互相替代。

## 18. 沿流向分段

\[
x^*=\frac{x-x_{min}}{x_{max}-x_{min}}.
\]

默认20段。孔按孔心唯一分配，喉按中点分段。每段记录饱和度、库存、界面面积、累计传质、气泡数、事件和左右界面通量。

收缩比：

\[
\chi_{ij}=\frac{\min(r_{p,i},r_{p,j})}{r_{t,ij}}.
\]

异常 `rt>min(rp)` 只报告，不自动截断。

## 19. 终止条件

支持主贯通域饱和度阈值：

\[
\bar S_{g,main}\le S_{g,target},
\]

以及残余自由气体阈值：

\[
N_{free}(t)\le f_{res}N_{free}(0).
\]

还保留最大物理时间、磁盘限制、手动STOP、数值失败和守恒失败。最大物理时间先达到时应报告为右删失。

## 20. 未确认的物理闭合

论文应明确以下内容仍需实验或文献标定：

- 固定 `kL`；
- 球形等效界面面积；
- Qin参数 `b`；
- 有效接触角；
- 孔半径作为毛细有效半径；
- 相对渗透指数 `m`；
- 扩散面积指数 `alpha`；
- tortuosity；
- constant-Z近似；
- 出口回流储层浓度；
- 网络对真实样品的代表性。

25°C物性可以是文献参考值，但没有来源和不确定度时不能称为实验确认值。

## 21. 论文建议报告

- 网络节点、边、孔隙体积和分段定义；
- 初始 `Sw/Sg` 分布；
- 压力、总流量、入口和回流浓度；
- `T, mu, rho, sigma, Hcp, Dm, kL, Z, theta, b, m, alpha, tau`；
- `mass_transfer_multiplier`；
- backward-Euler和Newton双门槛；
- 时间步范围和事件策略；
- 气泡消失阈值；
- 最大摩尔与体积守恒误差；
- restart一致性；
- 自由气体、溶解库存、累计传质和边界输出；
- right-censored算例；
- 未确认参数和敏感性分析。

## 22. 实现文件

```text
src/config.cpp        配置解析和校验
src/network.cpp       网络、cluster和massless junction
src/solver.cpp        通量、残量、Newton和active-set
src/output.cpp        checkpoint、VTK/VTP和审计
src/analysis.cpp      全局、分段和守恒历史
src/main.cpp          时间循环、终止和manifest
src/amgx_backend.cpp  AMGX线性后端
src/pardiso_backend.cpp  oneMKL PARDISO线性后端及固定结构符号分析复用
```
