# 气泡溶解求解器中文源码阅读导引

本文面向已经熟悉 `../trans-react.cpp`、但第一次阅读 `bubble_dissolution_solver` 的研究人员。目标是先看懂物理量、控制方程、自动微分和一次 Newton 迭代，再看工程性代码。

本文只解释当前实现，不提出或执行重构。文中行号对应 2026-09-04 前已经通过完整数值验证的源码版本。

## 1. 先建立一张总图

当前程序每个时间步做的事情可以概括为：

```text
旧状态
  ↓
决定本步有哪些未知量
  ↓
根据当前猜测选择迎风方向和传质开/关分支
  ↓
写出 Rl、Rd、Rg、Rc
  ↓
FADBAD++ 自动得到 Jacobian
  ↓
求解 J·delta = -R
  ↓
限制步长并更新 Pl、C、Pg、Sw
  ↓
同时检查残量和相对增量
  ↓
处理已经消失的气泡并检查守恒
  ↓
接受时间步，写诊断和 checkpoint
```

与 `trans-react.cpp` 最重要的共同点是：都用 `fadbad::B<double>` 直接书写残量，再用 `diff()` 和 `d()` 取得 Jacobian。

最重要的不同点是：当前程序还有动态自由度、zero-distance cluster、massless junction、气泡 active set、变量/残量缩放、线搜索、时间步回退和扩展精度残量。这些不是额外的物理模型，而是为了正确处理当前网络和保证收敛。

## 2. 推荐阅读顺序

### 第一步：`include/bubble/model.hpp`

先看这里建立词汇表，不需要一次看完所有配置字段。

- `model.hpp:142-184`：`Pore`、`Edge`、`ZeroDistanceCluster`、`Network`。
- `model.hpp:186-199`：一个节点的状态 `NodeState`。
- `model.hpp:201-215`：未知量和残量的编号表 `DofMap`。
- `model.hpp:226-262`：一次装配和一次时间步的结果。
- `model.hpp:271-312`：毛管压力、气液界面面积、自由气体摩尔数。
- `model.hpp:314-359`：`Solver` 的主要入口。

第一遍只要记住：状态是 `Pl/Pg/Sw/C`，方程是 `Rl/Rd/Rg/Rc`，`DofMap` 把它们放进线性系统。

### 第二步：`src/solver.cpp:15-92`

这里是最短、最接近公式的部分。

- `solver.cpp:15`：`using AD = fadbad::B<double>`。
- `solver.cpp:28-46`：一条边上的液体流量、对流和扩散。
- `solver.cpp:49-69`：缩放残量和缩放状态的无穷范数。

建议先独立读懂 `flux_values()`，再进入大函数。

### 第三步：`src/network.cpp`

这一步解释为什么当前自由度不能像 `trans-react.cpp` 那样按固定块简单排列。

- `network.cpp:110-185`：读取 pore/edge，建立邻接关系。
- `network.cpp:187-229`：检查体积、长度、边界标志和饱和度。
- `network.cpp:231-319`：用 union-find 建立 zero-distance cluster，并识别 massless junction。

### 第四步：按顺序阅读 `src/solver.cpp` 的五段

1. `Solver::initialize()`，`solver.cpp:100-223`：初始状态和初始小气泡处理。
2. `Solver::make_dof_map()`，`225-297`：给未知量和方程编号。
3. `Solver::assemble()`，`299-553`：构造残量和 Jacobian，是阅读重点。
4. `Solver::attempt_step()`，`596-927`：一次 Newton 时间步。
5. `Solver::advance_with_retry()`，`929-946`：失败后减小 `dt` 再试。

`compute_balance()`（`556-594`）可以放在 `attempt_step()` 之后读。

### 第五步：`src/main.cpp`

先读 `main.cpp:74-122` 和 `166-265`：配置、restart、时间循环、输出时机。

`main.cpp:128-142,188-240` 说明诊断列怎样从 `StepReport` 写入 TSV 和 manifest。

### 第六步：`src/output.cpp`

- `output.cpp:94-134`：VTK。
- `output.cpp:136-241`：退化几何审计。
- `output.cpp:268-365`：扩展精度 checkpoint 的写入和读取。

### 第七步：`src/config.cpp` 和 `CMakeLists.txt`

最后再看严格配置检查和构建规则。与 AD 最相关的只有：

- `CMakeLists.txt:14-16` 查找外部 FADBAD++；
- `CMakeLists.txt:83-87` 把 FADBAD++ include 和 Eigen 连接到核心目标；
- `model.hpp:3` 包含 `badiff.h`。

## 3. 状态变量和自由度归属

| 变量 | 物理意义 | 单位 | 自由度归属 |
|---|---|---|---|
| `Pl` | 液相绝对压力 | Pa | 每个内部 zero-distance cluster 一份；cluster 所有成员共享 |
| `C` | 液相中溶解组分的摩尔浓度 | mol/m3 | 每个内部 cluster 一份；cluster 所有成员共享 |
| `Pg` | 静止气泡的绝对压力 | Pa | 只属于一个仍含气泡且体积为正的成员；不在 cluster 成员间共享 |
| `Sw` | 液相饱和度，气相饱和度为 `1-Sw` | 无量纲 | 只属于一个仍含气泡的成员；无气节点固定为 1 |

边界节点是外部 reservoir。它们提供给定的 `Pl/C`，不分配内部存储自由度，也没有气泡自由度。

massless junction 的体积为零、`Sw=1`。它仍有 `Pl/C`，但没有 `Pg/Sw`；它的方程只表达流入与流出相等。

### 方程与未知量在矩阵中的配对

`make_dof_map()` 在 `solver.cpp:230-252` 使用同一个编号位置配对方程行和未知量列：

```text
Pl 列 ↔ Rl 行
C  列 ↔ Rd 行
Pg 列 ↔ Rg 行
Sw 列 ↔ Rc 行
```

这只是线性系统的排列方式。例如 `Rc` 是毛管闭合方程，不表示 `Sw` 自己有一条独立的时间演化方程。

## 4. 四类残量的物理意义

所有量都在新时间层评价，旧状态只出现在累积项。

### 4.1 `Rl`：液相体积守恒，单位 m3/s

```text
Rl = V·(Sw_new-Sw_old)/dt + 所有向外液体流量之和
```

`Rl=0` 表示液相体积变化与净流出相平衡。

无气节点的 `Sw_new=Sw_old=1`，所以累积项为零，`Rl` 退化为不可压缩液体的流量守恒。

### 4.2 `Rd`：溶解组分守恒，单位 mol/s

```text
Rd = [V·Sw·C 的新旧差]/dt
     + 所有向外的对流和扩散摩尔流率
     - ntr
```

`ntr>0` 表示自由气体进入液相，所以对液相来说是来源，在移到残量左边后写成 `-ntr`。

### 4.3 `Rg`：自由气体摩尔守恒，单位 mol/s

```text
Rg = (ng_new-ng_old)/dt + ntr

ng = Pg·V·(1-Sw)/(Z·R·T)
```

自由气泡不沿网络边移动，因此 `Rg` 没有相邻节点的气相通量。`ntr>0` 会消耗自由气体，所以在自由气体残量中是 `+ntr`。

### 4.4 `Rc`：毛管闭合关系，单位 Pa

```text
Rc = Pg - Pl - Pc(Sw)
```

它把气泡压力、液相压力和饱和度联系起来。

## 5. 四个控制方程在 `solver.cpp` 中的位置

### 公共本构

- `Q/Fadv/Fdiff`：`solver.cpp:28-46`。
- `Pc(Sw)`：`model.hpp:271-299`。
- 气液界面面积：`model.hpp:301-305`。
- 自由气体摩尔数：`model.hpp:307-312`。
- `ntr`：`solver.cpp:331-343`。

### `Rl`

AD 表达式在 `solver.cpp:347-355`：

```cpp
rl = p.volume * (sw - sw_old) / dt;
```

边通量在 `solver.cpp:468-478` 加到 cluster 的 `Rl` 行。对一条按 `a -> b` 定向的边，a 端加 `q`，b 端减 `q`，所以一条内部边在全网求和时自动抵消。

`solver.cpp:345-353` 还有 closed all-liquid component 的压力基准特例。它用一个压力基准方程替代冗余 `Rl`，消除不可压缩压力的任意常数。

### `Rd`

AD 表达式在 `solver.cpp:356-358`：

```cpp
rd = p.volume *
     (sw * (c - c_old) + c_old * (sw - sw_old)) / dt
     - ntr;
```

括号内与 `sw*c - sw_old*c_old` 完全等价，但避免直接相减两个很接近的库存数。

边上的 `Fadv+Fdiff` 在 `solver.cpp:469-478` 按方向加入 `Rd`。

### `Rg`

AD 表达式在 `solver.cpp:359-364`：

```cpp
rg = coefficient *
     ((1.0 - sw) * (pg - before.pg)
      - before.pg * (sw - sw_old)) / dt
     + ntr;
```

它与 `(ng_new-ng_old)/dt+ntr` 等价，也使用增量形式减小相消误差。

### `Rc`

AD 表达式在 `solver.cpp:365`：

```cpp
rc = pg - pl - capillary_pressure(sw, p.radius, config_.physics);
```

### cluster 聚合

一个 cluster 的所有成员都把本地 `Rl/Rd` 累加到相同的 `shared.rl/shared.rd` 行，见 `solver.cpp:316-379`。每个含气成员仍把 `Rg/Rc` 写入自己的 `d.rg/d.rc` 行，见 `381-395`。

零长度边在 `solver.cpp:426-428` 直接跳过；cluster 之间的正长度边仍只计算一次。

## 6. 从 `trans-react.cpp` 看当前 Jacobian 写法

### 6.1 原程序的典型写法

`trans-react.cpp:2873-2901` 的核心模式是：

```cpp
reverse_mode<double> Pi, Wi, F;
F = func_BULK_PHASE_FLOW_kong(Pi, Pjs, Wi, Wjs, i);
F.diff(0, 1);
B[row] = -F.val();
COO_A[diagonal].val = Pi.d(0);
```

含义是：

1. `Pi/Wi/Pjs/Wjs` 是带反向求导记录的数；
2. 普通算术构成残量 `F`；
3. `F.diff(0,1)` 把 `F` 设为第 0 个输出；
4. `Pi.d(0)` 取得 `dF/dPi`；
5. 导数写进 COO 矩阵；
6. AMGX 求解后，`trans-react.cpp:6694-6702` 用 `state += dX` 更新状态。

### 6.2 当前程序的相同核心

当前只是把类型别名缩短为：

```cpp
using AD = fadbad::B<double>;  // solver.cpp:15
```

节点最多同时有 `Rl/Rd/Rg/Rc` 四个输出。`solver.cpp:397-421` 先收集实际存在的方程，再给每个输出一个编号：

```cpp
equations[r].second->diff(r, outputs);
...
const double derivative = variable.second->d(r);
entries.emplace_back(row, column, derivative);
```

因此：

- `diff(r, outputs)` 对应原程序的 `F.diff(0,1)`；
- `variable.d(r)` 对应原程序的 `Pi.d(0)` 或 `Pjs[k].d(0)`；
- `Eigen::Triplet(row,column,value)` 对应原程序的 `COO_A[...]`；
- 当前 `rhs=-residual` 对应原程序的 `B[row]=-F.val()`。

一条边同时产生液体流量 `q` 和组分总通量 `total_f`，所以 `solver.cpp:480-506` 使用两个输出：

```cpp
flux.q.diff(0, 2);
total_f.diff(1, 2);
```

然后 `d(0)` 读取液体流量导数，`d(1)` 读取组分通量导数。

当前没有自写自动微分，也没有用有限差分生成生产 Jacobian。有限差分只存在于测试中，用来检查 FADBAD++ 得到的 Jacobian。

## 7. 一次 Newton 迭代的简化伪代码

下面故意省略容器和文件操作，只保留计算过程：

```text
给定旧时间层状态 old
给定当前猜测 trial

1. 用 trial 的压力判断每条边的流向
2. 用 trial 判断每个气泡是否仍处于“溶解传质开启”分支

3. 对每个节点：
     计算累积项、ntr、Rl、Rd
     若有气泡，再计算 Rg、Rc
     用 B<double> 得到这些节点项对 Pl、C、Pg、Sw 的导数

4. 对每条正长度边：
     只计算一次 Q、Fadv、Fdiff
     把通量以相反符号加入两端残量
     用 B<double> 得到边通量对两端变量的导数

5. 缩放 Jacobian 和残量
6. 解 J·delta = -R
7. 先缩小 delta，保证压力、浓度和饱和度不越界
8. 再尝试不同步长，使新残量充分下降
9. 得到新的 trial

10. 若总缩放残量和相对增量都小于容差：
      Newton 收敛
    否则：
      回到第 1 步
```

若 Newton、线搜索或守恒检查最终失败，`advance_with_retry()` 会恢复旧状态并减小 `dt`，不是放宽容差。

## 8. 三种对象在一次时间步中的数据流

### 8.1 一个不含气泡的内部节点

假设它所在的 cluster 只有它一个成员，并且连到正常网络。

1. **旧状态：** 有 `Pl/C`；`Sw=1`；`active_gas=false`；`Pg` 不参与计算。
2. **自由度：** `make_dof_map()` 只分配 `Pl` 和 `C` 两列，以及 `Rl` 和 `Rd` 两行；`solver.cpp:230-246`。
3. **冻结分支：** 它没有传质分支，`ntr=0`；边的迎风方向仍由压力决定。
4. **节点残量：** `Rl` 的饱和度累积为零；`Rd` 变成 `V(C_new-C_old)/dt`。
5. **边残量：** 相邻边把液体流量加入 `Rl`，把对流和扩散加入 `Rd`。
6. **Jacobian：** 节点和边的 AD 图只向实际存在的 `Pl/C` 列写导数。
7. **Newton 更新：** `solver.cpp:754-755` 更新 `Pl/C`；不会更新 `Pg/Sw`。
8. **接受后：** 状态仍无气泡，输出中 `Sw=1`、`active_gas=0`。

如果它处于完全封闭且全液体的连通分量，系统会选择一个节点的 `Rl` 行作为压力基准，避免所有压力同时加常数而不改变流量。

### 8.2 一个含气泡的内部节点

1. **旧状态：** 有 `Pl/C/Pg/Sw`，并且 `active_gas=true`。
2. **自由度：** cluster 共享 `Pl/C`；该气泡成员独立获得 `Pg/Sw`，对应 `Rg/Rc`；`solver.cpp:230-252`。
3. **冻结分支：** 根据 `Hcp*Pg-C` 判断本次 Newton 线性化是否计算 `ntr`；`solver.cpp:629-634`。
4. **节点残量：** 同时构造 `Rl/Rd/Rg/Rc`；`solver.cpp:326-395`。
5. **边残量：** 液体饱和度影响边导流率和有效扩散面积；该节点没有气相边通量。
6. **Jacobian：** `Pg/Sw` 会影响 `ntr`、气体库存、毛管压力和相邻液相边通量。FADBAD++ 自动把这些依赖写入对应列。
7. **Newton 更新：** `solver.cpp:754-757` 更新四个变量；fraction-to-boundary 保证 `Pg>0`、`Sw` 保持在允许范围内。
8. **收敛检查：** 残量和相对增量必须同时通过。
9. **气泡消失检查：** 若 `1-Sw` 或 `ng` 低于阈值，`solver.cpp:866-905` 把剩余自由气体保守地转入溶解库存，设 `Sw=1`，并永久删除以后的 `Pg/Sw` 自由度。
10. **接受后：** 守恒检查通过才提交新状态。

### 8.3 一个 zero-distance cluster

假设几个原 PNM 控制节点被零长度边连接。

1. **网络预处理：** `network.cpp:240-270` 用 union-find 把它们放入同一 cluster；原节点和零长度边仍保留用于映射和输出。
2. **共享状态：** cluster 所有成员共享一个 `Pl` 和一个 `C`。
3. **独立气泡：** 每个有正体积且含气泡的成员仍有自己的 `Pg/Sw`。不同成员的气泡体积、半径、界面面积不会混合。
4. **共享方程：** 整个 cluster 只有一行 `Rl` 和一行 `Rd`。每个成员的累积项和 `ntr` 都累加到这两行。
5. **零长度边：** 不计算 `1/L` 通量；`solver.cpp:426-428` 跳过它们。
6. **外部正长度边：** 仍在具体成员上计算其 `Sw` 对导流率和扩散面积的影响，但通量最终进入共享 `Rl/Rd`。
7. **成员气泡方程：** 每个气泡成员各有一对 `Rg/Rc`。
8. **Newton 更新：** cluster 成员的 `d.pl/d.c` 指向同一个 delta 编号，因此所有成员得到相同的 `Pl/C` 更新；各气泡的 `Pg/Sw` 更新彼此独立。
9. **输出：** VTK 仍逐个输出原节点，同时用 `zero_distance_cluster` 字段说明共享关系。

cluster 中若包含 `V=0,Sw=1` 的 massless junction，它对累积项贡献严格为零，但仍通过边通量参加共享守恒。

## 9. 为什么有 AD-double 和 long-double 两条数值路径

这是理解 `assemble()` 的关键。

### 9.1 AD-double 路径：负责 Jacobian

当前 AD 类型明确是：

```cpp
using AD = fadbad::B<double>;
```

`solver.cpp:326-365` 用 `AD` 构造 `Rl/Rd/Rg/Rc`，`397-421` 和 `480-506` 用 `diff()/d()` 取得偏导。

这条路径的职责是回答：

```text
当前状态发生一个很小变化时，各残量怎样变化？
```

所得 Jacobian 是 `double` 稀疏矩阵，交给 Eigen 求解。当前验证表明 double Jacobian 足以得到可靠 Newton 方向。

### 9.2 long-double primal 路径：负责判断真实残量大小

完整网络含有高压力、非常接近全液状态的气泡。某些有效更新比一个 binary64 状态的最小间隔还小。若状态和关键累积差只用 double，两个接近的大数相减后可能看不见仍然存在的小残量。

因此：

- `NodeState` 使用 `long double`，见 `model.hpp:186-199`；
- `solver.cpp:310` 创建 `residual_extended`；
- `solver.cpp:332-395` 用 long double 重算关键节点累积、传质、`Rg` 和 `Rc`；
- 边通量仍从 AD/double primal 读取，再累加到 long-double residual 容器，见 `461-479`；
- 最后在 `511-513` 转为 Eigen 的 double residual，用于缩放范数和线性代数。

这条路径的职责是回答：

```text
当前状态代入方程后，残量本身到底还有多大？
```

两条路径不是两套自动微分：

- AD-double：只负责导数/Jacobian；
- long-double primal：只负责高精度状态和残量数值；
- 没有 long-double Jacobian，也没有第二个 AD 引擎。

## 10. `ntr = 0` 为什么不是把物理传质率清零

`solver.cpp:397-413` 先用同一个 `ntr` 构造 `Rd` 和 `Rg`，再对所有输出调用 `diff()`：

```cpp
equations[r].second->diff(r, outputs);
ntr = 0.0;
```

FADBAD++ 的 `B<double>` 使用带引用计数的反向计算图。`Rd`、`Rg` 和局部变量可能共同保存对同一图节点的引用。输出已经播种后，需要释放多余的中间引用，导数才能完整传播到 `pl/c/pg/sw` 这些叶变量。

这里的 `ntr = 0.0` 只做一件事：

```text
让局部 AD 变量 ntr 不再持有旧计算图。
```

它不做以下事情：

- 不修改已经构造好的 `rd` 和 `rg`；
- 不修改 `out.transfer[i]`，该值已在 `solver.cpp:380` 保存；
- 不把物理传质率改为零；
- 不改变 long-double residual 中已经加入的 `ntr_value`。

边通量处的：

```cpp
flux.f_adv = 0.0;
flux.f_diff = 0.0;
```

在 `solver.cpp:483-486` 作用相同：释放 `total_f` 已经引用的中间 AD 图，随后才读取叶变量的 `d(1)`。这些赋零语句位于局部 Jacobian 提取阶段，不是物理裁剪。

## 11. 初次阅读时可以暂时跳过什么

为了先看懂物理主线，第一次阅读可以暂时跳过：

- `src/config.cpp` 的 JSON 字段逐项读取和错误文字；只要知道所有物性、尺度和容差都经过严格检查。
- `src/output.cpp:41-90` 的 SHA-256 实现。
- `src/output.cpp:94-266` 的 VTK 文本格式和 VTK 自检细节。
- `src/output.cpp:268-407` 的 checkpoint JSON、provenance hash 和 UTC 时间格式；先知道 `long double` 状态会以字符串保存即可。
- `src/main.cpp:276-432` 的成功/失败 manifest 组装。
- `src/main.cpp:128-159` 的 TSV 表头和初始 regularization 文件格式。
- `solver.cpp:781-845` 的长错误诊断字符串；先看线搜索判据 `740-780`。
- `solver.cpp:514-546` 的四类最大残量节点统计；它不改变方程或 Newton 更新。
- `solver.cpp:949-998` 的 restart 状态合法性检查。
- `tests/` 和 `scripts/` 的具体实现；读完主线后再用它们确认 Jacobian、cluster、massless junction、守恒和输出行为。

不要在第一遍跳过：`flux_values()`、`make_dof_map()`、`assemble()` 中的四类残量、`diff()/d()`、Newton 双收敛门，以及 active-set 的气泡消失处理。

## 12. 从控制方程到输出诊断的索引表

| 控制方程 | 相关函数/本构 | 残量写入 | Jacobian 取得 | Newton 更新变量 | 输出诊断 |
|---|---|---|---|---|---|
| `Rl` 液相体积 | `flux_values()`；`Solver::assemble()` | 节点项 `solver.cpp:347-355,368-378`；边项 `468-478` | 节点 `397-421`；边 `480-506` | `Pl`：`solver.cpp:754` | `Rl_scaled_inf/Rl_max_node_id`：`solver.cpp:514-530`，写出于 `main.cpp:196-198` |
| `Rd` 溶解组分 | `flux_values()`、`interface_area()`、`ntr` | 节点项 `solver.cpp:358,374-379`；边项 `469-478` | 节点 `397-421`；边 `480-506` | `C`：`solver.cpp:755` | `Rd_scaled_inf/Rd_max_node_id`：`solver.cpp:514-530`，写出于 `main.cpp:196-198` |
| `Rg` 自由气体 | `free_gas_moles()`、`interface_area()`、`ntr` | `solver.cpp:359-364,381-392` | 节点 `397-421` | `Pg`：`solver.cpp:756` | `Rg_scaled_inf/Rg_max_node_id`：`solver.cpp:532-546`，写出于 `main.cpp:196-198` |
| `Rc` 毛管闭合 | `capillary_pressure()` | `solver.cpp:365,393-395` | 节点 `397-421` | `Sw`：`solver.cpp:757` | `Rc_scaled_inf/Rc_max_node_id`：`solver.cpp:532-546`，写出于 `main.cpp:196-198` |

总残量、状态尺度、更新尺度、相对增量和两道收敛门在 `solver.cpp:635-655` 计算，并由 `main.cpp:188-219` 写入 `diagnostics.tsv` 和 `run_manifest.json`。

## 13. 建议的第一次实际阅读路线

如果只有半小时，可以按下面顺序打开代码：

```text
model.hpp:186-215          状态和 DOF
model.hpp:271-312          Pc、Agl、ng
solver.cpp:28-46           Q、Fadv、Fdiff
solver.cpp:225-297         谁拥有哪些未知量
solver.cpp:326-395         四类残量
solver.cpp:397-421         节点 Jacobian
solver.cpp:426-506         边残量和边 Jacobian
solver.cpp:635-780         Newton、双门和线搜索
solver.cpp:847-925         提交状态、气泡消失和守恒
network.cpp:231-319        回头理解 cluster 来源
main.cpp:166-265           完整时间循环和输出
```

读完这条路线后，再看 checkpoint、manifest、严格输入校验和全部测试，会容易很多。
