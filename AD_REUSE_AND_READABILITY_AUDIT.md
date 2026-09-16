# 原示例风格与自动微分复用审查

审查时间：2026-09-03T18:52:40Z（UTC）
审查性质：只读。本文没有伴随任何 C++、配置或既有输出改动。
审查对象：`../trans-react.cpp`、FADBAD++ 2.1 的 `fadbad.h/badiff.h/fadiff.h`，以及当前项目的 `CMakeLists.txt`、`include/bubble/model.hpp` 和全部 `src/*.cpp`。测试源码用于确认 AD 的验证边界，但不作为生产求导路径。

## 1. 结论摘要

1. 当前程序已经直接复用原示例所用的 FADBAD++ 反向自动微分类型和调用协议：`fadbad::B<double>`、重载算术/数学运算、在残量输出上调用 `diff()`、再从输入变量调用 `d()` 取得导数。它不是对 FADBAD++ 的仿写。
2. 当前程序没有另外实现一套自动微分、符号求导或生产用有限差分 Jacobian。生产 Jacobian 的唯一导数来源是 FADBAD++ `B<double>`。
3. 当前有一份与 AD 残量并行的 `long double` primal 残量计算，用于解决已实际观测到的 binary64 消减/状态 ULP 停滞。它重复了部分残量公式，但不计算导数，不属于第二套 AD。此重复是主要的可维护性风险，也是后续最值得小心整理的位置。
4. 与 `trans-react.cpp` 相比，当前物理表达式更短、数据所有权更清楚；复杂度主要集中在局部多输出反向图、zero-distance cluster 聚合、massless junction、active set、扩展精度 primal、缩放稀疏 Newton 和接受/回退逻辑。这些复杂性有数值或拓扑理由，不能为追求“示例风格”而删除。
5. 没有发现数值正确性所要求的立即源码修改。本审查建议的改动均应等用户确定范围后再实施，并以本文件末尾的 SHA-256 基准和已经通过的 13/13 Release、11/11 Sanitizer 为回归基线。

## 2. 原 `trans-react.cpp` 的自动微分和 Newton 数据流

### 2.1 头文件、核心类型和变量声明

- `trans-react.cpp:53-54` 同时包含 `badiff.h`（反向 AD）与 `fadiff.h`（前向 AD）。
- `trans-react.cpp:60-62` 使用 `fadbad` 命名空间并声明：

```cpp
template <typename T>
using reverse_mode = B<T>;
```

- 实际残量代码使用 `reverse_mode<double>`，例如 `kong_matrix_QIN()` 在 `trans-react.cpp:2873-2877` 声明本节点变量、残量和邻居数组：

```cpp
reverse_mode<double> Pi, Wi, F;
reverse_mode<double>* Pjs, *Wjs;
```

- 对传输方程，`trans-react.cpp:3038-3049` 为两个组分分别声明 `P1i/P2i/C1i/C2i/F1/F2` 及邻居变量。
- 全文件没有 `fadbad::F<...>`、`FTypeName` 或等价的前向 AD 实例；`fadiff.h` 只是被包含，没有参与实际 Jacobian。

FADBAD++ 的相关接口来自原始头文件：

- `fadbad.h:78-85` 用宏把 backward/forward/Taylor 类型名定义为 `B/F/T`。
- `badiff.h:255-318` 定义 `B<T>` 的值、反向图引用、`val()/x()`、`diff(idx,size)` 和 `d(i)` 接口。
- `badiff.h:204-231,278-286` 表明反向传播与引用释放相关：输出播种后，图节点释放时向叶变量传播导数。因此保留中间 AD 对象会影响读取导数的时机。
- `fadiff.h:38-110` 定义前向 `F<T,N>`；原示例和当前生产实现都没有使用它。

### 2.2 原示例的物理残量构造

原示例不是本气泡模型；它处理压力、CO2/CH4 输运、吸附等另一套方程。可复用的是“用 AD 标量直接书写物理残量”的风格，而不是具体公式。

- 压力/体积流动代表：`func_BULK_PHASE_FLOW_kong()`，`trans-react.cpp:2044-2113`。它循环邻接喉道，形成类似 `conductivity * (Pi-Pj)` 的通量和。
- 组分输运代表：`func_TRANSPORT_FLOW_kong()`，`trans-react.cpp:2296-2334`。其数据流直接表现为“时间累积 + 上游对流 + 扩散”：

```cpp
Return += Pb[Pore_id].volume * (Wi - Pb[Pore_id].C1_old) / dt;
...
Return += Wi * Conductivity * (Pi - Pjs[iCounter]);
Return += Conductivity_dis_co2 * (Wi - Wjs[iCounter]);
```

- 上游分支使用当前普通 `double` 压力决定，然后在选定分支内构造 AD 表达式；这是“冻结非光滑分支、对分支内连续公式求导”的早期形式。

### 2.3 原示例的 Jacobian 获取

`kong_matrix_QIN()`（`trans-react.cpp:2842-约3299`）展示了完整模式：

1. 把当前普通状态赋给 `reverse_mode<double>` 本节点和邻居变量（`2873-2885`）。
2. 调用物理残量函数得到标量 `F`（`2895`）。
3. 在残量输出上调用 `F.diff(0,1)`（`2896`），即播种一个输出分量。
4. 用 `-F.val()` 形成 Newton 右端（`2897`）。
5. 用 `Pi.d(0)`、`Pjs[k].d(0)` 取得对本节点和邻居未知量的偏导并直接写入裸 `COO_A`（`2899-2937`）。

传输方程对 `F1/F2` 各自调用 `diff(0,1)`（`trans-react.cpp:3062-3068`），再分别从压力和浓度叶变量读取 `d(0)`（`3070-3137`）。因此原示例的基本粒度是“一个标量残量对应一张局部反向图”。

COO 随后由 `Matrix_COO2CSR()` 排序并转换成 CSR（`trans-react.cpp:7139-约7204`）。这部分是手工稀疏存储，不是 AD 本身。

### 2.4 原示例的 Newton 求解和更新

主执行分支由 `main()` 在 `trans-react.cpp:7237-7261` 选择；Neumann transport 分支进入 `AMGX_solver_C_kong_PNM_Neumann_boundary()`。

- `trans-react.cpp:4586-4593` 每次内迭代重新构造 AD Jacobian、转换 CSR、调用 AMGX，并以 `norm_inf > eps` 控制循环。
- `AMGXsolver_subroutine_kong()` 在 `trans-react.cpp:6667-6680` 上传/替换矩阵系数，求解 `J dX = -F` 并下载 `dX`。
- `trans-react.cpp:6683-6690` 名为 `norm_inf` 的量实际是 `sqrt(sum(dX[i]^2))`，即更新二范数，不是无穷范数，也不是残量范数。
- `trans-react.cpp:6692-6722` 对压力和两个浓度执行无阻尼更新 `state += dX`。

原示例因此提供了清晰的 AD 教学骨架，但没有当前程序需要的未知量/残量缩放、fraction-to-boundary、Armijo、残量与更新双门、状态回滚、`dt` 重试、active-set 后处理和守恒接受门。不能把原 Newton 控制流原样移植到当前程序。

## 3. 当前程序实际复用的原框架组件

这里的“复用”分两层：FADBAD++ 库/API 是直接复用；`trans-react.cpp` 的具体类、函数、全局数组和物理方程没有被复制或链接。

| 复用项                   | 当前位置                                                             | 实际作用                                                            | 与原示例的关系                                 |
| ------------------------ | -------------------------------------------------------------------- | ------------------------------------------------------------------- | ---------------------------------------------- |
| 原始 FADBAD++ 2.1 头文件 | `CMakeLists.txt:6,14-16,83-86`                                     | 查找外部 `badiff.h` 并作为 SYSTEM include 暴露给目标              | 同一外部库，不 vendoring                       |
| 反向 AD 头               | `include/bubble/model.hpp:3`                                       | 包含 `badiff.h`                                                   | 原示例为 `trans-react.cpp:53`                |
| 反向标量类型             | `src/solver.cpp:15`：`using AD = fadbad::B<double>`              | 节点/边局部残量图的标量                                             | 等价于原 `reverse_mode<double>`              |
| FADBAD 算术与数学运算    | `include/bubble/model.hpp:271-312`；`src/solver.cpp:28-46`       | `Pc(Sw)`、界面面积、EOS、导流/对流/扩散共同支持普通标量和 AD 标量 | 复用 `B<T>` 运算符及 `exp/pow/cos`         |
| primal 值读取            | `src/solver.cpp:380,464-466`；`include/bubble/model.hpp:264-269` | 从 AD 图取得传质量和边通量数值                                      | 原示例主要使用 `.val()`；当前多使用 `.x()` |
| 输出播种                 | `src/solver.cpp:397-413`、`480-486`                              | 节点最多 4 个输出、边 2 个输出调用 `diff(component,count)`        | 同一 `diff()` 协议，但当前利用多输出分量     |
| 叶变量导数读取           | `src/solver.cpp:414-421`、`487-506`                              | 读取 `d(r)` 并生成 Jacobian triplet                               | 同一 `d()` 协议                              |
| 反向图释放语义           | `src/solver.cpp:410-413,483-486`                                   | 在读叶导数前释放共享中间量，触发 FADBAD++ 传播                      | 源于 `badiff.h` 的引用计数传播机制           |
| AD 冒烟验证              | `tests/test_ad_smoke.cpp:25-43`                                    | 把 FADBAD 结果与解析/中心差分比较                                   | 直接验证 `B/diff/d` API                      |

当前还沿用了原示例的两种高层思想，但不是代码级复用：

- 局部 stencil 构造 AD 图，再把非零导数放入全局稀疏矩阵。
- 求解 Newton 线性系统 `J delta = -R` 后更新状态。

当前稀疏容器是 `std::vector<Eigen::Triplet<double>>` 和 `Eigen::SparseMatrix`（`src/solver.cpp:313-314,547-551`），不是原示例的裸 `Acoo/ia/ja/a`。当前线性后端是 Eigen SparseLU 或 BiCGSTAB+ILUT（`src/solver.cpp:659-709`），没有复用 AMGX/CUDA 生命周期。

## 4. 是否存在重复的 AD 或残量求导机制

### 4.1 没有第二套求导机制

生产代码中没有以下任何机制：

- 自建 dual number、tape 或 expression-template AD；
- 手写 `dRl/dU`、`dRd/dU`、`dRg/dU`、`dRc/dU` 解析 Jacobian；
- 用有限差分生成生产 Jacobian；
- 使用 `fadbad::F` 前向模式与 `fadbad::B` 重复求导。

`tests/test_bubble_solver.cpp:144-约177,306-371` 的中心差分仅用于审计 FADBAD Jacobian，不进入生产求解。

### 4.2 存在有意的 primal 残量双路径

`Solver::assemble()` 同时构造：

- `AD rl/rd/rg/rc`（`src/solver.cpp:326-366`），只负责 Jacobian 图；
- `std::vector<long double> residual_extended` 及相应 `*_value`（`310,332-343,368-395`），负责最终残量数值；
- 边项目前从 AD 通量取 primal 后累加到 extended 容器（`461-479`）。

这意味着节点累积、传质、`Rg` 和 `Rc` 公式存在两份表达。原因不是架构偏好，而是完整网络曾在高压、近全液状态下出现小于 binary64 状态 ULP 的有效 Newton 修正；`NodeState` 因此也使用 `long double`（`include/bubble/model.hpp:186-199`）。

审查判断：这不是重复 AD，不能简单删除；但两份公式可能随未来物理修改而漂移。若后续授权整理，应该优先用少量、纯标量的小函数共享“单个物理项”的代数表达，而不是再建立一套抽象框架。

## 5. Rl、Rd、Rg、Rc 到当前代码的逐类映射

控制方程依据为 `../remote_codex_bubble_case/02_控制方程与数值规格.md:93-134`（规格 §5.1-§5.4；Word 原式组 Eq. (29)-(35)）。cluster 扩展依据为补充规格 `06_退化几何处理补充指令.md:45-74`。

| 方程   | 数学含义和单位                                                 | AD 构造                    | extended-primal 构造 | 通量/聚合/DOF                                                                                                                               |
| ------ | -------------------------------------------------------------- | -------------------------- | -------------------- | ------------------------------------------------------------------------------------------------------------------------------------------- |
| `Rl` | `V(Sw_new-Sw_old)/dt + sum Q_out = 0`；`m3/s`              | `src/solver.cpp:347-355` | `368-378`          | cluster 共享行由 `make_dof_map():230-245` 分配；边通量符号累加 `468-478`；closed liquid component 的 gauge 替代式见 `255-295,345-353` |
| `Rd` | `[V Sw C]_new-old/dt + sum(Fadv+Fdiff) - ntr = 0`；`mol/s` | `src/solver.cpp:358`     | `374-379`          | `ntr` 在 `331-343`；边组分通量在 `470-478`；cluster 成员自然加到共享 `rd` 行                                                        |
| `Rg` | `(ng_new-ng_old)/dt + ntr = 0`，无气相边通量；`mol/s`      | `src/solver.cpp:359-364` | `382-392`          | 仅 active gas member 有独立 `rg` DOF：`246-252`                                                                                         |
| `Rc` | `Pg-Pl-Pc(Sw)=0`；`Pa`                                     | `src/solver.cpp:359-365` | `393-395`          | `Pc` 定义在 `include/bubble/model.hpp:271-299`；仅 active gas member 有独立 `rc` 行                                                   |

公共本构和边项：

- `Q/Fadv/Fdiff`：`src/solver.cpp:28-46`。边只计算一次，并向 a/b 两端以相反符号装配（`426-479`）。
- `Agl`、EOS `ng`、Qin `Pc`：`include/bubble/model.hpp:271-312`。
- dissolution-only `ntr`：`src/solver.cpp:331-343`；分支在每次 Newton 迭代前冻结于 `629-634`。
- 最终 residual vector 从 extended 容器转换为 binary64 供缩放稀疏线性代数使用：`511-513`。
- 四类缩放残量和对应输入 PNM 节点 ID：`514-546`。

massless junction 没有另设特殊残量：其 `p.volume==0` 使 `Rl/Rd` 累积项自然为零，剩余外部边通量守恒决定 `Pl/C`。其合法性检查和分类在 `src/network.cpp:187-228,231-319`。

## 6. 原示例与当前实现的代表性写法对照

### 6.1 标量残量与单输出反向图

原示例（`trans-react.cpp:2895-2911`）：

```cpp
F = func_BULK_PHASE_FLOW_kong(Pi, Pjs, Wi, Wjs, i);
F.diff(0, 1);
B[i - inlet] = -F.val();
COO_A[i - inlet].val = Pi.d(0);
...
COO_A[index].val = Pjs[counter].d(0);
```

数据流是：一个节点标量残量 `F` -> 一个反向分量 -> 手算全局下标 -> 写裸 COO。

当前节点实现（`src/solver.cpp:397-421`）：

```cpp
std::vector<std::pair<int, AD*>> equations;
equations.emplace_back(shared.rl, &rl);
equations.emplace_back(shared.rd, &rd);
...
for (unsigned int r = 0; r < outputs; ++r)
  equations[r].second->diff(r, outputs);
ntr = 0.0;
...
entries.emplace_back(equations[r].first, variable.first, derivative);
```

数据流是：一个节点可有 `Rl/Rd[/Rg/Rc]` 多输出共享图 -> 每个输出占一个反向分量 -> 通过 DOF map 生成 Eigen triplet。`ntr=0.0` 不是物理清零，而是释放保留的共享图引用以完成反向传播。

### 6.2 边通量

原示例在 `func_TRANSPORT_FLOW_kong()` 中把邻接循环、上游选择、导流率调用和残量累加混在一起（`trans-react.cpp:2309-2331`）。

当前先用一个短的 scalar-generic 函数表达单边物理（`src/solver.cpp:28-46`）：

```cpp
const Scalar q = conductance * (pla - plb);
const Scalar cup = upwind_a ? ca : cb;
const Scalar adv = q * cup;
const Scalar diff = diffuse ? ... * (ca - cb) / e.length : Scalar(0.0);
return {q, adv, diff};
```

随后 `assemble()` 负责边界值、cluster 两端、符号和稀疏下标（`426-507`）。这种物理/拓扑分离比原示例更适合当前退化网络，应保留。

### 6.3 Newton 更新

原示例直接 `state += dX`（`trans-react.cpp:6692-6703`），用更新二范数控制循环（`6683-6690,4586-4593`）。

当前先做行列缩放并求解（`src/solver.cpp:659-715`），再计算 saturation/positivity 的 fraction-to-boundary（`716-739`），执行 Armijo 回溯并只提交候选状态（`740-851`）。收敛同时要求残量门和相对更新门（`643-655`）；整个步失败后由 `advance_with_retry()` 回滚并缩小 `dt`（`929-946`）。这是数值算法差异，不是单纯语法冗长。

## 7. 哪些可以更接近原示例，哪些必须保留

### 7.1 可以更接近原示例的部分

1. **让方程本体先出现。** 当前 `Rl/Rd/Rg/Rc` 已经相对集中在 `src/solver.cpp:347-366`，但 gauge、extended-primal 和 AD 生命周期细节包围了公式。可把四类“单成员局部项”移入同一小节或几个简单 helper，使阅读顺序重新成为变量 -> 本构 -> 残量 -> `diff/d`。
2. **使用具名物理项。** 例如 `liquid_accumulation`、`dissolved_accumulation`、`gas_accumulation`、`transfer_rate`，再写 `rd = dissolved_accumulation - transfer_rate`。这样既保留稳定增量式，又能一眼对应控制方程。
3. **统一 AD API 读法。** 原示例用 `.val()`，当前常用 `.x()`，`scalar_primal(const B&)` 甚至先复制再 `.x()`（`model.hpp:266-269`）。FADBAD `B` 已公开 `val() const`；未来可评估统一成 `.val()`，减少“x 是状态还是值”的认知成本。
4. **沿用一个明确类型名。** `AD` 很短，但文档上下文不足时不如 `Reverse` 或原示例的 `reverse_mode` 清晰。若改名，只做一个 `using`，不要引入 wrapper class。
5. **把 `diff/d` 的局部模式写成紧邻残量的短块。** 可以借鉴原示例的顺序，但继续使用 DOF map 和 Eigen triplet，不回退到手算数组偏移。

### 7.2 必须保留的复杂性

- **zero-distance cluster：** union-find 构造与校验 `src/network.cpp:56-77,231-319`；共享 `Pl/C/Rl/Rd`、成员独立 `Pg/Sw/Rg/Rc` 的 DOF map `src/solver.cpp:225-297`；cluster 聚合残量 `316-424`。
- **massless junction：** `V==0,Sw==1` 的合法代数节点及零累积行为；不能虚构体积或删掉守恒节点。
- **active set：** 初始化时的不可逆退休和守恒转移 `src/solver.cpp:146-209`，Newton 内冻结传质分支 `629-634`，收敛后的气泡消失与账本 `866-915`。
- **扩展精度：** `NodeState` 的 `long double`（`model.hpp:186-199`）和 extended-primal 残量（`solver.cpp:310,332-395`）。在找到同等严格且通过完整网络基准的替代前不能删除或降回 double。
- **FADBAD 共享图释放顺序：** `solver.cpp:410-413,483-486`。这由所用 FADBAD++ 版本的引用计数反向传播语义决定。
- **稀疏 Newton：** 动态 DOF、局部 triplet、行列缩放、SparseLU/BiCGSTAB+ILUT、线性残量审计、fraction-to-boundary、Armijo、双收敛门和 `dt` retry。
- **守恒和审计：** equation/adjusted balance、残量类别最大值、原 PNM ID、checkpoint 扩展精度状态及输出清单。

## 8. 最小化可读性改进方案

以下只是建议次序，不是本阶段实施内容。

### 第 1 步：仅整理命名和注释，零控制流变化

- 在 `Rl/Rd/Rg/Rc` 定义前标注 `Word Eq. (29)-(35) / spec §5.1-§5.4`，并逐项写 SI 单位：`m3/s, mol/s, mol/s, Pa`。
- 把 `rl/rd/rg/rc` 周围的中间量命名为 `liquid_accumulation`、`dissolved_accumulation`、`gas_accumulation`、`capillary_closure`。
- 明确注释稳定增量恒等式，例如 `Sw*(C-C_old)+C_old*(Sw-Sw_old) == Sw*C-Sw_old*C_old`。
- 明确 `ntr=0` 和 `f_adv/f_diff=0` 是 FADBAD 图生命周期操作，不是物理值修正；最好用一个小作用域让析构自然发生，避免“赋零”的误读。
- 评估把 AD primal 读取统一为 `.val()`；这是原 FADBAD API 直接复用，不引入封装。

### 第 2 步：拆分 `Solver::assemble()`，仍保持简单函数

建议只拆成 3 个内部函数，不建立类层次或新模板体系：

1. `assemble_node_terms(...)`：成员累积、`ntr`、`Rg/Rc` 以及节点局部 Jacobian；
2. `assemble_edge_terms(...)`：一次边通量、两端符号和边局部 Jacobian；
3. `finalize_residual_audit(...)`：extended -> double、四类最大值和 triplet finalize。

现有 `Solver::assemble()` 为 `src/solver.cpp:299-553`，拆分后顶层应只展示“节点项 -> 边项 -> finalize”三步。helper 参数应使用现有 `Network/Config/DofMap/Assembly`，不要增加 visitor、policy、expression tree 或宏。

### 第 3 步：拆分 `Solver::attempt_step()` 的数值控制流

现有函数为 `src/solver.cpp:596-927`。可按已有阶段抽出：

- `freeze_nonsmooth_branches()`（当前 `611-634`）；
- `solve_scaled_linear_system()`（`659-715`）；
- `choose_step_bound_and_line_search()`（`716-846`）；
- `retire_disappeared_bubbles()`（`866-915`）。

顶层仍保留 Newton 顺序和所有失败出口。不要改变容差、迭代次数语义或 `dt` retry。

### 第 4 步：减少 AD/long-double 公式漂移风险

只对纯代数“项”使用已有风格的 `template<class Scalar>` 小函数，例如 gas accumulation 或 dissolved accumulation，使 `AD` 和 `long double` 调用同一表达式。不要把整个 assembly 模板化，不要创建新的 AD facade，也不要强迫稀疏索引、active set 或输出审计进入模板。

这个步骤风险高于前三步，必须先为每个 helper 增加 double/long-double/AD 一致性测试，再运行完整网络收敛审计。

### 第 5 步：最后再整理驱动函数

`main()` 为 `src/main.cpp:74-432`，同时负责 CLI、restart、时间推进、TSV、manifest 和失败 manifest。可选地抽出 `RunAuditAccumulator` 或几个普通函数，但这与 AD 可读性无直接关系，优先级低于 `assemble()` 和 `attempt_step()`。不建议趁 AD 重构同时改变输出格式。

## 9. 修改优先级判定

### 必须修改

**当前为 0 项。** 没有发现重复 AD、错误 Jacobian 来源或与已验证控制方程不一致、因而必须立即修改的生产源码。数值任务刚通过严格基准，此时以“风格统一”为由立即动核心求解器得不偿失。

如果用户授权下一阶段重构，则以下是实施过程的强制约束，而不是本审查发现的缺陷：

- 每个提交保持方程、单位、符号、active-set、cluster/massless 和 extended-precision 语义不变；
- 不改变 `nonlinear.residual_tolerance=1e-8` 或接受逻辑；
- 每个小步对照本文件 SHA-256 基准审查 diff，并重跑相关 AD/Jacobian 测试；最终重跑完整 13/13、11/11 和 5 步全网审计。

### 可选改善

1. 增加方程编号、SI 单位和稳定恒等式注释。
2. 将 `AD` 改为更自解释的单一 alias，并统一 `.val()` 读取风格。
3. 用作用域代替“给中间 AD 对象赋零”来表达 FADBAD 图释放。
4. 按 §8 拆分 `assemble()` 和 `attempt_step()`。
5. 用少量现有风格的 scalar helper 共享 AD/long-double 物理项。
6. 最后才拆分 `main()` 的输出/manifest 逻辑。

### 不建议修改

1. 不建议复制 `trans-react.cpp` 的 `PNMsolver` 大类、全局状态、裸 `new[]/delete[]`、手算 COO 偏移或 CUDA/AMGX 生命周期。
2. 不建议为了“复用更多”重新包含或使用 `fadiff.h/F<double>`；当前局部残量的反向模式已经正确，混入前向模式只会形成第二条求导路径。
3. 不建议自建 AD wrapper、宏 DSL、方程注册表、visitor/policy 层或新的 expression-template 系统。
4. 不建议把 current multi-output reverse seeding 降回“每个残量重建一张图”；那会增加重复计算和代码量。
5. 不建议移除 extended-primal 残量、把 `NodeState` 降回 double，或把 Jacobian 全部升为 long double；前两者会重现已解决的收敛停滞，后者与当前 FADBAD/Eigen 路径和验证基准不匹配。
6. 不建议为了缩短函数删除 cluster、massless junction、active-set、守恒门、双收敛门、线搜索或 timestep retry。
7. 不建议在同一阶段同时改物理公式、稀疏布局、输出格式和可读性；应把纯整理与算法变化严格分开。

## 10. 已阅读的当前核心文件及职责

| 文件                         | 核心职责                                                         | 与 AD/可读性的关系                             |
| ---------------------------- | ---------------------------------------------------------------- | ---------------------------------------------- |
| `include/bubble/model.hpp` | 配置/网络/状态/DOF/报告类型，本构 scalar helper，Solver 接口     | FADBAD include、extended state、可复用本构表达 |
| `src/config.cpp`           | 严格 JSON 读取和前置验证                                         | 无 AD；保证尺度/容差/物性有效                  |
| `src/network.cpp`          | PNM 读取、拓扑验证、union-find cluster、massless 分类            | 解释动态 DOF 与聚合残量复杂性                  |
| `src/solver.cpp`           | 本构边通量、初始化、DOF、残量/Jacobian、Newton、守恒、active set | AD 审查主体                                    |
| `src/output.cpp`           | VTK、几何审计、checkpoint、SHA-256                               | 无求导；保留 extended state 和可复现性         |
| `src/main.cpp`             | CLI、restart、时间推进、诊断、manifest、失败报告                 | 无求导；控制流较长但 AD 重构优先级较低         |
| `CMakeLists.txt`           | 外部依赖、构建、源码树哈希、测试注册                             | 证明 FADBAD 是外部原始头文件而非复制实现       |

## 11. 数值验证通过时的源码 SHA-256 基准

捕获时间：2026-09-03T18:52:40Z。下列文件在生成本文前计算；本文自身不属于数值源码哈希。核心清单与 CMake 的 `BUBBLE_PROJECT_SOURCE_TREE_SHA256` 定义完全相同。

### 11.1 核心数值源码（权威基准）

| 相对路径                     | SHA-256                                                              |
| ---------------------------- | -------------------------------------------------------------------- |
| `CMakeLists.txt`           | `2f3f48d95237573fb6067762427266b46a6d12eb59691272781f4ff15850700b` |
| `include/bubble/model.hpp` | `60d17a13bd77fd15393a9e5ba95257a31b8d5748b47f51833c9c2be02273f877` |
| `src/config.cpp`           | `3e8a26adfd5fe842c9cf61e85b6b7ecff0be06159dabeb421319a94fcf6aea03` |
| `src/main.cpp`             | `99046a2898b053bbef67559f8f28a596310b0ebc61bc76fe68d3790894667757` |
| `src/network.cpp`          | `65edf8935ee6abe73da621ec7e7733962b2902c5b4a8ab317cd5dda5810dd734` |
| `src/output.cpp`           | `1ece9f7db710ea54dca5cd1fa1bf8e9d3d60642337bcd90ef4b298bba4999a8d` |
| `src/solver.cpp`           | `080e332fb3acf0cd45261491661496f47eeecafb540ad6dacd9d501cd148c38f` |

按 `relative_path + TAB + file_sha256 + LF`、上述排序拼接后的清单 SHA-256：

```text
1158c9bc2f4bfd3fd7b58312d5c5ec272adeaaf9f9c4d359bc7c123e052b37d8
```

该值与数值验证通过后 `output/full_network_smoke/run_manifest.json` 中记录的 `project_source_tree_sha256` 一致。

### 11.2 验证与辅助源码（补充基准）

| 相对路径                                    | SHA-256                                                              |
| ------------------------------------------- | -------------------------------------------------------------------- |
| `tests/test_ad_smoke.cpp`                 | `3bd81337554d0231b811e0eb9015d1931765d94650c23f5bcf117cafd2e55ae1` |
| `tests/test_bubble_solver.cpp`            | `fa34b474109967f0b0a49ca7d83819f59cb0e72357188fe3cebb6eb073825ae6` |
| `tests/audit_cli_output.py`               | `e01851986c4e7bfecc5c67c868f14a948d06a8da07e7facf17867f8f7d0da97c` |
| `tests/audit_correction_template.py`      | `29d4e1e23842c86812f42d944607e882cc4ccded2492a6fadb04fb0a3c44b51f` |
| `tests/audit_failure_manifest.py`         | `e77210435fcbcdb3f1a9cf701529c78fbb991a26502d59055b7594f75381b05e` |
| `tests/audit_full_network_output.py`      | `06e45de0233910708e8396cdf6dd39f55f04076779c7c580b0e02e300fcc3696` |
| `tests/expect_failure.py`                 | `8e63298d0f09ca89ef509c6f8982695e9d40da344209402472b40218d2c921df` |
| `tests/test_cli_restart.py`               | `e5acea1a88df6d355b8f026b6d80ef8167ed0ec0723afa57c63bc530185f4532` |
| `tests/test_failure_progress_manifest.py` | `be0d3b8976a0ccd669063bfc750aaf066a58af4bbc8b3cd0a90f3c80496b6f1b` |
| `tests/test_geometry_overlay.py`          | `174ac50b447bf5cda61b0fdbe83bf146cc42a48770305036cfceee14f694b818` |
| `tests/test_invalid_configs.py`           | `5411e5ac4f6b367ff8129d6a301f857776f9bc91219b00ca7ac767f7884487c2` |
| `scripts/audit_strict_geometry.py`        | `2506be4e3ce0c64343f2036b4ee0a877784726e4084764cf4e23c88eef2f91e5` |
| `scripts/geometry_correction_overlay.py`  | `82dc8dd4194cf405588dd2128ec9305cd0c3890f3da4ef36ab854bf3f989cd3b` |

### 11.3 只读参考和 AD 依赖身份

| 文件                                    | SHA-256                                                              |
| --------------------------------------- | -------------------------------------------------------------------- |
| `../trans-react.cpp`                  | `d70e9c44eb54e51418791a24d587ef3e12b75314bafcb6e1e3fcd0e92269444a` |
| `/opt/FADBAD++-2.1/FADBAD++/fadbad.h` | `a4b29e1d27f400b71b93e13f168b1b58ef2b08ae9b0305516340ae1721fc312c` |
| `/opt/FADBAD++-2.1/FADBAD++/badiff.h` | `e66d01fecfd5c47db7d2507021b89016b59e6d0c0c1da1a10e58f2b518ca0040` |
| `/opt/FADBAD++-2.1/FADBAD++/fadiff.h` | `ccf4d7b0f98eafdc49a07960d11adbc41e17e7baa5d24d7473cd44f95c97deb0` |

## 12. 建议的下一步确认点

如果进入实现阶段，建议用户只在下面三种范围中选一个，避免一次改动过宽：

1. **注释/命名范围：** 方程编号、单位、稳定恒等式、FADBAD 图生命周期说明；
2. **函数拆分范围：** 只拆 `assemble()`，不碰 Newton 和输出；
3. **精度双路径范围：** 在新增一致性测试后，共享 AD/long-double 的少数纯代数 helper。

在用户确认前，不应实施上述任何重构。
