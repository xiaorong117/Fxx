# 气泡溶解求解器甲方交付与操作手册

> 适用项目：`bubble_dissolution_solver_amgx`
> 文档用途：安装验收、建立算例、运行、监控、停止、恢复、导出数据，以及使用 Codex 协助操作。
> 科学模型的公式、离散和假设另见 [物理与数值实现](PHYSICS_IMPLEMENTATION_ZH.md)。

## 1. 软件身份和适用范围

本程序是 C++17 孔隙网络气泡溶解求解器。它使用同一套自动微分残量和
monolithic Newton 系统求解 `Rl/Rd/Rg/Rc`，支持：

- 液相压力驱动流动；
- 溶解气体对流、扩散和气液界面传质；
- 静止被困气泡的自由气体摩尔守恒；
- pressure inlet 或“入口总流量方程 + 出口压力”边界；
- Eigen、AMGX 和 oneMKL PARDISO 线性后端；
- 自适应 backward-Euler 时间步；
- active-set 气泡消失、保守型微小气泡批量退休；
- checkpoint/restart、守恒审计、20 段统计和压缩 VTP/PVD 输出。

当前模型不包含气泡迁移、气相网络流动、聚并、破裂、再侵入、多组分空气
分馏和热效应。因此应表述为“静止被困气泡溶解孔隙网络模型”，不能泛称为
完整动态两相流模型。

`mass_transfer_multiplier != 1` 的算例只允许用于加速可视化/数值预览，不能
用其物理时间预测真实实验。当前批量科研比较使用 `mass_transfer_multiplier=1`。

## 2. 交付目录和不可修改对象

推荐保持以下目录层次：

```text
bubble_dissolution_solver_amgx/
├── AGENTS.md                 Codex/代码代理操作约束
├── README.md                 项目入口
├── CMakeLists.txt            构建入口
├── include/bubble/           公共结构和接口
├── src/                      控制方程、Newton、后端和输出
├── configs/                  单算例及批量配置
├── scripts/                  输入转换、审计、调度和监控
├── tests/                    单元与回归测试
├── docs/                     用户和物理文档
└── output/
    ├── verification/         验收证据
    └── production/           正式运行结果
```

以下对象是只读来源或历史证据，禁止覆盖：

```text
/workspace/zz/bubble_dissolution_solver
/workspace/zz/bubble_dissolution_solver/output/full_network_air_water_relative_pressure_run01
/workspace/zz/zz_projects/35级配
output/verification/adaptive_dt/fixed_final4
output/verification/adaptive_dt/adaptive_final4
output/verification/backend_comparison
output/verification/batch_microbubble_retirement
```

原始 PNM 只能读取。转换、几何修订和单位归一化必须写入新的派生目录，并保存
源文件及派生文件 SHA-256。

## 3. 当前已验收生产版本

当前 35 级配 × 两种初始状态的正式批量版本使用：

```text
求解器：/workspace/zz/build-bubble-batch-retirement-release/bubble_solver
SHA-256：6ea4e9c1a3e00f3919e2415edeaa49e527bbd6f7c9042bbb98d7ee44695ebd04
后端：PARDISO
配置清单：configs/batch35_flow_pe555_pardiso_v3/batch_cases.json
派生输入清单：output/batch35_inputs/batch_input_inventory.json
生产输出：output/production/batch35_flow_pe555_pardiso/
```

该二进制启用了 PARDISO、未启用 AMGX。AMGX 的正确性测试已通过，但当前全网络
配置比 CPU 后端慢，不用于这批生产计算。请求一个未编译进二进制的后端会明确
失败，程序不会静默回退。

生产二进制与源码、配置、PNM 的哈希必须共同保存。不能在修改源码后把旧二进制
产生的数据描述为新版本结果。

## 4. 环境和构建

主要依赖：

- CMake 3.18 或更高；
- 支持 C++17 的编译器；
- Eigen3；
- FADBAD++ 2.1；
- nlohmann/json；
- zlib；
- 可选：oneMKL runtime（PARDISO）；
- 可选：CUDA + NVIDIA AMGX。

PARDISO Release 构建示例：

```bash
cd /workspace/zz
cmake -S bubble_dissolution_solver_amgx \
  -B build-bubble-client-release \
  -DCMAKE_BUILD_TYPE=Release \
  -DBUBBLE_ENABLE_AMGX=OFF \
  -DBUBBLE_ENABLE_PARDISO=ON \
  -DPARDISO_MKL_ROOT=/opt/conda/envs/Mesh
cmake --build build-bubble-client-release -j 8
```

构建完成后先保存：

```bash
sha256sum /workspace/zz/build-bubble-client-release/bubble_solver
```

再运行与本次修改有关的测试：

```bash
ctest --test-dir /workspace/zz/build-bubble-client-release \
  --output-on-failure
```

不能通过增大容差、把 `NOT_CONVERGED` 当成功或静默切换后端来“跑通”测试。

## 5. 配置文件的关键部分

### 5.1 算例和输出身份

每次新运行至少必须使用唯一的：

```json
{
  "case": {"name": "unique_case_name"},
  "io": {"output_directory": "/absolute/new/output/path"}
}
```

不得把新结果写入已有输出目录。

### 5.2 科研物理设置

当前批量对比的主要设置是：

```text
liquid                     pure water
gas                        effective dry air
temperature                298.15 K
background pressure        101325 Pa absolute
inlet species              degassed water, C_in=0
outlet liquid pressure     0 Pa relative
mass_transfer_multiplier   1
```

入口采用一个全局流量方程约束所有入口边的液体流量之和，入口共享压力由 Newton
反算，局部流量由网络导流能力自然分配。它不是给每个入口喉平均分配流量。

### 5.3 当前仍未实验确认的闭合参数

以下值虽在配置中有数值，但仍应标记为 `assumed` 或 `unconfirmed`：

- Henry 系数及有效空气组分近似；
- 液膜传质系数 `kL`；
- 等体积球形气液界面面积模型；
- Qin 毛细压力模型中的有效半径和参数 `b`；
- 接触角；
- 相对导流指数；
- 扩散面积指数与曲折度；
- 静止、不可再活化气泡假设；
- 入口流量范围与目标 Pe/Ca 的实验代表性。

这些参数未经实验标定前，结果适合数值对比和趋势分析，不应表述为已经验证的
绝对溶解时间预测。

## 6. 单算例启动

先确认没有同名输出：

```bash
test ! -e /workspace/zz/bubble_dissolution_solver_amgx/output/production/NEW_CASE
```

使用 tmux 和 `flock` 防止重复启动：

```bash
PROJECT=/workspace/zz/bubble_dissolution_solver_amgx
CONFIG=$PROJECT/configs/your_case.json
OUTPUT=$PROJECT/output/production/NEW_CASE
BINARY=/workspace/zz/build-bubble-client-release/bubble_solver

mkdir -p "$OUTPUT"
tmux new-session -d -s NEW_CASE \
  "cd /workspace/zz && exec flock -n '$OUTPUT/run.lock' \
   '$BINARY' --config '$CONFIG' >>'$OUTPUT/driver.log' 2>&1"
```

启动后必须确认：进程存在、配置路径正确、第一接受步写出、Newton 双门槛通过、
`epsilon_N/epsilon_V` 通过且 checkpoint 可读。

## 7. 70-case 批量计算

当前真实库存是 PNM-0.6 共 35 个、PNM-0.8 共 35 个，总计 70 个。输入准备、配置
生成和两个 `t=0` pilot 均有独立审计。正式调度命令形式为：

```bash
cd /workspace/zz/bubble_dissolution_solver_amgx
tmux new-session -d -s batch35_production_24 \
  "exec python3 scripts/run_batch35.py \
   --index configs/batch35_flow_pe555_pardiso_v3/batch_cases.json \
   --binary /workspace/zz/build-bubble-batch-retirement-release/bubble_solver \
   --jobs 24 --threads-per-job 8 \
   >>output/production/batch35_flow_pe555_pardiso/scheduler.log 2>&1"
```

一个清单只能有一个总调度器。不要再启动两个读取同一 `batch_cases.json` 的调度器。
`--jobs 24` 表示 24 个并发槽位，可理解为 3 组 × 8，但仍由一个进程统一分配。

批量进度：

```bash
jq '{status,jobs,running:(.running|length),completed:(.completed|length),failed}' \
  output/production/batch35_flow_pe555_pardiso/batch_progress.json
```

若出现一个失败，调度器会停止派发新任务并保留现场。不得删除失败目录后继续假装
全部完成。

## 8. 浏览器监控和远程访问

启动面板：

```bash
cd /workspace/zz/bubble_dissolution_solver_amgx
tmux new-session -d -s bubble_dashboard \
  "exec python3 scripts/case_dashboard.py \
   --host 127.0.0.1 --port 8765"
```

面板只监听 loopback。当前开发方 Mac 的 SSH 配置使用 `QIN-container` 别名，因此
转发命令是：

```bash
ssh -N -L 8765:127.0.0.1:8765 \
  -o ExitOnForwardFailure=yes QIN-container
```

浏览器访问 `http://127.0.0.1:8765/`。甲方必须建立自己的 SSH/Tailscale 账号、
密钥和 Host 别名，不能复制开发方私钥。`localhost` 地址只能由建立该转发的电脑
访问，不能作为公共网址发送。

如果启用控制模式，控制令牌必须为权限 0600 的独立文件。不要把令牌放入配置、
截图、聊天、Git 或交付压缩包。仅查看监控不需要控制令牌。

## 9. 正常状态判断

正常运行至少满足：

- `step` 和 `physical_time_s` 持续增加；
- `requested_dt_s`、`accepted_dt_s` 有限且位于上下限内；
- `residual_scaled_inf <= residual_tolerance`；
- `relative_update_norm <= update_tolerance`；
- `residual_criterion_passed=1` 且 `update_criterion_passed=1`；
- 没有 NaN、Inf、负气体摩尔数和长期连续 retry；
- `epsilon_N_adjusted`、`epsilon_V_adjusted` 小于配置硬门槛；
- `total_free_gas_moles` 总体下降；
- checkpoint 周期更新，且不落后 diagnostics 太多；
- `run_manifest.failed.json` 不存在。

微小气泡批量退休不会删除质量。剩余正自由气体摩尔数转入所属 cluster 的溶解
库存，退休气体体积进入 active-set 液体体积账本，然后仍需通过原守恒门槛。
`estimated_crossing_time_s` 是批量事件在接受步内的估计时间，不是重新求解得到的
严格事件时刻。

## 10. 安全停止和 checkpoint 恢复

单算例安全停止：

```bash
touch /absolute/output/CASE/STOP
```

程序在当前接受步结束后写 checkpoint、manifest 和终止原因
`safe_stop_requested`。确认进程退出后，把 STOP 移入归档目录，而不是删除证据：

```bash
mkdir -p /absolute/output/CASE/control
mv /absolute/output/CASE/STOP \
   /absolute/output/CASE/control/STOP.previous
```

恢复：

```bash
tmux new-session -d -s CASE_restart \
  "cd /workspace/zz && exec flock -n '/absolute/output/CASE/run.lock' \
   '/absolute/bubble_solver' --config '/absolute/case.json' \
   --restart '/absolute/output/CASE/checkpoint.json' \
   >>'/absolute/output/CASE/driver.log' 2>&1"
```

restart 会校验网络哈希、配置延续签名、节点数、压力表示、自适应状态和不可逆
active-set。不要为了绕过校验而手工修改 checkpoint。

批量运行时不要从网页单独停止某个子 case，因为总调度器会立即补位。整批维护
必须先阻止调度器继续派发，再让活动 case 安全停止，逐一验证 checkpoint，最后由
单一新调度器恢复。参考现有审计：

```text
output/production/batch35_flow_pe555_pardiso/
  concurrency_transition_8_to_24_20260910/handoff_audit.json
```

## 11. 输出文件及论文用途

| 文件 | 用途 |
|---|---|
| `diagnostics.tsv` | 每个接受步的 dt、Newton、线性残量、四残量、库存、边界流量和守恒 |
| `global_history.tsv` | 全局气体体积、自由/溶解摩尔数、总体和主贯通域 Sg |
| `segment_history.tsv` | 每步 20 个流向区段的 Sg、气泡数、气体体积和传质 |
| `segment_geometry.tsv` | 分段几何、孔隙体积和静态指标 |
| `gas_disappearance_events.tsv` | 每个气泡退休的节点、区段、模式、估计时间和转移量 |
| `balance_history.tsv` | 摩尔与液体体积守恒审计 |
| `checkpoint.json` | 完整 restart 状态和累计账本 |
| `state_*.vtp`, `states.pvd` | ParaView 空间场与时间序列 |
| `run.log`, `driver.log` | 求解器和外层启动日志 |
| `run_manifest.json` | 正常/安全结束后的配置、构建、输入、输出哈希和终止原因 |
| `run_manifest.failed.json` | 异常退出证据 |
| `termination_events.json` | 最终停止原因、目标值和最后事件 |

TSV/JSON 保留高精度用于统计；网页只缩短显示位数。论文数据必须直接读取机器文件，
不能从网页截图抄数。

正常完成应看到明确的 `termination_reason`，例如：

```text
main_real_gas_saturation_reached
residual_free_gas_fraction_reached
```

达到最大物理时间、磁盘保护或安全停止不等于达到科学终点，必须单独标记为删失或
中断数据。

## 12. 常见故障

### 找不到 PNM

核对配置中的绝对路径，或从配置预期的工作目录启动。不要复制一个同名错误文件来
绕过路径问题。

### dt 降到下限

从 `driver.log` 和 diagnostics 区分线性失败、Newton/Armijo、状态越界、守恒门槛
和气泡事件。不要首先降低 `dt_min` 或放宽容差。

### 面板显示“失败”

先区分“本批失败”和 production 根目录里的历史失败算例。权威批量状态来自
`batch_progress.json`。面板显示层不能代替 failure manifest。

### 面板无法访问

依次检查：

```bash
tmux has-session -t bubble_dashboard
ss -ltnp '( sport = :8765 )'
curl -fsS http://127.0.0.1:8765/api/cases
```

再检查客户端 SSH 隧道。不要把面板改为监听 `0.0.0.0` 后裸露到公网。

## 13. 用 Codex 操作本项目

让 Codex 的工作目录指向项目根目录，它会优先读取 [AGENTS.md](../AGENTS.md)。
推荐首次提示：

```text
先完整阅读 /workspace/zz/bubble_dissolution_solver_amgx/AGENTS.md、
docs/CLIENT_DELIVERY_GUIDE_ZH.md 和 docs/PHYSICS_IMPLEMENTATION_ZH.md。
只读审计当前 tmux、求解进程、batch_progress.json 和最后 diagnostics；
不要停止、重启、改配置或源码。报告当前版本、算例状态、守恒和异常。
```

要求 Codex 修改时，应明确：目标、允许修改的目录、禁止覆盖的输出、是否允许停止
进程、验收测试和最终报告字段。不要只说“优化一下”或“让它跑快一点”。

## 14. 甲方验收清单

- [ ] 项目、二进制、配置、PNM 和关键验证输出 SHA-256 已移交；
- [ ] Release 构建和相关回归测试通过；
- [ ] 代表性 `t=0` pilot 通过；
- [ ] backend 与 manifest 一致，无静默回退；
- [ ] Newton 双门槛和守恒门槛未放宽；
- [ ] checkpoint restart 连续性通过；
- [ ] 面板只能读取 production 根目录；
- [ ] 甲方使用独立 SSH 身份，不共享开发方私钥和控制令牌；
- [ ] 预览算例与科研基准明确分离；
- [ ] 未确认的物理闭合参数已列入报告；
- [ ] 正常终止、中断、失败和删失结果能明确区分。
