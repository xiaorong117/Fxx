# 气泡溶解求解器操作手册

本文面向需要自行编译、配置、启动、监视、停止和恢复计算的使用者。项目目录为：

```text
/workspace/zz/bubble_dissolution_solver_amgx
```

所有运行命令建议从 `/workspace/zz` 执行，否则配置中的相对 PNM 路径可能找不到。

## 1. 项目目录

```text
include/bubble/       公共数据结构和接口
src/                  网络、残量、Newton、输出和线性后端
configs/              JSON算例配置
tests/                单元测试和回归测试
scripts/              审计、比较和监控脚本
output/verification/  验证证据
output/production/    生产和预览算例
docs/                 操作和物理文档
```

不要覆盖原始工程和受保护基准：

```text
/workspace/zz/bubble_dissolution_solver
/workspace/zz/bubble_dissolution_solver/output/full_network_air_water_relative_pressure_run01
```

## 2. 构建

### 2.1 推荐的Eigen生产构建

```bash
cd /workspace/zz
cmake -S bubble_dissolution_solver_amgx \
  -B build-bubble-production-release \
  -DCMAKE_BUILD_TYPE=Release \
  -DBUBBLE_ENABLE_AMGX=OFF
cmake --build build-bubble-production-release -j2
```

程序位于：

```text
/workspace/zz/build-bubble-production-release/bubble_solver
```

当前全网络实测中 Eigen 比现有 AMGX 配置快；新加入的PARDISO后端在五步checkpoint
验收中又比Eigen获得约2.27倍进程墙钟加速。正式切换已有生产算例前仍应使用独立目录完成
对应配置的restart验收。

### 2.2 AMGX构建

```bash
cmake -S bubble_dissolution_solver_amgx \
  -B build-bubble-amgx-release \
  -DCMAKE_BUILD_TYPE=Release \
  -DBUBBLE_ENABLE_AMGX=ON \
  -DAMGX_ROOT=/opt/amgx_install
cmake --build build-bubble-amgx-release -j2
```

请求AMGX而构建中不可用时会明确失败，不会静默回退Eigen。

### 2.3 PARDISO构建

```bash
cmake -S bubble_dissolution_solver_amgx \
  -B build-bubble-pardiso-release \
  -DCMAKE_BUILD_TYPE=Release \
  -DBUBBLE_ENABLE_AMGX=OFF \
  -DBUBBLE_ENABLE_PARDISO=ON \
  -DPARDISO_MKL_ROOT=/opt/conda/envs/Mesh
cmake --build build-bubble-pardiso-release -j2
```

配置中的线性段使用：

```json
{
  "backend": "pardiso",
  "relative_tolerance": 1e-10,
  "max_iterations": 2000,
  "ilut_drop_tolerance": 1e-8,
  "ilut_fill_factor": 20,
  "pardiso_threads": 16
}
```

ILUT字段为向后兼容保留，PARDISO不使用它们。求解器对固定2×2块图只进行一次符号分析，
每个Newton迭代仍重新组装Jacobian并执行数值LU分解。若自由度或结构改变，会安全地重新
进行符号分析。请求PARDISO但构建未启用时会明确失败。

微小气泡批量退休默认关闭。经过对应算例验收后可在新算例中显式启用：

```json
"active_set": {
  "batch_microbubble_retirement": {
    "enabled": true,
    "maximum_individual_free_gas_fraction": 1e-10,
    "maximum_batch_free_gas_fraction_per_step": 1e-8
  }
}
```

不得删除原有 `Sg_off/ng_off`、负摩尔检查或守恒门槛。逐泡事件记录位于 `gas_disappearance_events.tsv`；批量总量超限时程序自动回到精确事件定位。

### 2.4 ASan/UBSan构建

```bash
cmake -S bubble_dissolution_solver_amgx \
  -B build-bubble-asan \
  -DCMAKE_BUILD_TYPE=Debug \
  -DBUBBLE_ENABLE_AMGX=OFF \
  -DBUBBLE_ENABLE_SANITIZERS=ON
cmake --build build-bubble-asan -j2
ctest --test-dir build-bubble-asan --output-on-failure
```

## 3. 运行测试

```bash
ctest --test-dir build-bubble-production-release --output-on-failure
```

常用CPU回归：

```bash
ctest --test-dir build-bubble-production-release --output-on-failure \
  -R 'bubble_small_cases|adaptive_dt_small_regression|adaptive_dt_end_time_roundoff|cli_restart_adaptive_equivalence'
```

主要测试：

- `bubble_small_cases`：AD Jacobian、四残量、虚拟边界、cluster、massless junction、active-set和守恒；
- `adaptive_dt_small_regression`：自适应时间步；
- `adaptive_dt_end_time_roundoff`：终止时间浮点尾数；
- `cli_restart_adaptive_equivalence`：连续运行和restart一致；
- `accelerated_preview_small_regression`：传质倍率、气泡消失和1%终止；
- `case_dashboard_provenance`：只读监控与路径限制。

## 4. 常用配置

### 科研基准方向

```text
configs/bubble_solver.air_water_adaptive_long_eigen.json
```

采用Eigen、自适应时间步和 `mass_transfer_multiplier=1`。物理闭合参数仍需实验确认。

### 压力边界并注入脱气水

```text
configs/bubble_solver.dirichlet_degassed_kLa1.json
```

```text
入口液相相对压力 = 100 Pa
出口液相相对压力 = 0 Pa
C_in = 0
C_backflow = 0
mass_transfer_multiplier = 1
```

### 总流量入口并注入脱气水

```text
configs/bubble_solver.transport_preview_Q2_C0.json
```

```text
sum(Q_inlet_edges) = 2.0e-9 m3/s
入口共享压力由Newton反算
出口液相相对压力 = 0 Pa
C_in = 0
mass_transfer_multiplier = 1
```

入口各喉流量由局部导流能力和压力自然分配，不是人为平均分配。

### 加速预览

```text
configs/bubble_solver.accelerated_preview.json
configs/bubble_solver.accelerated_fast_kLa1000.json
```

`mass_transfer_multiplier` 不等于1的结果只能用作数值或可视化预览，不能解释为真实物理时间。

## 5. 新建算例

复制最接近目标的配置，至少修改：

```json
{
  "case": {"name": "new_unique_case_name"},
  "io": {
    "output_directory":
      "bubble_dissolution_solver_amgx/output/production/new_unique_case_name"
  }
}
```

启动前确认输出目录不存在：

```bash
test ! -e /workspace/zz/bubble_dissolution_solver_amgx/output/production/new_unique_case_name
```

改变物性、边界类型、入口浓度或传质倍率后必须从 `t=0` 开始，不能续接不兼容checkpoint。

## 6. 直接运行

```bash
cd /workspace/zz
./build-bubble-production-release/bubble_solver \
  --config /workspace/zz/bubble_dissolution_solver_amgx/configs/your_case.json
```

## 7. 使用tmux长期运行

```bash
CASE_NAME=new_unique_case_name
CONFIG=/workspace/zz/bubble_dissolution_solver_amgx/configs/your_case.json
OUTPUT=/workspace/zz/bubble_dissolution_solver_amgx/output/production/new_unique_case_name
BINARY=/workspace/zz/build-bubble-production-release/bubble_solver

mkdir -p "$OUTPUT"
touch "$OUTPUT/run.lock"
tmux new-session -d -s "$CASE_NAME" \
  "cd /workspace/zz && exec flock -n $OUTPUT/run.lock \
   $BINARY --config $CONFIG >$OUTPUT/driver.log 2>&1"
```

检查：

```bash
tmux ls
pgrep -af bubble_solver
tmux attach -t new_unique_case_name
```

从tmux脱离但不停止程序：按 `Ctrl-b`，再按 `d`。

## 8. 查看进度

```bash
tail -f output/production/new_unique_case_name/diagnostics.tsv
tail -f output/production/new_unique_case_name/run.log
tail -2 output/production/new_unique_case_name/global_history.tsv
jq '{step_number,time_s,next_dt_s,gas_disappearance_events}' \
  output/production/new_unique_case_name/checkpoint.json
```

重点检查：

```text
step, time_s, dt_s, retries
newton_iterations, linear_iterations
residual_scaled_inf, relative_update_norm
epsilon_N_adjusted, epsilon_V_adjusted
overall_Sg, main_connected_component_Sg
total_free_gas_moles, total_dissolved_moles
gas_disappearance_events
```

## 9. 浏览器监控

服务器启动：

```bash
cd /workspace/zz/bubble_dissolution_solver_amgx
tmux new-session -d -s bubble_dashboard \
  'python3 scripts/case_dashboard.py --host 127.0.0.1 --port 8765'
```

Mac建立SSH转发：

```bash
ssh -N -L 8765:127.0.0.1:8765 用户名@服务器地址
```

访问：

```text
http://127.0.0.1:8765/
```

监控程序只读，且只扫描：

```text
/workspace/zz/bubble_dissolution_solver_amgx/output/production/
```

面板会递归识别该目录内的批量算例，但拒绝访问此根目录之外的路径。左侧
“全部生产算例”显示运行中、排队中、已完成和失败数量，并可按算例名、
`PNM-0.6`/`PNM-0.8` 与运行状态筛选。列表接口只读取末行诊断和调度状态；
完整曲线、20 段统计、守恒记录及 PNM 溯源仅在选中某个算例后加载，因此
可用于数十个网络的长期批量监控。批量调度目录在控制面板中保持只读，
避免单独停止或重启某个子任务而绕过批量 `run.lock`。

## 10. 安全停止

### 10.1 在监控面板启用受限控制

默认dashboard保持只读。需要控制按钮时，先建立只允许当前用户读取的随机令牌：

```bash
cd /workspace/zz/bubble_dissolution_solver_amgx
umask 077
openssl rand -hex 32 > .dashboard-control-token
chmod 600 .dashboard-control-token
```

以控制模式启动：

```bash
python3 scripts/case_dashboard.py \
  --host 127.0.0.1 \
  --port 8765 \
  --enable-control \
  --control-token-file /workspace/zz/bubble_dissolution_solver_amgx/.dashboard-control-token \
  --solver-binary /workspace/zz/build-bubble-flow-release/bubble_solver
```

通过SSH在服务器终端读取令牌，并粘贴到网页的“控制令牌”输入框：

```bash
cat /workspace/zz/bubble_dissolution_solver_amgx/.dashboard-control-token
```

不要把令牌写入聊天、截图、配置、manifest或版本库。

控制面只开放：

- 创建STOP进行安全停止；
- 从最新checkpoint恢复；
- 暂存 `dt_max_s`；
- 暂存checkpoint间隔；
- 暂存VTK间隔。

参数修订不会覆盖原配置，而是写入：

```text
output/production/CASE/control/config_revision_XXXX.json
output/production/CASE/control/staged_config.json
```

修订只在下一次restart生效。物性、边界条件、PNM、传质倍率、active-set阈值和Newton容差不允许热修改，必须复制为新算例并从 `t=0` 启动。

全部成功和失败操作记录在：

```text
output/production/CASE/control/control_history.jsonl
```

控制服务仍只监听loopback，并验证token、请求来源、算例名、production路径、配置输出目录、checkpoint、固定二进制和run.lock。网页不能提交任意shell命令或任意文件路径。

### 10.2 命令行安全停止

```bash
touch output/production/new_unique_case_name/STOP
```

求解器在接受步后保存checkpoint和最终输出，终止原因写为 `safe_stop_requested`。不要优先使用 `kill -9`。

## 11. 从checkpoint恢复

确认原进程已经退出，移走STOP：

```bash
mv output/production/new_unique_case_name/STOP \
   output/production/new_unique_case_name/STOP.previous
```

恢复运行：

```bash
tmux new-session -d -s new_unique_case_name_restart \
  "cd /workspace/zz && exec flock -n \
   /workspace/zz/bubble_dissolution_solver_amgx/output/production/new_unique_case_name/run.lock \
   /workspace/zz/build-bubble-production-release/bubble_solver \
   --config /workspace/zz/bubble_dissolution_solver_amgx/configs/your_case.json \
   --restart /workspace/zz/bubble_dissolution_solver_amgx/output/production/new_unique_case_name/checkpoint.json \
   >>/workspace/zz/bubble_dissolution_solver_amgx/output/production/new_unique_case_name/driver.log 2>&1"
```

restart会校验配置延续签名、节点数、压力表示、PNM哈希、自适应状态和active-set不可逆状态。流量入口还会恢复共享入口压力。

## 12. 时间步配置

```json
"time_step_control": {
  "enabled": true,
  "dt_min_s": 1e-6,
  "dt_max_s": 0.1,
  "growth_factor": 1.25,
  "shrink_factor": 0.5,
  "easy_newton_iterations": 3,
  "normal_newton_iterations": 7,
  "hard_newton_iterations": 12,
  "gas_timescale_fraction": 0.1,
  "gas_timescale_percentile": 0.5,
  "gas_timescale_min_mole_fraction": 1e-5,
  "concentration_change_limit_enabled": false,
  "max_consecutive_growth": 3
}
```

不要为“跑通”而随意放宽Newton或守恒容差。

## 13. 科学终止

```json
"termination": {
  "primary": "residual_free_gas_fraction",
  "target_gas_saturation": 0.01,
  "residual_free_gas_fraction": 0.01,
  "maximum_physical_time_s": 1000000,
  "maximum_wall_seconds": 0,
  "maximum_output_bytes": 12884901888
}
```

`maximum_wall_seconds=0` 表示没有墙钟限制。最大物理时间先达到时，结果属于右删失，不等于稳定残余状态。

## 14. 输出文件

- `diagnostics.tsv`：逐接受步Newton、残量、流量、事件和守恒；
- `global_history.tsv`：全局论文作图数据；
- `segment_history.tsv`：逐步逐段空间统计；
- `segment_geometry.tsv`：分段静态几何；
- `balance_history.tsv`：摩尔和液体体积守恒；
- `checkpoint.json`：重启状态；
- `state_*.vtp` 与 `states.pvd`：ParaView时间序列；
- `run.log` 与 `driver.log`：运行日志；
- `run_manifest.json`：正常或安全结束后的可复现清单；
- `run_manifest.failed.json`：异常退出清单；
- `termination_events.json`：终止和阈值事件。

机器文件保留高精度，浏览器仅压缩显示位数。

## 15. ParaView

打开 `states.pvd`，推荐字段：

```text
Pl_relative_Pa
Pg_relative_Pa
Sw
Sg
C
active_gas
permanently_inactive
gas_moles
dissolved_moles
segment_id
```

绝对压力关系为：

```text
P_absolute = P_background + P_relative
```

## 16. 常见故障

### 找不到PNM

`cannot open geometry file` 通常表示工作目录错误。请从 `/workspace/zz` 启动。

### dt下降到下限

查看 `driver.log` 的完整原因，区分线性失败、Armijo失败、状态边界、消失事件和守恒失败，不要直接减小 `dt_min`。

### manifest缺失

运行中通常没有最终manifest；正常完成或安全停止后才生成。

### dashboard显示中断

```bash
pgrep -af bubble_solver
stat output/production/CASE/diagnostics.tsv
tail output/production/CASE/driver.log
```

## 17. 生产前检查

- JSON可被 `python3 -m json.tool` 解析；
- case name和输出目录唯一；
- 边界类型、入口浓度和物性符合目标；
- 科研算例 `mass_transfer_multiplier=1`；
- 使用Eigen后端；
- 从 `t=0` 开始或checkpoint完全兼容；
- PNM路径和哈希正确；
- 时间步、终止和磁盘限制明确；
- 回归测试通过；
- tmux、PID、第一接受步和守恒均已确认。
