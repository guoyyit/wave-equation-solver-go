# wave-equation-solver-go

二维标量波动方程的数值求解与恢复服务。输入完全在内存中构造（库 API 或本地 CLI 标志），不依赖外部网络或真实凭据；输出为指定时刻与空间点上的离散解及诊断信息（JSON）。

## 求解问题

在矩形区域 `[0,Lx] x [0,Ly]` 上求解二阶波动方程

```
u_tt = C^2 (u_xx + u_yy)
```

带初始位移 `u(x,y,0)`、初始速度 `u_t(x,y,0)`，每条边可独立选择**固定边界**（Dirichlet，`u=0`）或**吸收边界**（一阶 Mur）。

## 离散格式

- **空间**：规则网格 `Nx x Ny`，五点中心差分（二阶精度）。
- **时间**：显式 leapfrog（中心时间差分，二阶精度）：
  `u^{n+1} = 2u^n - u^{n-1} + (C·Dt)^2 · ∇²_h u^n`。
- **启动步**：由初始条件用 Taylor 展开构造 `u^{-1} = u^0 - Dt·v^0 + (Dt^2/2)·C^2·∇²u^0`，保证从第一步起整体二阶精度。
- **固定边界**：边界点恒置 0（含 t=0 与 u^{-1}，保证自洽）。
- **吸收边界**：一阶 Mur 条件，例如左边界
  `u^{n+1}_{0,j} = u^n_{1,j} + (r-1)/(r+1)·(u^{n+1}_{1,j} - u^n_{0,j})`，`r = C·Dt/Dx`。
  CFL 满足时系数模 `< 1`，边界更新本身不会引入不稳定。
- **角点**：任一边为固定则角点置 0；两边均吸收则取两个一维 Mur 更新的平均——混合边界组合自洽。

## 稳定性约束

显式格式受 **CFL 条件**约束：

```
r = C · Dt · sqrt(1/Dx^2 + 1/Dy^2) <= 1
```

`Config.Validate()` 在求解前检查；`r > 1` 的配置被**拒绝并说明原因**（给出当前 r 与最大允许 Dt）。此外 `Run` 每步检查场值有限性，出现非有限值立即报错，绝不产生无界发散的"成功"结果。

## 初始条件契约

初始位移与初始速度必须**显式提供**：要么给出长度 `Nx*Ny` 的场，要么设置对应 `Zero` 标志声明其为零。缺失、长度不符、或"既给场又置 Zero 标志"的矛盾输入都会报错，绝不静默采用默认值。两者同为零的平凡问题同样被拒绝。

## 检查点与恢复

- 设置 `CheckpointPath` 后，`Run` 每 `CheckpointEvery` 步及结束时提交检查点。
- 提交是**原子的**：先写 `*.tmp`，`fsync` 后 `rename` 覆盖旧检查点——崩溃只会留下旧或新的完整检查点，不会有部分状态。
- 检查点包含 leapfrog 所需的两个时间层与步数，并带配置**指纹**（FNV-1a，覆盖所有影响数值结果的字段；`Workers`、`CheckpointEvery` 等调度字段除外）。指纹不匹配即拒绝恢复，防止用错误输入续算。
- 恢复语义：新建 `Solver` 自动加载最近提交的检查点并从该步继续；`Run` 推进到 `Steps` 为止，已完成的求解器上调用 `Run` 是空操作——不重复推进。相同输入下，中断-恢复与不间断运行的结果**逐位一致**。

## 可复现与并行

- 内部更新按行分块并行（`Workers`，`<=0` 表示 `NumCPU`）。每个输出单元只依赖上一时间层，无归约顺序问题，因此**单线程与多线程结果逐位一致**（满足"给定误差内一致"，实际为零误差）。
- 固定参数（含 `Seed`）重复运行结果逐位一致。
- 性能随并行度可观测提升（M1 Pro，512x512x20 步：1 worker 16.4ms → 4 workers 7.5ms）：
  `go test ./internal/wave/ -run xxx -bench BenchmarkRun`

## 配置项（`wave.Config`）

| 字段 | 含义 | 约束 |
|---|---|---|
| `Nx, Ny` | 网格点数 | `>= 2` |
| `Lx, Ly` | 区域边长 (m) | `> 0` |
| `C` | 波速 (m/s) | `> 0` |
| `Dt` | 时间步长 (s) | `> 0` 且满足 CFL |
| `Steps` | 总步数 | `>= 1` |
| `Boundaries` | 四边类型：`BoundaryFixed` / `BoundaryAbsorbing` | 枚举值合法 |
| `Init` | 初始位移/速度 + `Zero` 标志 | 见"初始条件契约" |
| `Seed` | 预留随机种子 | 计入检查点指纹 |
| `Workers` | 并行度 | `<=0` 表示 NumCPU |
| `CheckpointPath` | 检查点文件路径 | 空则禁用检查点 |
| `CheckpointEvery` | 提交间隔（步） | 启用检查点时 `>= 1` |

## CLI 使用

```sh
# 基本运行：中心高斯脉冲，四边吸收，输出诊断与采样点
go run ./cmd/wavesim -nx 64 -ny 64 -steps 200 -dt 0.008 -sample "32,32;10,10"

# 混合边界 + 检查点；中断后重跑同一命令即从最近检查点恢复
go run ./cmd/wavesim -boundary-left fixed -boundary-top fixed \
  -checkpoint /tmp/wave.ckpt -checkpoint-every 20 -steps 1000
```

输出为 JSON：`Resumed`（是否从检查点恢复）、`Diagnostics`（能量、最大幅值、已完成步数）、`Samples`（指定网格点的解值）。所有标志见 `go run ./cmd/wavesim -h`。

## 测试

```sh
go test ./...
```

覆盖：配置校验与稳定性拒绝（`config_test.go`）、16 种边界组合自洽性（`boundary_test.go`）、吸收边界反射能量对比（`absorbing_test.go`，残余能量降至全固定参考的 ~0.04%）、崩溃恢复/原子提交/指纹校验（`checkpoint_test.go`）、可复现性与并行一致性及扩展性基准（`parallel_test.go`）。

## 目录结构

```
cmd/wavesim/        本地 CLI（内存构造输入，JSON 输出）
internal/wave/
  config.go         配置、边界类型、校验（CFL、初始条件契约）
  grid.go           规则网格存储
  solver.go         leapfrog 求解器、并行步进、Run/恢复语义
  boundary.go       固定/Mur 吸收边界与角点规则
  checkpoint.go     原子检查点提交与加载、配置指纹
  diagnostics.go    能量、最大幅值、有限性检查
```
