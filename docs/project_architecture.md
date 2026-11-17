# 项目分层与 Workspace 迁移方案

日期：2026-11-17。状态：基础库与 PLT 已完成 Workspace 源码迁移；验证与路径映射见 [迁移记录](workspace_migration.md)。

## 1. 建议采用的结构

保留一个仓库，使用 Cargo Workspace 管理各层。按职责拆成独立 package/crate，每个 crate 内继续使用模块和多个文件。
当前已分离基础库与 PLT；C、D 在开始实现相应功能时建立，不提前创建空壳。

| 层次 | 建议 package | 负责的内容 | 状态 |
| --- | --- | --- | --- |
| 原始 A | 上游公开提交与本地历史备份 | 原始实现和历史对照 | 来源可定位，见下文 |
| 维护后的 A | `qfall-tools` | 原有基础工具、MP12 与必要的随机源/错误处理/数值修复 | 已位于 `crates/qfall-tools` |
| B | `plt` | 统一 Setup、PTrapGen、PSampPre、Rec，门限配置、参数评估 | 已位于 `crates/plt` |
| C | `protocol-core`（暂名） | 基于 PLT 的上层密码功能和协议核心 | 未实现，具体职责待设计 |
| D | `vehicle-auth`（暂名） | 车联网场景接入、应用配置、运行入口与系统实验 | 未实现 |

依赖方向：`vehicle-auth -> protocol-core -> plt -> qfall-tools`。
这是功能依赖方向，不要求把 serde 等共享依赖也逐层转发；各 crate 可以显式声明自己使用的通用依赖。
基础库不依赖 PLT，PLT 不依赖上层协议。每层都不复制下层源代码。

当前目录与后续扩展位置（C、D 尚未建立）：

```text
tools/
├── Cargo.toml                 # 虚拟 Workspace；依赖版本与构建配置
├── Cargo.lock                 # 整个工作空间的锁文件
├── README.md                  # 项目入口和各层导航
├── LICENSE
├── crates/
│   ├── qfall-tools/           # 维护后的基础库
│   │   ├── Cargo.toml
│   │   ├── src/
│   │   └── benches/           # 原有 PSF 基准
│   ├── plt/                   # B：PLT 算法
│   │   ├── Cargo.toml
│   │   ├── src/
│   │   ├── tests/
│   │   ├── examples/          # PLT 驱动、审计、公开求解筛查
│   │   └── benches/
│   ├── protocol-core/         # C：实现时创建
│   └── vehicle-auth/          # D：实现时创建
├── docs/                      # 分层说明、验收清单与历史实验记录
├── scripts/                   # 跨 crate 的运行与独立验证脚本
├── experiments/               # 需要系统实验时再建立配置和场景
└── target/                    # 构建与可丢弃的运行产物
```

Cargo Workspace 共享锁文件和构建目录，并可按 package 或整个工作空间执行命令；虚拟 Workspace 需要显式指定 resolver。
本项目已使用 edition 2024，建议明确写 `resolver = "3"`。参见 [Cargo Workspace 文档](https://doc.rust-lang.org/cargo/reference/workspaces.html)。

当前仅登记已经实现的两个 package；下面是关键配置片段：

```toml
[workspace]
members = ["crates/qfall-tools", "crates/plt"]
resolver = "3"

[workspace.package]
edition = "2024"
rust-version = "1.88"

[workspace.dependencies]
qfall-tools = { path = "crates/qfall-tools" }
qfall-math = "=0.2.1"
# 其余依赖保留当前约束，按各成员的实际使用情况继承。

[profile.bench]
debug = true
```

完整配置见 [根 Cargo.toml](../Cargo.toml)。Workspace 中声明依赖不会自动将其加入每个成员；成员仍须显式继承。
criterion 等仅用于测试/基准的依赖应在核查使用位置后归入相应成员的 dev-dependencies。
本次迁移未升级外部依赖，也未改变算法参数。

## 2. A 与 B 的来源及路径归属

迁移前的 Git 基线：`851bacc722cfb52f664127c441f5f658729213fd`，提交标题 `Update testing pipeline`。
原始仓库为 `https://github.com/qfall/tools.git`；项目 origin 为 `https://github.com/wydnzdgz/P-R-S.git`。
项目已改用独立开发历史，并移除本地 upstream 引用；原提交作为来源对照保留在文档与本地备份中，见[历史整理记录](project_history.md)。

源码迁移前，许多新增文件曾是 untracked，当时仅比较提交历史会漏掉这些 B 层实现；迁移检查点包含了这些文件。
当前各层代码已纳入项目提交。今后的迁移仍应同时核对已跟踪差异与未跟踪文件，避免备份遗漏。

| 迁移前路径 | 归属与迁移目标 | 理由 |
| --- | --- | --- |
| `src/compression*`、`src/utils*` | 原基础库，进入 `crates/qfall-tools/src/` | 与 PLT 业务无关 |
| `src/primitive/psf*` | 基础库及其维护修复 | 包含原 GPV/MP12，不属于 PLT 算法本体 |
| `src/primitive/psf/mp_perturbation/checked.rs` | 维护后的基础库 | 新增的共享 fallible MP12 实现，供旧 PSF 和 PLT 调用 |
| `src/sample/randomness.rs` | 维护后的基础库 | 新 OS 熵及采样随机源支持，属于采样器职责 |
| `src/sample/g_trapdoor/gadget_classical.rs`、`trapdoor_distribution.rs` | 基础库的必要修改 | 三元采样和错误传播修复，需保留 |
| `src/primitive/plt.rs`、`src/primitive/plt/` | `crates/plt/src/` | PLT 的统一算法、类型、配置、适配和评估 |
| `src/primitive.rs` | 基础库保留，移除 PLT 模块登记 | 基础库不再包含 B |
| `src/sample.rs` | 基础库保留维护改动 | 继续登记底层随机源模块 |
| `tests/plt_correctness.rs`、`tests/plt_sampling.rs` | `crates/plt/tests/` | PLT 的正确性及采样集成测试 |
| `examples/plt*.rs`、`examples/support/` | `crates/plt/examples/` | 现有例程均服务于 PLT 实验/审计，包括直接调用基础采样器的参照实验 |
| `benches/plt.rs` | `crates/plt/benches/`，增加独立基准入口 | PLT 基准独立运行 |
| `benches/psf.rs` | `crates/qfall-tools/benches/` | 原基础采样器基准 |
| `benches/benchmarks.rs` | 按归属拆成两个成员的基准登记 | 不再由基础库登记 PLT 基准 |
| `scripts/` | 根目录保留 | 扫描、独立整数验证与证据复核；调整构建命令和源码收集路径 |
| `docs/` | 根目录保留 | 文档可跨层引用，历史证据保留原样 |
| 根 Cargo.toml、Cargo.lock、CI、README | 项目管理与复现配置 | 不算某个算法层的实现 |

核心区别是：**原始 A 是历史对照，实际被 B 调用的是带有必要修复的 A。**
不应为了让目录看起来“原样保留”，把已修复的随机源和数值处理退回旧实现；也不应在 B 再复制一套 MP12。
后续可用单独的基础库维护说明记录相对上述基线的每项修复。移动文件时保留已有版权头和 LICENSE。

## 3. B 内部继续按职责组织

每个 crate 可以包含很多文件。PLT 本轮整体提取，现有目录为：

```text
crates/plt/src/
├── lib.rs               # 类型、统一算法、公开接口和模块登记
├── params.rs            # PltConfig 与参数合法性
├── mp12_adapter.rs      # PLT 到基础采样器的接口与宽度约定
├── vandermonde.rs       # 门限目标分解
├── assessment.rs        # 参数与数值诊断
└── acceptance_tests.rs  # 访问内部状态的故障注入测试
```

后续 lib.rs 增长时，可再提取 types.rs、error.rs 或 algorithms.rs；本次没有为搬迁重新组织算法主体。
不为每个函数强制创建文件，也不因文件增加就增加一个 crate。
跨层边界对应 crate；层内按相关职责使用 module。
`lib.rs` 对外导出统一算法及必要类型；中间计算、内部状态和故障注入接口继续限制可见性。
当前实现已有足够多的公开 MP12 fallible 入口，拆分时应先复用它们，避免机械地把所有内部函数改为 pub。

提取后调用路径从：

```rust
use qfall_tools::primitive::plt::{PltConfig, setup};
```

变为：

```rust
use plt::{PltConfig, setup};
```

`mp12_adapter.rs` 中原先指向 `crate::primitive::psf` 和 `crate::sample` 的引用已改为基础库依赖路径。
不要让基础库重新导出 PLT 来保留旧路径，那会形成基础库与 PLT 的循环依赖。
若未来确实需要兼容旧路径，再讨论独立的兼容入口；当前研究项目优先明确依赖边界。

上层通过配置和算法接口使用 PLT，不在协议消息处理代码中读取 R、U、Gram 逆或重新实现采样公式。
研究审计需要的秘密导出接口可以在 PLT 审计入口使用，不应自然成为协议消息格式。

## 4. 如何理解“下层写好后不再改”

目标是稳定公开接口、控制依赖方向和记录版本。下层仍可能需要缺陷修复、参数调整和性能改进。
修复应落在负责该行为的层，并回归依赖它的上层；不在 C、D 各自打补丁。
功能验收、接口稳定与安全参数验收要分别记录。

当前 PLT 已有通用接口和较多工程测试，但 B3 尚未通过，所测运行实例存在公开原像反例。
因此可以先整理结构、开发上层接口和测试框架，尚不能将 B 标记为已完成安全验收的冻结版本。
具体证据仍以 [验收清单](acceptance_checklist.md) 和 [参数评估](plt_parameter_exploration.md) 为准。

版本管理建议：一份 Git 仓库记录来源与变更；Workspace 表达构建依赖；各 crate 的模块表达职责。
分支不替代模块。阶段 tag 应区分流程原型、接口版本和安全验收状态，不用一个“完成”标签概括所有性质。

## 5. 迁移过程与后续验收原则

1. **建立迁移前检查点。** 包含当前未提交/未跟踪的算法、测试、脚本和需保留的证据；明确排除构建产物及无关本地配置。记录 Git 对照点、锁文件和当前验证命令。
2. **只分离两个现有 crate。** 先整体移动 PLT，再整理基础库的位置；迁移过程中保持算法公式、宽度、范数界、随机源和拒绝行为。
3. **更新调用方与工具。** 示例、测试和基准使用新的 package/import 路径；扫描脚本显式指定 `-p plt`，源码哈希收集覆盖成员 src；CI 使用 Workspace 命令并保留锁文件复现流程。现有 CI 的 `cargo update` 不应混入固定依赖的验收任务。
4. **更新文档路径。** 修复指向移动源文件的相对链接；历史报告中的旧命令/旧路径作为历史记录保留，补充新的复测入口。
5. **验收边界。** 原 PSF 与 PLT 的单元/集成/示例/文档测试、格式检查和 Clippy 通过；基准可构建，独立 Python 验证通过。静态检查普通依赖方向，基础库不引用 PLT。
6. **补充迁移后证据。** 重跑小规模 PLT 输出并独立核验，再运行代表性采样/数值诊断。保留历史 JSON 与源码指纹，新增记录说明路径和代码版本变化；不把随机样本必须逐字一致作为验收条件。
7. **再按需建立 C、D。** 每层各有测试、示例和基准；端到端测试属于最上层。稳定接口后再考虑替换采样后端所需的 trait，不提前搭建复杂抽象。

现在可使用的常用命令：

```bash
cargo test --locked -p qfall-tools
cargo test --locked -p plt
cargo test --locked --workspace --examples
cargo test --locked --workspace --doc
cargo fmt --all -- --check
cargo clippy --locked --workspace --all-targets -- -D warnings
cargo bench --locked -p plt --no-run
cargo run --locked -p plt --example plt -- --threshold 3 --participants 5
```

本次实现和验证结果见 [迁移记录](workspace_migration.md)。历史原像反例与 B3 拒绝状态保留；代码分层不构成新增安全结论。
