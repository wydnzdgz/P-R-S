# 源码迁移记录

日期：2026-11-17。状态：基础库与 PLT 已拆分，迁移后的本地验证通过。

## 1. 现在从哪里找代码

| 内容 | 当前目录或入口 |
| --- | --- |
| 根构建配置 | [Cargo.toml](../Cargo.toml)，虚拟 Workspace |
| A：维护后的基础库 | [crates/qfall-tools](../crates/qfall-tools/README.md) |
| 基础采样器实现 | [MP12](../crates/qfall-tools/src/primitive/psf/mp_perturbation.rs)、[fallible 核心](../crates/qfall-tools/src/primitive/psf/mp_perturbation/checked.rs) |
| 基础库来源与修复 | [baseline_maintenance.md](baseline_maintenance.md) |
| B：PLT | [crates/plt](../crates/plt/README.md) |
| PLT 类型与统一算法 | [lib.rs](../crates/plt/src/lib.rs) |
| PLT 参数、采样适配、评估 | [params.rs](../crates/plt/src/params.rs)、[mp12_adapter.rs](../crates/plt/src/mp12_adapter.rs)、[assessment.rs](../crates/plt/src/assessment.rs) |
| PLT 集成测试 | [正确性](../crates/plt/tests/plt_correctness.rs)、[采样与参数](../crates/plt/tests/plt_sampling.rs) |
| PLT 例程与实验 | [examples](../crates/plt/examples/) |
| 分层基准 | [基础采样器](../crates/qfall-tools/benches/)、[PLT](../crates/plt/benches/) |
| 跨层实验脚本 | [scripts](../scripts/) |

依赖为 `plt -> qfall-tools`。根目录共享一个 Cargo.lock 与 target；C、D 尚未建立。
PLT 本轮整体提取，已有内部模块继续使用，没有为迁移重新拆写算法。

## 2. 变更范围与保留内容

- 移动 43 个既有文件，其中 31 个移动前后逐字节一致。其余涉及跨 crate 的导入、模块登记、示例命令说明或基准说明。
- 基础库的 MP12、随机源、错误传播和数值检查实现保留。PLT 的算法公式、参数推导、采样宽度、范数界和拒绝策略保留。
- 外部依赖的 113 个锁定 package 的名称、版本和来源与迁移前一致；Cargo.lock 更新了两个本地 package 的依赖关系。
- 库依赖与测试/例程依赖分开：基础库不再携带 PLT；PLT 的 criterion、JSON 驱动、公开求解随机尝试等依赖由其自身 dev-dependencies 管理。
- 6 个既有历史 JSON 保持原始字节和哈希。文档源码链接及当前复测命令适配新路径；历史基准报告的原命令保留为历史记录。
- 没有新建 Git 提交或改写历史，基线对照仍为 `851bacc722cfb52f664127c441f5f658729213fd`。

逐文件的旧路径、新路径、移动前后 SHA-256、A/B 归属，以及基础库文件与上游的对应关系见
[迁移清单](workspace_migration_manifest.json)。运行及验证摘要见 [结果记录](workspace_migration_results.json)。

## 3. 调用方法变化

PLT 从基础库中独立出来，使用新的 crate 路径：

```rust
use plt::{PltConfig, setup, ptrap_gen, psamp_pre, rec};
```

旧的 `qfall_tools::primitive::plt` 路径移除；基础库的 `qfall_tools::primitive::psf` 等路径保留。
不在基础库中反向导出 PLT，以保持单向依赖。

从工作空间根目录执行：

```bash
# 分别测试 A、B
cargo test --locked -p qfall-tools
cargo test --locked -p plt

# 全部测试及例程测试
cargo test --locked --workspace
cargo test --locked --workspace --examples

# 运行 PLT，门限仍通过参数配置
cargo run --locked -p plt --example plt -- \
  --threshold 3 --participants 5 --output target/plt_3_of_5.json
python3 scripts/verify_plt_report.py target/plt_3_of_5.json

# 分别运行基准
cargo bench --locked -p qfall-tools --bench benchmarks
cargo bench --locked -p plt --bench benchmarks -- PLT3of5

# 分层与代码检查
python3 scripts/verify_workspace_layout.py
cargo fmt --all -- --check
cargo clippy --locked --workspace --all-targets -- -D warnings
```

`plt_2_of_3` 兼容例程仍可通过 `-p plt --example plt_2_of_3` 使用；它继续调用统一驱动。

## 4. 迁移后的实际验证

| 验证 | 结果 |
| --- | --- |
| 基础库单元测试 | 90 项通过 |
| PLT 单元测试 | 7 项通过 |
| PLT 集成测试 | 9 项通过 |
| 例程测试 | 7 项通过 |
| 文档测试 | 基础库 50 项、PLT 1 项通过；合计 Rust 测试仍为 164 项 |
| 构建与静态检查 | 全目标 check、Clippy、格式检查及禁止文档警告的 doc 构建通过 |
| 两层基准 | `cargo bench --locked --workspace --no-run` 通过；本轮不重新生成正式性能结论 |
| 全集合整数报告 | 2-of-3、3-of-5、5-of-10，每集合 3 个目标，共 795 次恢复，经 Python 独立核验四条同余式、整数乘积/求和及范数 |
| 重复运行 | 每配置两套密钥、每套 20 次，三组配置共 120 次恢复通过，未记录失败或跳过 |
| 兼容入口 | `plt_2_of_3` 运行和 3 次恢复的独立核验通过 |
| 完整 B 类工程筛查复跑 | 27,000 次 PLT 部分采样、20,000 次底层离散参照采样及 54 个实际矩阵诊断通过；独立 Python 核验通过 |
| 参数扫描复跑 | n=1,q=257,2-of-3 与 n=2,q=2147483647,3-of-5 均完成；新报告记录两个 crate 的源码与 manifest 哈希 |
| 公开格筛查 | 两个复测参数点均再次找到公开原像；Python 独立核验通过 |
| 脚本单元测试 | 两项已知答案/篡改检测测试通过 |
| 历史证据 | 原 B 类、19 点参数扫描与 12 个公开格见证的独立核验通过 |
| CI 配置 | 三份 workflow YAML 本地解析通过，构建/测试使用 Workspace 与 locked；不再在固定版本验收任务中执行 cargo update |

CI 增加了例程测试、基准构建和分层检查；本轮没有运行托管的 GitHub Actions 或远程覆盖率任务。

完整实验文件与日志生成于 `target/workspace_migration/`，其哈希保存在结果记录中。
额外把这些原始数据归档到本地 `.migration-backups/2026-11-17-workspace-validation.tar.gz`，以免 cargo clean 清除唯一副本。

已保存的迁移时点还可执行：

```bash
python3 scripts/verify_workspace_layout.py --migration-snapshot
```

普通布局检查验证依赖方向与历史证据；上述额外选项还要求当前源码和锁文件匹配本次迁移清单。
后续有意修改代码后，迁移快照可以继续作为历史记录，不能把旧哈希不匹配误认为算法失败。

## 5. 迁移前检查点与后续工作

迁移前已将 76 个文件保存到本地：

```text
.migration-backups/2026-11-17-before-workspace.tar.gz
.migration-backups/2026-11-17-before-workspace-manifest.json
.migration-backups/2026-11-17-before-workspace.diff
```

检查点包含当时未跟踪的 PLT、测试、脚本和证据；不包含 target、Git 元数据、编辑器配置和论文 PDF。
该目录由 .gitignore 排除，保留在本机，不受 cargo clean 影响；它不替代后续的正常版本提交。

后续修改按职责进入相应 crate。底层采样器修复在 qfall-tools，PLT 选参和算法在 plt，上层功能开始开发时再建立新 crate。
本次源码迁移完成了组织和调用边界调整。**B3 仍未通过，PLT 参数批准入口继续拒绝，尚未实现的上层协议没有新增验收结论。**
