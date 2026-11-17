# PLT 与车联网协议研究工作空间

项目仓库：[wydnzdgz/P-R-S](https://github.com/wydnzdgz/P-R-S)。

本项目在 [qFALL-tools](https://github.com/qfall/tools) 基线上实现可配置门限的整数格 PLT。
基础采样器与 PLT 已分离为两个 Cargo package；后续上层功能与车联网协议按需增加独立 package。

| 层次 | 目录 | 负责内容 |
| --- | --- | --- |
| A：维护后的基础库 | [crates/qfall-tools](crates/qfall-tools/README.md) | MP12、GPV、陷门生成、压缩与工具函数，保留随机源及数值修复 |
| B：PLT | [crates/plt](crates/plt/README.md) | 通用 Setup、PTrapGen、PSampPre、Rec，门限配置与参数评估 |
| C、D：后续层 | 尚未创建 | 上层功能和车联网协议，实现时分别建 crate |

依赖方向为 `plt -> qfall-tools`。各层的源码、测试、示例和基准由各自 package 管理；根目录保留文档和跨层实验脚本。
原始基线与维护改动见 [基础库来源记录](docs/baseline_maintenance.md)。
Git 历史从本项目的起始提交记录，原作者声明与许可证保留；详见[历史整理说明](docs/project_history.md)。

## 构建与测试

使用 Linux/macOS 或 WSL，Rust 1.88+，以及构建 FLINT 所需的 m4、C 编译器和 make。
Python 独立验证脚本使用标准库；Workspace 布局检查需要 Python 3.11+。

```bash
# 分层测试
cargo test --locked -p qfall-tools
cargo test --locked -p plt

# 整体检查
cargo test --locked --workspace
cargo test --locked --workspace --examples
cargo fmt --all -- --check
cargo clippy --locked --workspace --all-targets -- -D warnings
python3 scripts/verify_workspace_layout.py
python3 scripts/check_repository_files.py
```

## PLT 实验

一套算法支持不同 t-of-k 参数：

```bash
cargo run --locked -p plt --example plt -- \
  --threshold 2 --participants 3 --output target/plt_2_of_3.json
python3 scripts/verify_plt_report.py target/plt_2_of_3.json

cargo run --locked -p plt --example plt -- \
  --threshold 3 --participants 5 --output target/plt_3_of_5.json

cargo bench --locked -p plt --bench benchmarks -- PLT3of5
```

Rust 调用路径已改为独立 crate：

```rust
use plt::{PltConfig, setup};

let config = PltConfig::new(3, 5)?;
let (public, master) = setup(&config)?;
```

基础库原有接口仍使用 `qfall_tools::primitive::psf` 等路径。
Workspace 共享一个 Cargo.lock 和 target 目录；依赖升级与固定版本的复现实验分开进行。

## 文档与当前范围

- [结构说明与分层规则](docs/project_architecture.md)
- [源码迁移记录、检查点与验证结果](docs/workspace_migration.md)
- [仓库上传规则与敏感/大文件检查](docs/repository_file_policy.md)
- [PLT 接口、尺寸及运行方法](docs/plt_2_of_3.md)
- [PLT 工具与车联网协议总体验收清单](docs/acceptance_checklist.md)
- [A 类工程验收记录](docs/plt_acceptance_A.md)
- [B 类采样器与参数审计](docs/plt_sampling_B.md)
- [参数运行范围及公开原像反例](docs/plt_parameter_exploration.md)

当前实现支持流程与工程实验。已测试的参数存在公开短原像反例，B3 未通过，协议参数批准入口继续拒绝。
源码分层不会改变这一结论。上层车联网协议尚未实现。

历史实验 JSON 保留原内容及源码指纹；迁移后的复测使用新的输出文件。
本项目保留原有 MPL-2.0 版权声明，见 [LICENSE](LICENSE)。
