# 基础库来源与维护改动

日期：2026-11-17。

`crates/qfall-tools` 是本项目维护的基础库，来源为 `https://github.com/qfall/tools.git`。
基线对照为上游提交 [851bacc](https://github.com/qfall/tools/commit/851bacc722cfb52f664127c441f5f658729213fd)（Update testing pipeline）。
当前仓库使用独立的项目开发历史；上游历史另存本地备份，不再作为 main 的祖先或本地分支展示，见[历史整理记录](project_history.md)。
当前工作空间中的 `plt` 是独立新增层，不属于上游基础库。

| 范围 | 相对原始基线的维护内容 |
| --- | --- |
| `primitive/psf/mp_perturbation.rs` 与 `mp_perturbation/checked.rs` | 增加 fallible MP12 接口，原接口委托共享核心；传播采样错误并检查矩阵/数值条件 |
| `sample/randomness.rs` | 关键操作获取新 OS 种子，熵失败返回错误 |
| `sample/g_trapdoor/trapdoor_distribution.rs` | 三元采样使用新种子；增加 fallible 采样入口 |
| `sample/g_trapdoor/gadget_classical.rs` | 传播陷门分布采样错误 |
| 依赖配置 | 锁定采样相关版本；继承根 Workspace 的版本约束与 Cargo.lock |

随机源审计依据见 [B 类记录](plt_sampling_B.md)。这些修改是先前审计形成的维护工作，本次源码迁移保留其行为。
基础库的压缩、工具函数、GPV、MP12 及其原有测试继续归属于本 crate；PLT 通过公开采样接口使用它。

迁移后基础库保留 `qfall_tools::primitive::psf` 等原有路径，不再登记 `primitive::plt`。
PLT 调用方改用 `plt::{...}`。基础库不依赖 PLT，也不包含协议代码。

新增底层修复时，在此记录来源、原因、接口影响和验证证据；保持已有版权头和根 LICENSE。
源码迁移的旧路径、新路径及内容哈希见 [迁移清单](workspace_migration_manifest.json)。
