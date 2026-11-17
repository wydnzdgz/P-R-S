# 项目开发历史

日期：2026-11-17。

本仓库从当前工程状态建立独立的起始提交，后续只记录 P-R-S 项目的开发。
源码、参数、测试与验收数据继续保留。原基础库的版权头、LICENSE、来源和维护改动说明保留。

## 原项目来源

- 原始仓库：<https://github.com/qfall/tools>
- 基线提交：[851bacc722cfb52f664127c441f5f658729213fd](https://github.com/qfall/tools/commit/851bacc722cfb52f664127c441f5f658729213fd)
- 本项目维护内容：[基础库来源与维护改动](baseline_maintenance.md)

原先继承的 Git 提交不再是 main 的祖先。本地的上游 dev 分支、upstream 远程引用以及上游版本标签也一并移除，避免继续在开发图中显示。
origin 继续指向 `https://github.com/wydnzdgz/P-R-S.git`。

## 本地恢复备份

整理前的历史完整打包到：

```text
.migration-backups/2026-11-17-before-project-history-reset.bundle
.migration-backups/2026-11-17-before-project-history-reset.json
```

bundle 已用 `git bundle verify` 验证；JSON 记录整理前的提交、引用、远程地址和文件哈希。
该目录被 .gitignore 排除，不上传到 GitHub。需要检查旧历史时，可以从 bundle 克隆到另一目录查看。

更新远程使用带精确预期提交的 `--force-with-lease`；如果远程在检查后发生变化，推送会拒绝覆盖。
这次整理只改变开发历史与来源说明，历史实验 JSON 及算法实现保持原样。
