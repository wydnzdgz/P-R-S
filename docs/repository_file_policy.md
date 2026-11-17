# 仓库文件上传规则

日期：2026-11-17。

## 保留与排除

- 提交源码、配置、Cargo.lock、文档，以及已审查的小型验收 JSON。
- 忽略环境凭据、私钥、云认证配置、编辑器设置、本地源码备份、构建缓存、日志和临时输出。
- 论文 PDF、压缩包、原始数据集和实验运行目录保留在本地。`experiments/` 下的配置仍可提交，不整体忽略实验目录。
- `.env.example`、`.env.sample`、`.env.template` 可以作为无真实凭据的配置样例；它们仍需通过内容检查。
- 默认单文件上限为 **1 MiB（1,048,576 字节）**。Git 不支持在 .gitignore 中按文件大小匹配，因此由独立检查脚本执行。

当前论文 PDF 已从首次发布提交中移除，文件仍留在本地。
当前保留的最大验收文件为 `docs/plt_parameter_exploration_results.json`，431,793 字节。
既有验收 JSON 和算法源码没有因上传规则而删除。

## 推送前检查

```bash
# 检查索引，即下一次提交实际包含的文件
python3 scripts/check_repository_files.py

# 检查 HEAD 可达的历史文件，避免旧提交仍带有不应上传的内容
python3 scripts/check_repository_files.py --history

# 检查典型敏感路径、凭据识别与公共数据误报
python3 -m unittest discover -s scripts -p 'test_repository_files.py'
```

脚本检查被 .gitignore 排除但仍被 Git 跟踪的路径、大小上限、私钥头、常见服务令牌和含凭据的 URL。
发现可疑内容时只打印文件名和类别，不回显匹配值。这是有限格式筛查，不能保证识别所有自定义秘密。
CI 使用完整历史检出并执行相同检查。

`.gitignore` 对已跟踪文件不生效。移除追踪时应保留本地文件；如果文件已经进入历史，必须同时检查历史。
本次首次推送前整理了尚未发布的本地提交，避免把论文 PDF 随旧提交一起上传，没有使用强制推送。

GitHub CLI 的本机凭据保存在项目外的用户配置目录；当前环境没有系统凭据库，CLI 使用本地凭据文件，检查时权限为 0600。
没有将凭据复制到项目或远程 URL 中，Git 的认证通过凭据助手处理。

新增 CI 文件检查属于迁移之后的发布配置改动；[迁移清单](workspace_migration_manifest.json) 保留原迁移时点的哈希。
日常仍可运行 `python3 scripts/verify_workspace_layout.py`；额外的迁移快照匹配选项只用于对照原迁移时点。
