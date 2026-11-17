# PLT 2-of-3 性能实测（2026-11-12）

CPU：Intel Core i7-12700，Microsoft 虚拟机环境；rustc 1.99.1 (58cf8f9c3 2026-10-13)。
使用 Cargo bench 优化构建（保留调试信息）。每项预热 1 秒、30 个统计样本，目标测量时间 2 秒；部分项目实际采集稍长。
下表使用 Criterion 的 mean 点估计和 95% 置信区间，统一以微秒表示。

参数：q=257，n=1，m=10，t=2，k=3，底数 2，r=4，主采样宽度 256。
在线基准固定一套随机部分陷门，sigma²=18697451264，sigma≈136738.624，目标 y=(12,19)。
在线宽度按陷门自适应选择，下次运行可能不同。离线基准每轮重新随机生成陷门。

| 阶段 | 均值（µs） | 95% 置信区间（µs） |
| --- | ---: | ---: |
| `online_sequential_T01` | 433.396 | [428.556, 438.513] |
| `online_sequential_T02` | 434.659 | [431.021, 438.240] |
| `online_sequential_T12` | 436.320 | [432.559, 440.028] |
| `psamp_pre_T01_party0` | 214.914 | [212.204, 217.763] |
| `psamp_pre_T01_party1` | 212.078 | [209.806, 214.375] |
| `psamp_pre_T02_party0` | 213.176 | [211.498, 214.869] |
| `psamp_pre_T02_party2` | 216.120 | [213.750, 218.876] |
| `psamp_pre_T12_party1` | 214.336 | [212.970, 215.777] |
| `psamp_pre_T12_party2` | 214.378 | [212.365, 216.545] |
| `ptrap_gen_all_three` | 15396.396 | [15309.234, 15495.281] |
| `ptrap_gen_one` | 5053.930 | [5023.005, 5090.297] |
| `rec_checked` | 0.502 | [0.499, 0.505] |
| `rec_integer_sum_only` | 0.134 | [0.133, 0.135] |
| `setup` | 408.363 | [404.161, 413.373] |

`psamp_pre` 包括目标分解、非球形参数矩阵构造、Cholesky 分解、MP12 采样和整数矩阵乘法。
`online_sequential` 是两个参与者依次采样后调用 rec，未测并行和网络。
`rec_checked` 包含 API 的上下文及 A*x=y 检查；`rec_integer_sum_only` 仅做两个 40 维整数向量求和。
未测峰值内存和通信开销；未模拟安全参数下的大尺寸问题。在线测量只覆盖一套随机陷门和一个目标。
这些结果说明当前小参数实现的运行开销，不表示方案达到任何安全级别。

复现命令：

```bash
cargo bench --locked --bench benchmarks -- PLT2of3 --sample-size 30 --warm-up-time 1 --measurement-time 2
```

原始结果：`target/criterion/PLT2of3/`；参数记录：`target/criterion/plt_context.json`；
本次均值汇总：`target/criterion/plt_summary.json`。
后续运行会更新 target 下的结果，本文件记录本次测量。
