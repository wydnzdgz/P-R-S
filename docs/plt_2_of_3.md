# 可配置整数格 PLT 原型

2026-11-17 Workspace 迁移说明：本文源码链接与复测命令已适配独立 `plt` package；历史 JSON、数值与源码指纹保持原样，见[迁移记录](workspace_migration.md)。

本版本用一套通用的 `Setup`、`PTrapGen`、`PSampPre`、`Rec` 实现 Albrecht 等人
《Partial Lattice Trapdoors: How to Split Lattice Trapdoors, Literally》第 5.2 节图 2 的流程。
[作者公开版](https://olapiha.github.io/assets/pdf/TTrap.pdf)。
门限规模由 `PltConfig` 决定，2-of-3、3-of-5、5-of-10 不存在各自独立的算法分支。
本说明沿用早期文件名以保留原有链接。后续验收状态见[总体验收清单](acceptance_checklist.md)。

## 使用统一接口

```rust
use plt::{PltConfig, setup, ptrap_gen, psamp_pre, rec};
use qfall_math::integer_mod_q::MatZq;

// 切换规模只需修改这两个参数；默认 n=1、q=257。
let config = PltConfig::new(3, 5).unwrap();
let (public, master) = setup(&config).unwrap();
let partials: Vec<_> = (0..config.participants())
    .map(|id| ptrap_gen(&public, &master, id).unwrap())
    .collect();
drop(master);
let width = partials.iter().map(|p| p.minimum_sigma_squared()).max().unwrap();
let target = MatZq::sample_uniform(config.target_rows(), 1, config.modulus());
let recovery_set = vec![0, 2, 4];
let shares: Vec<_> = recovery_set.iter()
    .map(|&id| psamp_pre(&public, &partials[id], &recovery_set, &target, width).unwrap())
    .collect();
let x = rec(&public, &shares.iter().collect::<Vec<_>>()).unwrap();
assert_eq!(public.a() * x, target);
```

`PltConfig::with_lattice(n, q, t, k)` 还允许配置整数格行数和素数模数。
构造器检查 `n>0`、`q` 为素数、`1<=t<=k<q` 及尺寸计算不溢出。
字段不可直接修改；公钥、主陷门、部分陷门和份额均绑定配置，避免跨配置混用。
配置校验保证代数兼容性，不检验密码学安全条件，也不保证任意巨大尺寸能在本机运行。

## 运行同一个实验程序

```bash
cargo test --locked -p plt --test plt_correctness

cargo run --locked -p plt --example plt -- --threshold 2 --participants 3 --output target/plt_config_2_of_3.json
python3 scripts/verify_plt_report.py target/plt_config_2_of_3.json

cargo run --locked -p plt --example plt -- --threshold 3 --participants 5 --output target/plt_config_3_of_5.json
python3 scripts/verify_plt_report.py target/plt_config_3_of_5.json

# 可选行数、模数和每个集合的目标数
cargo run --locked -p plt --example plt -- --n 2 --modulus 17 --threshold 2 --participants 3 --trials 3 --output target/plt_blocks.json
```

默认每个集合处理 8 个目标，包括确定性目标、零目标、边界目标和随机目标。
全部集合由同一组合迭代器生成：2-of-3 共 3 个集合/24 次恢复，3-of-5 共 10 个集合/80 次恢复。
组合枚举是惰性的；遍历全部组合的时间和生成的报告大小仍会随二项式系数增长。
`--trials` 控制每个集合的目标数，必须为正；`--help` 显示所有选项。

旧命令 `cargo run --locked -p plt --example plt_2_of_3 -- target/plt_2_of_3.json` 仍可使用，
它只是通用实验程序的兼容入口，不包含另一套算法；省略门限参数时默认 2-of-3。

报告保存配置、目标、各参与者的 U/C、整数 d/x 和范数。每次重新随机生成陷门和样本。
Python 工具从报告参数推导尺寸与合法集合，用独立整数运算核验，不为特定门限另写算法。
没有 `trials_per_set` 字段的历史报告按早期 8 个目标兼容读取。
报告包含实验用部分陷门材料，只用于本地实验；不是线上报文或正式协议日志。

## 参数、尺寸与接口

设 `ell=ceil(log2(q))`，统一推导 `m=n*(ell+1)`，参与者编号为 `0..k`，点为 `1..=k`。
素数模数及 `k<q` 保证这些点非零、任意两点之差可逆。

| 对象 | 通用尺寸 | 默认 2-of-3 | 默认 3-of-5 | 默认 5-of-10 |
| --- | --- | --- | --- | --- |
| A | nt × 2mt | 2 × 40 | 3 × 60 | 5 × 100 |
| C_j | n × m | 1 × 10 | 1 × 10 | 1 × 10 |
| U_j | 2mt × m | 40 × 10 | 60 × 10 | 100 × 10 |
| y | nt × 1 | 2 × 1 | 3 × 1 | 5 × 1 |
| z_j | n × 1 | 1 × 1 | 1 × 1 | 1 × 1 |
| d_j | m × 1 | 10 × 1 | 10 × 1 | 10 × 1 |
| x_j、x | 2mt × 1 | 40 × 1 | 60 × 1 | 100 × 1 |

公共接口位于 `plt`：

- `setup(&config)`：返回 `Result<(PublicParameters, MasterTrapdoor), PltError>`；生成公钥、检查陷门关系并返回主陷门。
- `ptrap_gen(&public, &master, id)`：生成一个参与者的部分陷门，对矩阵目标逐列独立采样并检查有理列秩。
- `decompose_target(&config, &set, &y)`：计算 `(V_T^-1 tensor I_n)y`，返回与输入 set 顺序一致的 t 个 n 行分量。
- `participant_image(&config, id, &matrix)`：构造 `(1,a_j,...,a_j^(t-1))^T tensor matrix`。
- `psamp_pre(&public, &partial, &set, &y, &sigma_squared)`：set 必须有恰好 t 个不同合法编号并包含当前参与者。
- `rec(&public, &share_refs)`：恰好收集 t 个不同参与者的份额，检查每份的模等式与范数，求和后再次检查完整原像。
- `verify_preimage(&public, &target, &x, &sigma_squared)`：供接收方检查外部整数原像的形状、范数及最终模等式。
- `verify_partial_preimage(&public, id, &set, &target, &x_j, &sigma_squared)`：检查部分原像的公共模等式和范数，不证明其属于秘密子格，也不认证发送者。
- `config.recovery_sets()`：惰性生成所有 t 子集，供实验和测试使用；算法本身不需要预先枚举全部集合。

恢复集合和份额使用动态切片/向量，不再使用固定长度二元数组。
集合和份额顺序可以变化；不同配置、公钥、目标、集合或宽度的份额会被拒绝。
输入不足或超过 t 份也会被拒绝；超过 t 个在线代理时，调用方先选定同一 t 人集合后采样。
验证方的目标和宽度必须来自批准的协议参数，不能接受对方自行扩大宽度。
这些一致性检查不是恶意参与者的份额验证协议。

## 正确性证据

同一套测试函数遍历所有配置的合法集合，检查：

1. `A U_j = v_j tensor C_j (mod q)`。
2. `C_j d_j = z_j (mod q)`。
3. `A x_j = v_j tensor z_j (mod q)`。
4. `A x = y (mod q)`。

同时检查精确整数等式 `x_j=U_j*d_j`、`x=sum_j x_j`，以及
`||x_j||² <= (2mt)*sigma²`、`||x||² <= t²*(2mt)*sigma²`。
这些范数界是概率性的。当前实现对发行列、部分原像和完整原像执行运行时范数拒绝，
越界返回 `NormExceeded`，不重采样、不放宽参数。通过验证的输出分布是成功条件下的分布；
拒绝率及其对协议/分布的影响仍需在 B 类验收中解释。

当前覆盖 2-of-3 和 3-of-5 的全部集合、每组 8 个目标；5-of-10 的全部 252 个集合、每组 1 个随机目标；
以及 t=1、t=k、不同素数和 n>1。独立已知答案覆盖二次以下 Vandermonde 分解和集合排列，
错误路径测试覆盖非法配置、错误集合大小、重复编号、错误模数、目标及混用份额。

## 底层采样适配

底层统一复用 `PSFPerturbation`，保留 gadget 底数 2、舍入参数 r=4 和三元陷门分布。
矩阵列数通过 `m_bar=columns-rows*ell` 显式设置。

初始 s 不再固定：从 64 开始，选择满足 `s² > 5*worst_case(||[R;I]||_F²)+1` 的最小 2 的幂。
对于三元 R，最坏 Frobenius 范数平方为 `gadget_columns*(m_bar+1)`。
因此调整门限后仍有保守正定余量；默认 n=1、q=257 时主陷门的 s 分别为 64（t=2）和 128（t=3、5）。
初始球形宽度为 r*s，这条规则保证采样计算的正定性，不是安全参数选取规则。

在线仍以精确有理矩阵构造 `Sigma_j=sigma²*(U_j^T*U_j)^-1`，
向底层传入 `Sigma_j/r²`，由原实现处理 r² 及 1/(2*pi) 换算。
平方宽度对应密度 `exp(-pi*x^T*Sigma^-1*x)`，不能直接等同于普通协方差。

共同在线平方宽度建议取各参与者下界的最大值：

```text
L_j = 4 r² ||U_j||_F² (5 ||B_j||_F² + 1),  B_j=[R_j;I]
sigma_squared = max_j L_j
```

该下界未替代论文的平滑参数、统计距离和安全性条件。底层仍含 f64 分解和连续采样，
采样器可信性及参数适用性继续按总清单核查。暂未实现环格版本。

## 通用性能基准

```bash
cargo bench --locked -p plt --bench benchmarks -- PLT2of3 --sample-size 30 --warm-up-time 1 --measurement-time 2
cargo bench --locked -p plt --bench benchmarks -- PLT3of5 --sample-size 30 --warm-up-time 1 --measurement-time 2
cargo bench --locked -p plt --bench benchmarks -- PLT5of10 --sample-size 30 --warm-up-time 1 --measurement-time 2
```

`crates/plt/benches/plt.rs` 的同一个函数依次读取三组配置，每组使用首个合法恢复集测量相同阶段。
在线测试复用该组的一套随机陷门和一个目标；单次调用仍重算分解和 Cholesky 因子。
分别测量 Setup、一个/全部部分陷门生成、每个参与者的采样、t 次串行采样加恢复、带检查恢复和纯整数求和。
编译、文件输出和网络传输不计时；未模拟在线并行。

每组参数记录在 `target/criterion/PLT{t}of{k}/context.json`，避免不同规模互相覆盖。
`--test` 可以仅执行基准代码检查，不产生正式耗时结论。
[2026-11-12 报告](plt_benchmark_2026-11-12.md) 是旧版固定 2-of-3 实现的历史记录；
不能把它当成此次通用实现或其他门限规模的新测量。

## A 类工程验收与失败记录

本轮证据与边界见 [A1–A4 验收记录](plt_acceptance_A.md)。

```bash
cargo test --locked -p plt --lib
cargo test --locked -p plt --test plt_correctness
cargo run --locked -p plt --release --example plt_acceptance -- --rounds 1000 --keysets 3 --output target/plt_acceptance.json
```

公共 `SamplingProfile` 记录实际使用的陷门分布、尺寸、舍入参数和初始宽度，
Setup/PTrapGen 会检查 `A[R;I]=G` 和 `AU_j=v_j tensor C_j`。
单次操作的采样失败、范数越界和数值错误向调用方返回 `Result`；Setup 也已采用此接口，
调用方需使用 `?` 或显式处理错误。底层 PSF 的旧接口为兼容仍保留 panic 语义，
PLT 改用同一 MP12 核心的 `try_trap_gen`、`try_compute_sqrt_sigma_2` 和 `try_samp_p`。

重复实验记录每个阶段的尝试数、错误原因、范数失败、恢复成功/失败和跳过次数，
不会把失败后重采样伪装为成功。不涵盖内存耗尽、进程被终止或 C 库主动 abort。
这份记录检查工程不变量与失败率，不是 B4 的分布统计验证或 B3 的安全参数验收。

## B 类审计与参数限制

[B 类记录](plt_sampling_B.md) 包含随机源修复、条件数/换算检查、重复采样与小维离散参照。
生产 MP12 路径已改为每次操作读取 OS 熵、256 位种子、ChaCha12 和 StandardNormal，
关键路径不再使用 probability/Xorshift。采样相关依赖和 Cargo.lock 已固定。

`partial.sampling_diagnostics(&sigma_squared)` 可查看当前数值诊断。
`public.assess_parameters(&sigma_squared)` 返回参数评估；当前实验参数允许公开求解通过原像检查，
不能作为独立门限不可伪造性的依据。`public.require_protocol_parameters(&sigma_squared)` 对当前实例返回拒绝。
B1/B2/B4 有限工程检查通过没有改变这一结论。

```bash
cargo run --locked -p plt --release --example plt_sampling_audit -- target/plt_sampling_B.json
python3 scripts/verify_sampling_audit.py target/plt_sampling_B.json
```

新结果不覆盖 A 类历史快照；变更随机源和数值守卫后的性能需重新测量。
