# B 类采样器与参数审计记录

2026-11-17 Workspace 迁移说明：本文源码链接与复测命令已适配独立 `plt` package；历史 JSON、数值与源码指纹保持原样，见[迁移记录](workspace_migration.md)。

日期：2026-11-13。结论：**B1 源码核查及随机源修复完成，B2 数值与接口检查在测试范围内通过，B4 有限筛查通过；B3 当前参数不通过。B 类整体没有通过。**

本轮没有把“修复随机源”和“统计实验通过”当作当前 PLT 参数不可伪造的证据。
参数反例足以否定当前实例作为独立门限不可伪造原语的安全依据；上层若采用其他授权机制，需要在 C 类中单独说明。

2026-11-17 后续：[参数运行范围评估](plt_parameter_exploration.md) 已按“PLT 自身承担门限不可伪造性、先不设安全位数”开展。
19 个参数点中 12 点完成算法流程，均找到公开 LLL 原像反例；B3 仍未通过。本页原有数据保留为历史审计证据。

## 1. B1：调用链与实际修复

### 发现的问题

原 MP12 连续扰动路径：

```text
PLT -> MP12 -> qfall-math::sample_d_common_non_spherical
    -> sample_gauss_same_center -> probability::source::default(u64)
    -> random::default -> Xorshift128Plus
```

本地锁定版本的 `random` 源码明确将默认源定义为 Xorshift128Plus，并从一个 u64 构造初态。
即使那个 u64 来自 CSPRNG，这条路径也只有 64 位种子空间，而且生成器本身没有密码学不可预测性保证。
因此原路径不能作为本项目生产采样随机性的依据。这是本次审查发现并修复的实现问题。

### 修复后的关键路径

| 阶段 | 当前来源与处理 | 代码 |
| --- | --- | --- |
| 系统熵 | 每次操作读取 32 字节 OS 熵；失败返回错误，无固定种子/弱随机回退 | [randomness.rs](../crates/qfall-tools/src/sample/randomness.rs) |
| A_bar | 新 OS 种子交给 qfall-math 的 seeded 均匀矩阵接口 | [checked.rs](../crates/qfall-tools/src/primitive/psf/mp_perturbation/checked.rs) |
| 三元 R | 两个独立新种子分别生成 X、Y，再计算 X-Y；错误向上传递 | [trapdoor_distribution.rs](../crates/qfall-tools/src/sample/g_trapdoor/trapdoor_distribution.rs) |
| 连续扰动 | 新 OS 种子 -> rand 0.10.3 StdRng（该版本为 ChaCha12）-> rand_distr StandardNormal | [checked.rs](../crates/qfall-tools/src/primitive/psf/mp_perturbation/checked.rs) |
| 随机舍入 | 另取新 OS 种子，调用 qfall-math seeded 离散舍入接口，内部逐坐标派生子种子 | 同上；依赖中的 mat_q/rounding.rs |
| Gadget 条件离散采样 | 另取新 OS 种子，调用预计算 GSO 的 seeded SampleD，内部逐坐标派生子种子 | [mp_perturbation.rs](../crates/qfall-tools/src/primitive/psf/mp_perturbation.rs) |

MP12 生产路径不再使用上述 probability/Xorshift 连续分支。该依赖仍被 qfall-math 的其他 API 使用；
本轮没有审计或改写其所有函数，也不为任意自定义 TrapdoorDistribution 背书。

本项目不缓存跨操作的用户态 RNG；正常进程重启或 fork 后，后续进入采样操作会重新读取 OS 熵。
不宣称可对抗任意进程内存复制、VM 熵故障或恶意 OS。测试覆盖熵读取失败无回退、并发调用的新种子不重复；
这些测试检查接线与错误处理，不证明 OS 熵质量。

依赖中的整数离散高斯是有限区间拒绝采样：默认 tailcut=6，接受权重按 `exp(-pi*(x-c)^2/s^2)` 以 f64 计算。
GSO 通过精确有理内积更新中心，根号和指数包含数值近似。这些截断/舍入误差没有被包装成“精确无限支持高斯”。
不在运行中修改依赖的全局 TAILCUT；没有完成常数时间、侧信道或内存清零审计，后者属于部署要求。

### 版本与证据

直接固定 qfall-math=0.2.1、flint3-sys=3.6.1、rand=0.10.3、rand_distr=0.6.0，并保留 Cargo.lock。
依赖来源及审查文件指纹见 [dependency audit JSON](plt_B_dependency_audit.json)。
测试的固定输入只用于数值/统计参照，不作为生产随机种子。源码哈希与工具链记录在结果快照中。

## 2. B2：矩阵、高斯约定与数值范围

对整数格，`U^dagger=U^T`。当前实现使用精确有理 Gram 逆构造
`Sigma=sigma²*(U^T U)^-1`，向现有 MP12 输入 `Sigma/r²`。
连续扰动因子满足：

```text
L L^T = [Sigma - r²*5*B B^T - r²*I] / (2*pi),  B=[R;I], base=2
```

独立数值诊断将 `2*pi*L L^T + r² I + r²*5*B B^T` 重构回目标 Sigma，
避免只检查局部 Cholesky 因子却漏掉 r² 或 2*pi 换算。
同时精确检查 `Gram*Sigma=sigma² I`，这对应 U 映射下的欧氏高斯权重关系。

运行时额外计算扰动矩阵的无穷范数条件数估计，拒绝
`kappa_inf*epsilon_f64>1e-8`、非有限数、非正主元和超限重构残差。
该估计通过 Cholesky 三角求解获得；审计诊断则使用独立的精确有理逆后计算范数。
条件数很大的正定矩阵也会被拒绝，不能仅因为它在代数上正定就接受数值结果。

单元/集成测试覆盖默认多种门限、n>1、不同素数、良态/近奇异/尺度失衡/非对称/非正定/非有限情形。
这是对所测工程范围的检查，不是任意精度误差证明；新参数仍需重新评估数值范围和 B3 条件。

## 3. B3：当前参数明确不通过

### 可公开构造的反例

Setup 已验证 A 满行秩。对模素数的满行秩 A，仅用公钥就可以高斯消元求解 Ax=y：
选主元、令自由变量为零，随后将坐标取中心剩余类。
所得整数向量至多 nt 个非零坐标，因此：

```text
||x_public||² <= nt*floor(q/2)^2
```

当前接受界为 `beta²=t²*(2mt)*sigma²`，远大于这个公开上界。
本次三组配置均实际构造了不依赖主陷门/部分陷门/份额的公开解，并调用 `verify_preimage` 得到接受结果。
这说明“满足同余式且在当前范数界内”不能作为已获得门限授权的证据。
它没有攻击或否定满足条件的原论文构造，而是否定当前玩具参数的直接安全用途。

| 配置 | q,n,m | 公开解范数平方上界 | 本次公开解通过验证 |
| --- | --- | ---: | --- |
| 2-of-3 | 257,1,10 | 32,768 | 是 |
| 3-of-5 | 257,1,10 | 49,152 | 是 |
| 5-of-10 | 257,1,10 | 81,920 | 是 |

反例函数只接收 A 和 y；实验在求解前释放所有陷门对象。
快照保存 A、y、公开解及验证结果，Python 用独立整数运算重新核验。

### 与论文参数条件的对应

以下是与 [作者公开版](https://olapiha.github.io/assets/pdf/TTrap.pdf) 第 3.5 节和 Table 2 的必要对照，
不是声称已完成全部理论参数实例化。Table 2 的对应还依赖采用其 Lemma 13 底层方案。

| 条件 | 当前对应 | 结论 |
| --- | --- | --- |
| m=n*(ceil(log2(q))+1)，t<=k，点差可逆 | 按配置推导，素数 q，k<q | 代数条件满足。 |
| Table 2 中 S=SCor，smax(SCor)<=q^(1/g)/sqrt(m) | 整数格 g=1，q/sqrt(m)≈81.27；当前发行宽度 256 或 512 | 已知必要上界不满足。代码记录平方上界 66049/10，避免浮点比较。 |
| Table 2 的 m-k 条件及含 1/(m-k) 的表达式 | 5-of-10 时 m=k=10 | 不满足非退化条件；其他组也没有实例化 lambda 对应的渐近条件。 |
| 底层分布和公钥接近均匀 | 仍使用三元 R | 没有据 Lemma 13 认定引用适用。 |
| 平滑参数、误差预算、sigma 的理论下界 | 当前下界保证数值正定 | 未完成；不能将 omega 表达式随意替换成常数。 |
| beta>=t*sigma*sqrt(2mt) | 实现采用该范数界并运行时检查 | 功能条件满足，不推出困难性。 |
| 具体困难问题和攻击代价 | 没有安全参数族/目标强度评估，且有上述公开反例 | 当前协议安全实例拒绝。 |

`public.assess_parameters(sigma_squared)` 返回这些检查和未建立的条件；
`public.require_protocol_parameters(sigma_squared)` 是可选的明确拒绝入口，当前不批准任何未经评估的实例。
普通流程实验继续可运行，不会把参数构造器成功误当成协议安全批准。
收紧 beta、只增加 n 或通过统计测试，都不自动完成 B3。

## 4. B4：预先登记的有限实验

实验计划先写入 [plt_B_experiment_plan.md](plt_B_experiment_plan.md)，再运行；本次没有按结果修改阈值。
每组 3 套密钥，每套按规则选一个参与者与合法集合，对零目标、固定非零目标、固定随机目标各采样 1,000 次。
使用 `sqrt(2*pi)/sigma * chol(U^T U)^T d` 做连续近似下的粗略尺度筛查；非零陪集不强制零均值。

数值诊断覆盖 54 个实际辅助采样矩阵，`Gram*Sigma=sigma² I` 全部精确成立。

| 指标 | 最小值 | 最大值 |
| --- | ---: | ---: |
| `gram_condition_inf` | 4.02457 | 20.3963 |
| `covariance_condition_inf` | 4.02457 | 20.3963 |
| `perturbation_condition_inf` | 4.03036 | 20.5201 |
| `total_squared_width_relative_residual` | 1.14224e-16 | 3.66505e-16 |

27 个重复采样组共 27,000 次 PLT 部分采样，错误数为 0，全部通过预先登记的有限筛查：

| 指标 | 本次观测 | 预先登记标准 |
| --- | ---: | --- |
| 坐标方差范围 | [0.884086, 1.140851] | [0.7,1.3] |
| 零目标最大绝对坐标均值 | 0.094176 | <=0.15 |
| 最大绝对坐标相关 | 0.106263 | <=0.2 |
| 最大绝对 lag-1 相关 | 0.087406 | <=0.2 |
| 完整向量重复数 | 0 | 本次大宽度实验要求为 0 |

小维目标 0：10,000 个样本；四个投影的最大 CDF 差 0.007281（阈值 0.04）；截断相对质量上界 3.92e-10（要求 <1e-8）。

小维目标 1：10,000 个样本；四个投影的最大 CDF 差 0.011419（阈值 0.04）；截断相对质量上界 3.92e-10（要求 <1e-8）。


小维参照使用 q=3、三个坐标和固定非对角 Sigma，枚举有限盒中满足模条件的整数点，
使用离散权重构造 CDF；这与直接套连续正态的 CDF 不同。
遗漏质量界以无条件高斯上界控制，并除以盒内分区函数下界；exp 和累加仍有浮点误差，未称作符号精确结果。
四个投影通过不代表完整联合分布通过，更不代表长期秘密不泄露。

本次共 47,000 次样本：27,000 个 PLT 部分原像和 20,000 个底层 SampPre 样本。
没有观察到所检查的明显异常。成功输出的条件分布、拒绝策略和安全统计距离仍不能由有限实验认证。

## 5. 文件、命令和状态

- [结果快照](plt_sampling_B_results.json)：数值指标、投影参照/样本直方图、公开解反例、版本和源码 SHA-256。
- [审计程序](../crates/plt/examples/plt_sampling_audit.rs)：同一配置化 PLT 实现的实验驱动。
- [独立离散参照](../crates/plt/examples/support/sampling_reference.rs)：小维盒枚举。
- [核验脚本](../scripts/verify_sampling_audit.py)：重算 CDF 差、截断界和公开解的模等式/范数。
- [数值与参数测试](../crates/plt/tests/plt_sampling.rs)、[参数诊断接口](../crates/plt/src/assessment.rs)。

```bash
cargo test --locked --workspace
cargo fmt --all -- --check
cargo clippy --locked --workspace --all-targets -- -D warnings
cargo run --locked -p plt --release --example plt_sampling_audit -- target/plt_sampling_B.json
python3 scripts/verify_sampling_audit.py target/plt_sampling_B.json
```

本次全项目 157 个测试通过（97 单元、9 集成、51 文档），格式、Clippy 和文档构建检查另行完成。
审计命令成功退出表示工程实验完成；JSON 的 `stage_b_fully_passed` 和 `protocol_parameters_passed` 仍为 false，
不能只凭退出码把 B3 当作通过。数据快照的源码指纹是本轮版本，A 类历史快照不被覆盖。

若 PLT 自身承担门限不可伪造性，下一步必须明确协议攻击模型与参数要求后继续 B3；
若它仅为流程/性能工具，独立的代理认证/授权机制仍需在 C 类中验证，不能将当前 PLT 计入安全保证。
