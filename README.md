# MoonCopula（月依存）

纯 MoonBit 的二维 Archimedean Copula 建模库：把配对观测转换为秩分布，拟合并比较依赖模型，计算联合事件概率，再将条件场景还原为原始量纲。版本 `0.1.0`，模块 `CYang-dep/mooncopula`，MIT 许可证。

它面向需要联合超越概率的环境分析者、需要相关寿命场景的可靠性工程开发者，以及需要可解释模型比较的风险工具开发者。三个完整、固定输入的示例见下文；合成数据仅用于展示计算流程。

## 已实现能力

- 独立、Clayton、Gumbel、Frank 四类二维模型；0/90/180/270 度旋转。
- 闭方形 CDF、开方形密度/对数密度、条件 CDF、带残差及迭代次数的条件反解。
- 正向/逆向 Rosenblatt 变换；模拟均匀数由调用方提供，无隐藏随机状态。
- 稳定平均秩、顺序秩及拒绝并列策略，Kendall tau-b、Spearman rho、有限阈值尾部统计。
- 经验联合分布与矩形概率；无并列数据的 Checkerboard 平滑和均匀边缘。
- 有界伪极大似然拟合、收敛/边界标志、AIC/BIC、显式候选集比较与 AIC 权重。
- 联合/任一超越、混合象限概率、条件超越、重现期、渐近同尾系数。
- 经验边缘的阶梯广义逆与独立的 type-7 线性分位数，原始量纲场景生成。
- 网格差异、Rosenblatt 残差、KS/Cramer-von Mises 描述统计、阈值事件校准与 Wilson 区间。
- 调用方提供索引计划的成对重采样：逐次重新排名和拟合、百分位汇总、边界拟合计数；失败的重复不会被悄悄丢弃。

不提供高维/vine、Gaussian/t、自动边缘分布选择、缺失值清洗、时间序列建模、隐含随机数发生器、显著性 p 值或因果推断。与生态中面向金融 Gaussian/t 相关随机向量的 MonteCarlo 包有不同入口和核心任务；查重依据与限度见 [选题与查重](docs/competition/duplicate-check.md)。

## 环境和安装

固定验证工具链：`moonc v0.10.18+6c83a7f62`，`moon 0.1.20261028 (e1e1f15)`，高于验收要求的 0.10.16。使用官方 MoonBit 工具链；JS 运行需要 Node.js。支持 `wasm-gc`、`wasm`、`js` 以及 Linux `native`，CI 针对四个后端分别执行检查、构建、测试和全部示例。Windows native 需要另外安装 MSVC C++ 工具和 Windows SDK；本机缺少这些组件，未声称在 Windows 上运行 native。

在一个新的目录创建消费者，显式指定自己的模块身份：

```sh
moon new copula-consumer --user guest --name copula-consumer
cd copula-consumer
moon add CYang-dep/mooncopula@0.1.0
```

在 `cmd/main/moon.pkg` 中导入库：

```moonbit
import { "CYang-dep/mooncopula" }
options("is-main": true)
```

将 `cmd/main/main.mbt` 替换为下面的快速开始，然后执行 `moon check`、`moon build`、`moon run cmd/main`。本库运行时仅依赖 MoonBit core，不需要 Python、网络、token 或本地登录。

## 快速开始

```moonbit
///|
fn main raise {
  let model = @mooncopula.Copula::new(@mooncopula.Clayton, 2.0)
  println(model.cdf(0.5, 0.5))
  let inverse = model.inverse_conditional(0.3, 0.9)
  assert_true(inverse.residual <= 1.0e-8)
  println(inverse.value)
}
```

首行 CDF 约为 `0.3779644730092272`。反解同时提供 `value`、`residual`、`iterations`，不能只看结果值忽略收敛情况。示例使用 `fn main raise` 传播错误，调用系统可以捕获 `CopulaError`，区分非法输入、资源预算和数值失败。

## 三个可运行工作流

从源码复现：

```sh
git clone https://github.com/CYang-dep/mooncopula.git
cd mooncopula
git checkout v0.1.0
moon run examples/environment
moon run examples/reliability
moon run examples/risk
```

| 用户与输入 | 示例的处理 | 可观察结果 |
|---|---|---|
| 环境分析：12 组模拟降雨/流量 | 平均秩、Gumbel 拟合、联合超越、经验对照、网格诊断 | theta 约 9.80595，阈值 0.8 的联合超越约 0.187031，经验频率 1/12 |
| 可靠性开发：6 组模拟寿命、4 对固定均匀数 | Clayton 逆变换、经验寿命边缘恢复、残差检查 | 场景 `(100,130)`、`(240,410)`、`(400,360)`、`(500,410)` |
| 风险工具：8 组配对观测 | 11 候选 AIC 比较、校准、Rosenblatt 诊断、4 个固定重采样计划 | 最优 Frank 正分支 AIC 约 -6.31681，联合预测约 0.181204，输出 4 个重复结果 |

浮点末位和近似并列模型的顺序可因平台不同变化。风险示例的四个计划仅验证执行机制，不能据此声称百分位区间具有实际覆盖率；实际应用应由调用方提供符合数据抽样设计的充分重复次数。

## 核心 API 与语义

```moonbit
let uniforms = @mooncopula.pseudo_observations(raw)
let fitted = @mooncopula.fit(@mooncopula.Clayton, uniforms)
assert_true(fitted.converged)
let joint = fitted.model.risk(0.8, 0.9)
let scenarios = fitted.model.scenarios(raw, supplied_uniforms)
```

这段为主函数内的接口片段，`raw` 为配对原始观测，`supplied_uniforms` 为调用方提供的开方形均匀对。完整可运行变量与错误处理见示例。完整公开签名见 [接口文件](pkg.generated.mbti)。

- `ranks` 已经返回秩除以 `n+1` 的值，勿重复归一化。默认平均秩；并列可显式拒绝或按输入次序打破。常数坐标的相关系数无定义，会抛出错误。
- `cdf` 接受 `[0,1]²` 并有精确边缘；密度/似然/模拟输入要求 `(0,1)²`。`conditional(u,v)` 为 `∂C(u,v)/∂u`，不是 `C(u,v)/u`。
- `risk` 的下侧为 `<=`，上侧为 `>`；`rectangle(lo_u,hi_u,lo_v,hi_v)` 对应 `(lo_u,hi_u]×(lo_v,hi_v]`。经验 CDF 有离散边缘；只有无并列的 Checkerboard 平滑保证连续均匀边缘。
- Clayton 参数 `[0.01,30]`，Gumbel `[1,20]`，Frank `[-30,-0.01]∪[0.01,30]`；精确独立使用 `independence()` 或 `Independence,0`。
- Frank 拟合默认正分支；负分支须显式给出负上下界，不能跨零拟合。`standard_candidates()` 明确包含两个分支与旋转模型。AIC 权重是相对拟合分数，不是模型正确概率。
- 拟合先扫描 25 个点，再在最优邻域黄金分割；保留端点最优、报告预算耗尽，不保证全局最优。`at_boundary` 提示结果受参数域约束。二维连续边缘被经验排名替代后的似然属于伪似然。
- `return_period(p,period_length)` 返回 `period_length/p`；零概率返回 `None`，解释为平均事件周期需等长、独立且平稳的观测期假设。
- 诊断均为描述量；同样本拟合、排名与序列依赖会影响显著性或区间覆盖。Wilson 区间描述观测 Bernoulli 频率，未包含参数估计不确定性。Bootstrap 百分位并不自动保证覆盖率。
- 每组最多 4096 个观测、候选最多 32 个、网格最多 63×63、重采样 2–128 次且最多 262144 个总抽取。Kendall/tie 统计与经验网格的最坏复杂度为二次；调用者应按工作量选择预算。
- 不接受 NaN、无穷、缺失值；接近奇点时密度可下溢，反解可显式报告数值失败。`log_density` 优先于 `density`，不能将测试容差外推为全域精度保证。

## 开发与验证

```sh
moon info
moon fmt --check
moon check --deny-warn
moon build --deny-warn
moon test --deny-warn
moon test --deny-warn --target wasm
moon test --deny-warn --target js
moon test --deny-warn --target native
```

31 个测试用例包含独立数值预期、边界、并列、旋转、预算、拒绝无效输入、所有核心流程及函数不变量。statsmodels / SciPy / 80 位 Decimal 的独立向量可重新生成：

```sh
python -m pip install -r scripts/oracle-requirements.txt
python scripts/generate_oracle.py --check
python scripts/audit_source.py
```

Python 仅用于开发验证，要求 3.12。数据、版本与容差来源见 [独立预期](docs/validation/independent-oracle.md)。测试通过不代表数学全域精度、统计假设或无 Bug 的证明。CI 定义见 [.github/workflows/ci.yml](.github/workflows/ci.yml)。

## 结构、许可和 AI 使用

根目录 13 个生产 `.mbt` 文件按概率引擎、变换、秩统计、经验分布、风险、拟合、选择、边缘、诊断、重采样、校准拆分；测试与生成测试向量独立，三个 `examples/` 目录为真实消费者入口。设计见 [架构](docs/architecture.md)，发布变更见 [CHANGELOG](CHANGELOG.md)。

本项目为原创方程与算法实现，采用 [MIT](LICENSE)。没有复制 statsmodels、SciPy 或相邻生态包的实现。依赖、参考来源和许可证边界见 [THIRD_PARTY.md](THIRD_PARTY.md)。AI 参与选题研究、实现、测试与文档，参见 [AI 使用说明](AI_USAGE.md)。输入与状态的边界见 [SECURITY.md](SECURITY.md)。
