# MoonCopula 选题、需求与语义查重

日期：2026-11-11。采用固定九项工程验收清单及 builder 附加规模要求，不代表官方接受参赛资格。

## 需求门禁
成熟参照是 [statsmodels copula API](https://www.statsmodels.org/stable/distributions.html#copulas) 与 [copulae](https://copulae.readthedocs.io/en/latest/api_reference/copulae/archimedean/index.html)。它们说明跨生态存在将边缘分布与依赖结构分离的真实工具需求；不声称采访过用户。MoonBit 当前已检查的金融向量模拟接口不能直接输出联合超越概率、Archimedean 条件分布与样本秩拟合。这里的用户需求是有成熟工具支持的推断，非既有生产用户或安全认证。

三类用户与闭环：环境分析员输入两列降水/流量观测→秩变换与候选拟合→联合超越概率和条件情景；可靠性工程师输入传感器/寿命配对样本→经验边缘与条件逆变换→保持边缘尺度的相关情景；风险开发者输入均匀样本或参数→似然比较/尾依赖/矩形概率→可复算的联合风险结果。库负责数学模型，不负责因果推断、时间依赖或业务决策。

核心完成标准：独立/Clayton/Gumbel/Frank、90/180/270度旋转、CDF/密度/条件与逆条件、秩与经验 copula、参数拟合/模型比较、经验边缘模拟和明确错误/预算。每项必须有不变量、参考方程或独立实现对照。无 vine、高维 Gaussian/t 模型、实时服务或图形界面。

## 三领域候选与评分（每列1–5）
|候选|广度|需求|复用|低重合|深度|验证|可行|维护|结论|
|---|---|---|---|---|---|---|---|---|---|
|系统发育树/简约重建|4|5|4|1|5|5|4|4|硬失败，BioSeqs 直接覆盖|
|调查抽样权重/校准|5|5|5|1|4|5|4|4|硬失败，moon-causal/moon-rake 直接覆盖|
|二元 Archimedean 依赖建模|4|5|5|4|5|5|4|4|36分，边界内通过并原子预留|

## 实际检索
已读取并执行本机 osc2026-guide 的 Project Research Guide：moon search --limit100 → MoonCakes docs → 相邻源码/测试；不是仅援引 skill 名。MoonBit0.10.18，Moon20261028。
Moon 搜索：phylogen,newick,parsimony,copula,copulas,dependence,survey,sampling,poststrat,clayton,gumbel,frank,kendall,rosenblatt,archimedean,联合分布,rank correlation,empirical copula,tail dependence,pseudo observations,联合概率,FrankCopula,clayton_cdf,h_function。dependence 饱和100；sampling83、survey8、copula3、copulas1，其他存在模糊匹配，不能从零结果推断唯一性。每条实际原始结果留在本地研究记录。
GitHub repository 搜索 per_page100：copula/phylogen/survey sampling/clayton/rosenblatt language:MoonBit 均0。code搜索 ClaytonCopula、empirical copula、rosenblatt、Horvitz 均0，copula仅1个无关 typst 文本；索引遗漏已知 MonteCarlo 源码，明确不视为全面覆盖。实际用树/contents读取已知相邻项目补足。检查 moonbitlang/awesome-moonbit 完整 README，未检出相关直替代。
永久登记表全标题/能力扫描，并读 MoonBench、MoonContract、MoonRecur、MoonShard、MoonIsotonic、MoonCoDA、MoonConformal、MoonChoice、MoonIRT、MoonQueueLab、MoonRunoff、MoonXport 对应指纹；这些分别是性能/HTTP/日历/分块/PAVA/组成数据/保形集合/离散选择/测验/排队/损失发展/XPORT，不以联合分布依赖为主循环。

## 相邻能力矩阵
|项目与证据|输入与核心输出|重合分类/决策|
|---|---|---|
|[jdjjttr/montecarlo@0.1.0](https://mooncakes.io/docs/jdjjttr/montecarlo)，[源码](https://github.com/jdjjttr/moonbit-montecarlo/blob/main/lib/multi_asset/copula.mbt)、multi_asset_test、mbti|Cholesky+RNG→相关正态/t向量与GBM价格路径；源接口 sample_gaussian_copula/sample_t_copula|同方向不同核心工作流；无 Archimedean CDF/条件密度/秩拟合。独立统计模型比向定价引擎追加整套分析层更清晰，不复制接口或源码。Apache-2.0，公开未归档，最近push2026-08-24|
|[IvanAXu/BioSeqs@0.1.14](https://mooncakes.io/docs/IvanAXu/BioSeqs)，[源](https://github.com/paipai-Studio/BioSeqs)|Newick树、距离、Fitch/Sankoff 简约；src/phylo.mbt、test/moonbit/parsimony_test.mbt 实读|生物候选直接覆盖，淘汰|
|[qjgqjg/moon-causal](https://mooncakes.io/docs/qjgqjg/moon-causal)，README/生成API|Horvitz-Thompson/Hajek、survey_design、weighted_quantile、rake_integer_cells|调查候选直接覆盖，淘汰；不是依赖分布模型|
|[mengshi-12/moon-rake@0.1.1](https://github.com/mengshi-12/moon-rake)，README、raking_test.mbt|分类边际→有界校准权重/残差|调查候选直接覆盖，淘汰|

MoonCakes 三个详情页成功200；statsmodels/copie文档成功200。初次猜测 tests/phylo_test.mbt 与 rake_test.mbt 返回404，随后用真实树恢复正确测试路径，未掩盖访问失败。对核心接口 FrankCopula/clayton_cdf/rosenblatt/empirical copula 做反搜。精确目标 CYang-dep/mooncopula@0.1.0 预检404。

结论：截至记录日期，在以上有界检索和实际源码范围内未发现成熟维护的直接替代。不是绝对唯一性或官方裁决。实现前使用同目录排他锁、重读和原子替换完成登记。发布前需再次检索并记录结果。

## 发布前反向复核（2026-11-11）

实现完成后重新执行 `moon search --limit 100`：copula 返回 3 项（MonteCarlo、MoonBallot、Symbit），Archimedean 与 Rosenblatt 无直接结果，checkerboard 返回 5 项，均为绘图/计算机视觉棋盘。原有 MonteCarlo 重合边界不变。GitHub `copula language:MoonBit` 仓库搜索及 `"bootstrap_raw" language:MoonBit` 代码搜索均返回 total_count=0、incomplete_results=false；这些仍受索引覆盖限制，不据此宣称绝对无重复。共享注册表再次检索得到本项目自己的 reserved 项与对应 external-checked 项，未出现新的同工作流冲突。

精确 GitHub 目标仓库及 MoonCakes 0.1.0 manifest 重新直查均返回 HTTP 404；两端隔离凭据服务端复验均为 CYang-dep。研究中的相邻包实际源码、测试与接口比较继续作为主要差异依据。当前实现兑现了定义的三条用户链路，新增诊断/重采样/校准均服务于拟合结果解释与联合事件评估，没有改成其它已预留项目的工作流。原始去秘密证据留本地验收目录；公开文档不含私人申报资料。
