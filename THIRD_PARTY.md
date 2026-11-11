# 来源与许可边界

生产源码为本项目原创实现，MIT。没有移植其它项目源码。数学定义本身与实现表达区分对待，不将参考项目的完整能力归于本项目。

| 来源 | 用途 | 许可与处理 |
|---|---|---|
| moonbitlang/core | 唯一运行时依赖：数组、数学、Debug 等 | Apache-2.0，随官方工具链提供，不将其源码作为本项目生产行数 |
| statsmodels 0.14.5 | 开发期 Copula 方程与独立测试生成 | BSD-3-Clause，pip 安装使用，不复制/分发实现；公开文档 https://www.statsmodels.org/stable/distributions.html#copula |
| SciPy 1.16.2 / NumPy 2.3.3 | 独立统计预期、优化与固定种子合成数据 | BSD-3-Clause，不随生产包分发，开发环境保留其原许可证 |
| Python Decimal | 80 位精度独立算术参考 | PSF license，不分发标准库实现 |
| pandas 3.1.1 及 statsmodels 传递依赖 | 开发参考运行环境 | 依各安装包许可证，不嵌入项目源码 |
| jdjjttr/montecarlo | 查重边界分析 | Apache-2.0，阅读其 Gaussian/t 源码与接口，未复制；不宣称与其兼容 |

固定种子样本由生成器自行创建，没有使用第三方私人观测。生成的测试数据文件有明确标记，不作为生产实现统计。旋转、条件 CDF、秩变换与 Wilson/type-7 定义使用通行数学公式，以手算不变量和不同软件实现交叉验证。
