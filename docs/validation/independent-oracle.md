# 独立数值预期

测试向量由 `scripts/generate_oracle.py` 生成，固定 Python 依赖见 `scripts/oracle-requirements.txt`。生产库不依赖 Python。CI 重新生成并与提交向量比较；生成文件不计入生产源码规模。

- statsmodels 0.14.5：Clayton/Gumbel CDF 与密度，固定种子 Clayton 样本。
- Python Decimal 80 位十进制运算：Frank CDF、密度、条件概率，正负分支以及 1e-10 附近边界；用于避开普通双精度参考自身的消减误差。
- SciPy 1.16.2：Kendall tau-b、Spearman、平均秩（本库除以 n+1）、KS、Cramer-von Mises；有界最优化 statsmodels 对数似然。
- 独立手算：独立模型事件分解、Checkerboard 边缘、Wilson 区间、type-7 分位数、经验计数和所有旋转的微分/积分不变量。

来源：
- https://www.statsmodels.org/stable/distributions.html#copula
- https://docs.scipy.org/doc/scipy/reference/stats.html
- https://docs.python.org/3/library/decimal.html

数学公式原创实现；没有移植上述软件源码。固定合成数据由本项目脚本创建，无第三方真实个人数据。参考通过不代表任意极端输入均有相同误差保证；公开 API 在数值失败时返回明确异常，未返回伪造的收敛结果。向量容差按具体量设置，CDF 绝对误差 2e-9、log-density/Frank conditional 2e-8；不将这些测试范围外推为全域定理。
