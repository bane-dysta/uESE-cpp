# uESE 的公式、数值实现与适用范围

本程序实现 Vyboishchikov 与 Voityuk 的 **2021 uESE 模型**[1]。模型把溶剂化自由能写成

\[
\Delta G_{\mathrm{solv}}=E_{\mathrm{elst}}+\Delta G_{\mathrm{corr}}.
\]

前者由固定气相原子电荷的 COSMO 方程计算；后者是拟合的校正，包含空腔、色散、短程相互作用及未显式计算的极化效应。因此，不能把全部校正项严格理解为纯非静电能，也不能用 ESE-PM7 或 ESE-EE 的参数替换它。[1]

## 1. 输入电荷与表面电荷是不同的量

输入包含气相分子几何 \(\mathbf R_A\) 和 **CM5 原子电荷** \(Q_A\)。原文以 B3LYP/def2-TZVP 气相结果拟合参数。本程序读取这些电荷，不执行量子化学计算，不把其他布居分析电荷自动转换成 CM5，也不重新归一化电荷总和。[1]

每个暴露的表面片有位置 \(\mathbf r_i\)、面积 \(S_i\) 及待求的导体极限诱导电荷 \(q_i\)。按所属原子汇总：

\[
S_A=\sum_{i\in A}S_i,\qquad
q_A=\sum_{i\in A}q_i,\qquad
\sigma_A=q_A/S_A.
\]

大写 \(Q_A\) 是输入的原子电荷；小写 \(q_A\) 是该原子表面上的诱导电荷。校正项使用后者。完全埋藏的原子不进入 \(\sigma_A\) 的极值计算，以避免 \(0/0\)。

这里的 \(q_i\) 已经是一个表面片上的总电荷，不是电荷密度。因此，汇总 \(q_A\) 或计算静电能时都不能再乘一次 \(S_i\)。面积只在对角自相互作用、面积校正和密度定义等规定位置出现。电荷总和应作为诊断量输出；有限网格的诱导电荷总和与 \(-\sum_AQ_A\) 不完全相等时，本程序不强行缩放以抹去差异。

## 2. 静电项：原文式 (4)–(7)

先计算气相电荷产生的表面势，再解一个线性方程组：

\[
V_i=\sum_A\frac{Q_A}{|\mathbf R_A-\mathbf r_i|},\qquad
\mathbf A\mathbf q=-\mathbf V,
\]

\[
A_{ij}=\frac{1}{|\mathbf r_i-\mathbf r_j|}\quad(i\ne j),\qquad
A_{ii}=1.0694\sqrt{\frac{4\pi}{S_i}}.
\]

\[
E_{\mathrm{elst}}=\frac{K}{2}f(\epsilon)\sum_iq_iV_i,
\qquad
f(\epsilon)=
\begin{cases}
(\epsilon-1)/\epsilon,&\text{A、B、C 类},\\
(\epsilon-1)/(\epsilon+1/2),&\text{D 类}.
\end{cases}
\]

坐标以 Å、电荷以基本电荷 \(e\) 的数值表示时，程序使用单位换算常数 \(K=332.063713299\)，得到 kcal/mol。原文式 (4) 中的物理单位换算在这里显式写出；它不是拟合参数。介电因子只乘静电能，**不再次乘入用于校正项的 \(q_i\) 或 \(q_A\)**。[1]

“非迭代”表示不反复更新溶质电子密度。代码优先用 Cholesky 分解；离散矩阵不正定时改用带选主元的 LU 并报告警告，最后检查原方程残差。不正定矩阵的求解成功不等于该几何的离散误差可以忽略。

## 3. 四类溶剂的校正：原文式 (8)–(12)

下式的 \(V\) 是空腔体积，与表面势 \(V_i\) 不同；\(S=\sum_AS_A\)。

| 类别 | 溶剂类型 | \(\Delta G_{\mathrm{corr}}\) | 原文 |
| --- | --- | --- | --- |
| A | 水 | \(\sum_A\kappa_AS_A+\sum_Ap_Aq_A^2+\mathrm{SRC}\) | 式 (8) |
| B | 非水极性质子溶剂 | \(\sum_A\kappa_AS_A+\sum_Ap_Aq_A^2+\mathrm{SRC}+\zeta V\) | 式 (10) |
| C | 极性非质子溶剂 | \(\sum_A\kappa_AS_A+\xi_{\mathrm{solv}}S+\sum_Ag_Aq_A+\sum_Ap_Aq_A^2+\zeta V\) | 式 (11) |
| D | 非极性溶剂 | \(\sum_A\kappa_AS_A+\xi_{\mathrm{solv}}S+\sum_Ap_A\lvert q_A\rvert+\zeta V\) | 式 (12) |

原文把介电常数低于 9 的溶剂归入 D 类，即使它可能形成氢键。调用内置溶剂时使用作者程序的分类，不根据溶剂名称自行推断。体积系数依次为 A：0，B：\(-0.070\)，C：\(-0.236\)，D：\(-0.166\)。**C、D 的负号不能省略。**[1]

为避免参数名称与原子指标、溶剂类别混淆，将原文式 (9) 的 \(A,B,C\) 参数记为 \(a_s,b_s,c_s\)：

\[
\mathrm{SRC}=c_s\operatorname{erf}
\left(a_s|\sigma_{\min}|^3+b_s|\sigma_{\max}|^3\right),
\quad
\sigma_{\min}=\min_{A:S_A>0}\sigma_A,\quad
\sigma_{\max}=\max_{A:S_A>0}\sigma_A.
\]

| 溶剂类别 | \(a_s\) | \(b_s\) | \(c_s\)，kcal/mol |
| --- | ---: | ---: | ---: |
| A | 68500 | 87700 | -8.87 |
| B | 186000 | 78500 | -2.36 |

由于两个密度项都非负，误差函数的自变量也非负；结合 \(c_s<0\)，此处 SRC 在 \(c_s\) 与 0 之间。它在较强局部电荷密度下饱和，而非随密度无限降低。这一界限也是检查短程项符号或公式括号是否写错的简单办法。

参数 \(\kappa_A\) 和 \(\xi_{\mathrm{solv}}\) 的单位为 kcal mol\(^{-1}\) Å\(^{-2}\)，\(\zeta\) 为 kcal mol\(^{-1}\) Å\(^{-3}\)；以电子电荷的数值代入时，\(p_A,g_A\) 按表中 kcal/mol 数值使用。原 Table 1 将 \(a_s,b_s\) 标为 Å\(^{-6}\)，但这与 \(\sigma^3\) 的量纲存在疑点。实现严格保留表中数值，使用 \(q_A/e\) 与 \(S_A/\mathrm{Å}^2\) 的数值计算密度，不据单位标注擅自改变系数。[1]

## 4. 空腔离散及 NDIV 的含义

原文使用不缩放的 Bondi 半径、GEPOL93 算法，且不添加非原子中心的平滑球。[1–3] 本程序独立构造球面网格：32 个顶点形成 60 个初始三角面；细分后剔除其他原子球内部的采样点，再把每个初始面剩余的子点合并成一个表面电荷点，并投影回原子球面。这里使用等权子点面积，不是重新定义为精确球面三角形积分。

原文文字写“32 tesserae per atom”，本实现与作者程序核对所采用的几何则是 **32 顶点、60 初始面**。两者在术语上有差异，不能将“32”直接理解成 32 个最终表面电荷。默认 **NDIV=3**：每个初始面有 \(4^{3-1}=16\) 个遮蔽采样子点，孤立原子仍合并为 60 个电荷点。

增加 NDIV 主要细化遮蔽判断、面积和合并位置，**并非完整的 COSMO 电荷网格收敛**，也未重新拟合参数。有限离散还会产生方向误差；体积用 \(V\approx\frac13\sum_iS_i\mathbf r_i\cdot\mathbf n_i\) 计算，残余的闭合误差可能导致坐标平移敏感性。`--center` 可显式居中，但会改变与未居中参考计算的对照，应记录此选项。

## 5. 参数出处、适用范围与标准态

正文式 (4)–(12) 及 Table 1 已根据用户上传的完整 2021 年发表版本核对（参数见 [table1_parameters.json](../data/table1_parameters.json)）。**Supporting Information 第 2 页 Table S1 已逐项核对 77 个 C/D 类溶剂的 \(\xi_{\mathrm{solv}}\) 值，全部一致**，原始核对表见 [table_s1_xi.tsv](../data/table_s1_xi.tsv)。92 种溶剂的 \(\epsilon\) 仍来自**作者 uESE v1.2 的命名运行时参数**，因为 Table S1 不提供介电常数。见 [solvents.json](../data/solvents.json)、[solvent_provenance.json](../data/solvent_provenance.json) 及 [论文核对报告](PAPER_CROSSCHECK.zh-CN.md)。[1,4]

Table 1 的 A/B/D 类支持 H、C、N、O、F、P、S、Cl、Br、I；C 类只支持其中除 P、I 外的八种元素。B、Au 均无已发表校正参数。缺项与数值为零不同：本程序拒绝为缺项输出完整 uESE 自由能；`--electrostatic-only` 只给静电贡献，不能当作 \(\Delta G_{\mathrm{solv}}\)。

程序输出上述模型公式之和，**不自动添加 1 atm 与 1 mol/L 标准态之间的转换项**。将结果与气相热化学校正组合时，应先核对目标过程的标准态，避免重复或漏加转换；静电项和校正项的复现测试本身不能证明标准态已处理完毕。

## 6. 公式与源码对应

| 内容 | 实现位置 |
| --- | --- |
| 式 (4)–(7)、线性方程与残差 | [src/model.cpp](../src/model.cpp) |
| 式 (8)–(12)、各项能量分解 | [src/model.cpp](../src/model.cpp) |
| Table 1、半径及溶剂参数入口 | [src/parameters.cpp](../src/parameters.cpp) |
| 空腔、面积、体积和归属 | [src/surface.cpp](../src/surface.cpp) |

## 参考文献

1. S. F. Vyboishchikov, A. A. Voityuk. *Fast non-iterative calculation of solvation energies for water and non-aqueous solvents*. J. Comput. Chem. **42** (2021), 1184–1194. [DOI: 10.1002/jcc.26531](https://doi.org/10.1002/jcc.26531)。本文件主要对应印刷页 1185–1187。
2. J. L. Pascual-Ahuir, E. Silla, I. Tuñón. *GEPOL: An improved description of molecular surfaces. III. A new algorithm for the computation of a solvent-excluding surface*. J. Comput. Chem. **15** (1994), 1127–1138. [DOI: 10.1002/jcc.540151009](https://doi.org/10.1002/jcc.540151009)。
3. A. Bondi. *van der Waals Volumes and Radii*. J. Phys. Chem. **68** (1964), 441–451. [DOI: 10.1021/j100785a001](https://doi.org/10.1021/j100785a001)。
4. S. F. Vyboishchikov. [uESE 作者程序与说明](https://github.com/vyboishchikov/ESE)。本项目为独立 C++ 实现，不附带论文全文或作者可执行程序。
