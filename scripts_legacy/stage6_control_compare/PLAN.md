# Stage 6：对照组比较 / 正常基线（待数据）

**触发条件**：`samples.tsv` 出现 `group=control`（或 `normal_adjacent`）样本。

## 目标
用真实对照/正常样本建立**正常 CNA 基线（Panel of Normals, PoN）**，替代当前"数据内猜 pseudo-normal"，让肿瘤 CNA 对比更可靠，并做肿瘤 vs 对照的差异 CNA。

## 输入
- 对照样本：与肿瘤同流程跑到 Stage-3（同芯片设计、同参考）。
- 肿瘤样本 Stage-3 输出。

## 计划步骤
1. **对照过流程**：control 样本走 stage1→2→3，得到其 bin 级覆盖矩阵。
2. **建 PoN**：多个正常/对照 spot 聚合 → 每 bin 的正常期望（中位 + 离散度），作为二倍体基线。
3. **肿瘤归一化到 PoN**：肿瘤 bin 覆盖 / PoN 期望 → log2 CNA（比数据内 pseudo-normal 更稳，尤其低纯度）。
4. **批次/技术校正**：若对照与肿瘤不同批，先按 GC/mappability + 批次因子对齐（避免把批次当 CNA）。
5. **差异 CNA**：肿瘤 vs 对照逐 bin 检验（置换 + Wilcoxon），定位肿瘤特异扩增/缺失；标 PDAC 基因。

## 与现有的接口
- 对照的 bin 矩阵 → 作为 `stage3_*` 的外部参考（替换/补充当前 `process()` 里的 per-bin 全 spot 中位参考）。
- 复用 stage3 的 GC/mappability/blacklist 校正与绘图。

## 产出
PoN 基线、对照锚定的肿瘤 CNA 谱（更干净）、肿瘤特异差异 CNA 表与图。
