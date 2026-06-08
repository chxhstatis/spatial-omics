# Stage 4：同片 RNA 整合（待数据）

**触发条件**：`samples.tsv` 某样本 `rna=yes`（同片或相邻片空间 RNA 已返回）。

## 目标
把空间转录组与 Stage-3 的空间 DNA/CNA 联合，回答"基因表达由克隆遗传（内因）还是微环境（外因）驱动"，并交叉验证肿瘤区。

## 输入
- RNA：count 矩阵 + 空间坐标（Slide-seqV2/Stereo-seq/DBiT-RNA；或 fastq 经 Cell Ranger/SAW/Stickels 流程）。
- DNA：Stage-3 输出（`bulk_cna_profile`、`spot_cna_burden`、克隆标签）。

## 计划步骤
1. **RNA 预处理/聚类**：sctransform 归一化 → PCA/UMAP → Seurat 聚类 → 细胞类型注释（标志基因）。
2. **细胞类型反卷积**：RCTD（spacexr）按参考做 per-spot 细胞类型组成；得到肿瘤/间质/免疫密度。
3. **RNA 端独立推断 CNV**：inferCNV / SlideCNA（表达感知空间分箱）→ 与 DNA 的 bulk CNA 谱比对（Pearson），交叉验证。
4. **跨模态共配准**：DNA 与 RNA 阵列用 MATLAB `imregister`/手动 landmark 配准；把 DNA 的肿瘤/克隆标签传到 RNA spot（最近邻）。
5. **方差分解**（核心）：对每个基因，把表达方差拆为「克隆/CNA 身份 + 肿瘤密度（微环境代理）+ 未解释」（逐步回归），找内因 vs 外因驱动基因集 → 基因集富集（hypeR/MSigDB）。

## 与现有的接口
- 用 `stage3_*/spot_cna_burden.tsv`（肿瘤含量代理）作为 RNA 端肿瘤定位的 DNA 锚。
- 用 `stage3_*/bulk_cna_profile.tsv` 与 RNA 推断的 CNV 比对验证。

## 产出
RNA 细胞类型空间图、DNA↔RNA CNV 一致性、内因/外因驱动基因表、跨模态肿瘤区一致性图。
