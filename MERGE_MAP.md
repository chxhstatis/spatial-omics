# 代码合并地图（三套 → 一套，2026-06-08）

把之前并存的代码合并成**唯一规范仓库** `spatial-omics/`（import 名 `spatial_omics`）。

## 合并前的实际情况

| 原代码 | 性质 | 去向 |
|---|---|---|
| `spatial_omics_pipeline/`（13 个 stage 脚本，真实跑过 4 样本 + 100G） | 成熟的 stage1–6 脚本 + 内置参考 | stage1–2 → `workflow/`；stage3–6 → `scripts_legacy/`；算法逐步上提为 API |
| `spatial-omics`（handoff 包，AnnData 原生 pp/tl/pl，1129 行） | 产品形态的下游分析包（合成数据 ARI 0.96） | **作为本仓库的底座** `src/spatial_omics/` |
| `spacecna`（计划附录提到的"P0 内核"） | **盘上不存在**——当时那个内核即后来的 handoff 包，仅名字变过 | 无需处理（同一份东西） |

> 即"三套"实为"两套真代码 + 一个别名"。

## 合并后的结构

```
spatial-omics/
├── src/spatial_omics/          # ★下游分析包（产品中心，AnnData 原生）
│   ├── io.py                   #   from_pipeline / from_matrix / read_h5ad
│   ├── datasets/               #   simulate()（合成 ground-truth，无需病人数据）
│   ├── pp/                     #   bin_qc / correct_bias / normalize / add_reference
│   ├── tl/                     #   dual_smooth / call_clones / copy_number / permutation_significance
│   │   └── rigor.py            #   ★严谨层（差异化护城河，见下）
│   ├── pl/                     #   spatial_clones / copy_number / significance 等画图
│   ├── _ref/                   #   内置参考（bin_gc / k100.umap / blacklist，合并成一份）
│   └── _constants.py / _bins.py
├── workflow/                   # 上游（你芯片专属的 "Cell Ranger 等价物"）
│   ├── stage1_preprocess_qc/   #   FASTQ → 注入条码 → 比对 → 去重 → fragment
│   └── stage2_spatial_map/     #   16bp 条码 → 50×50 真实坐标（X/Y 192 校正表）
├── scripts_legacy/             # stage3–6 原脚本 = 参考实现（API 未覆盖处仍可直接用，不丢）
│   ├── stage3_cna_clone/       #   bulk_cna / per_spot_cnv / spatial_clones / heterogeneity /
│   │                           #   pick_normal_spots / compare_samples_cna / compute_bin_gc / fig1f
│   ├── stage4_rna_integration/ #   PLAN（待 RNA）
│   ├── stage5_wsi_registration/#   register / purity_he（H&E 颜色解卷积）
│   └── stage6_control_compare/ #   PLAN（待对照）
├── docs/ (mkdocs-material) · tests/ · examples/ · pyproject.toml · CITATION.cff
```

## ★严谨层 `tl/rigor.py`（这是与 inferCNV/CopyKAT/官方 MATLAB 的核心差异）

把"防伪克隆"的混杂控制从 stage3 脚本上提为一等 API。低纯度+稀疏数据上，naive 聚类会把覆盖密度/平滑核/通道条纹错当克隆；这一层在报告克隆前先证伪：

| API | 做什么 | 来源 |
|---|---|---|
| `tl.morans_i(values, x_id, y_id)` | 8 邻接 Moran's I + 真随机置换零分布 | stage3_spatial_heterogeneity |
| `tl.spatial_heterogeneity(adata)` | CNV 负荷 Moran vs **覆盖基线** → `cna_exceeds_coverage` 判据 | 同上 |
| `tl.clone_diagnostics(adata)` | 克隆 CNA **可区分度** + **CH 边界**诊断 → `clones_likely_artifact` | call_clones + 报告分析 |
| `tl.detect_channel_stripes(adata, key)` | 检测网格对齐的行/列条带 = 微流控通道伪影 | 100G 区域显著性分析 |

验证：合成数据（有真克隆）→ 判"真异质/非伪影"；真实低纯度 PDAC → 判阴性/伪影（与报告 1/2/3 一致）。**同一套函数，干净数据上确认克隆、脏数据上拦截伪克隆。**

## 待办（API 化路线，不阻塞合并）

- ✅ `pp.normal_anchor`（+ 信号塌缩自检）和 `pp.pick_normal_spots` —— 已上提为 API。
- ✅ `tl.he_purity`（H&E 颜色解卷积）—— 已上提（`[he]` 可选依赖 cv2/skimage，惰性导入）。
- `tl.cohort_compare`（跨样本 + chrX 性别排除，compare_samples_cna）→ 上提。
- `tl.variance_decomposition`（stage4，待同片 RNA）= multi-modal 灵魂。
- `workflow/` 容器化（nextflow/snakemake + Docker）= 上游一键化。

## 旧仓库处置

- `spatial_omics_pipeline/`：保留（git 历史 + scripts_legacy 的来源），标记为**被本仓库取代**；SSD 镜像仍是历史提交主线。
- `spatial-omics_handoff.tar.gz`：已归档（其内容即本仓库底座）。
