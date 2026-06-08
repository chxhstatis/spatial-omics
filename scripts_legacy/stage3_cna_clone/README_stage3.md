# Stage 3：空间克隆 / 拷贝数分析（自包含，可在 HPC 直接运行）

胰腺癌（PDAC）空间 DNA-seq（DBiT 50×50，50 µm，hg38）。
**无需联网、无需下载**——GC/mappability/blacklist 参考已内置在 `ref/`。

---

## 0. 这个包能跑出什么

1. **（主，稳健）bulk 肿瘤 CNA 谱** + **空间 CNA-负荷/肿瘤分布图** —— `stage3_bulk_cna.py`
2. 空间克隆检测（PCA+KMeans，已做 GC + 深度混杂校正）—— `stage3_spatial_clones.py`
3. **逐 spot CNV**（每个 spot 一条完整全基因组拷贝数谱，2500 spots × ~2646 bins）—— `stage3_per_spot_cnv.py`

> ⚠️ **本数据的实测结论（务必先读）**：520 在 50 µm、当前覆盖与**低肿瘤纯度**下，**检测不到稳健的空间亚克隆**；但 **bulk 层面有干净、PDAC 一致的 CNA**（8q/MYC、12p/KRAS 增益；17p/TP53、18q/SMAD4、6q 缺失，幅度被间质稀释到 ~±0.2）。所以**主交付物是 bulk CNA 谱 + 空间负荷图**，克隆检测作为参考。

---

## 1. 依赖
```
python3 + pandas numpy scipy scikit-learn matplotlib
# 你的 conda 环境一般已有；否则：
conda install -c conda-forge pandas numpy scipy scikit-learn matplotlib
```

## 2. 输入（Stage-2 输出）
每个样本需要这两个文件（脚本自动识别 `<sample>/` 或 `<sample>_xy_map/` 目录）：
```
<stage2_dir>/<sample>_xy_map/top50x50_fragments_with_cb.tsv.gz   # chr,start,end,CB,mapq
<stage2_dir>/<sample>_xy_map/matched_spot_barcodes.tsv           # barcode_obs -> x_id,y_id
```
> ✅ 在 HPC 直接跑可避开传输截断问题（之前 525 的 `top50x50_fragments_with_cb.tsv.gz` 因 Stage-2 导出中断而截断，需先重跑 525 的 Stage-2 得到完整文件）。

## 3. 运行
```bash
# 一行跑全部（bulk + clones + per-spot CNV）：
bash run_stage3.sh <stage2_dir> stage3_out 520_520 525_525

# 或分别跑：
python3 stage3_bulk_cna.py      --stage2_dir <stage2_dir> --samples 520_520 525_525 --ref_dir ./ref --outroot stage3_out_bulk
python3 stage3_spatial_clones.py --stage2_dir <stage2_dir> --samples 520_520 525_525 --gc_file ./ref/bin_gc_1mb.tsv --outroot stage3_out_clones --bin_size 1000000 --spatial_win 7
python3 stage3_per_spot_cnv.py   --stage2_dir <stage2_dir> --samples 520_520 525_525 --ref_dir ./ref --outroot stage3_out_perspot --spatial_win 7 --cn_chroms chr8,chr17,chr18
```

## 4. 输出与解读
`stage3_out_bulk/<sample>/`：
- **`bulk_cna_profile.png/pdf/tsv`** ⭐ —— bulk 肿瘤 CNA 谱（红=增益、蓝=缺失，已标 KRAS/MYC/TP53/SMAD4/CDKN2A/GATA6）。看 8q/12p 是否红、17p/18q/6q 是否蓝。
- **`spatial_cna_burden_map.png`** —— 每个 spot 与 bulk 肿瘤 CNA 的相似度（肿瘤含量代理）。黄/红=肿瘤富集区。本数据偏弱/弥散。
- `bulk_chrom_summary.tsv` —— 各染色体平均 log2 + corr(CNA,GC)（应 ≈0，证明 GC 已校正干净）。
- `spot_cna_burden.tsv` —— 每 spot 的 (x,y, total_frags, cna_burden)。

`stage3_out_clones/<sample>/`：
- `spatial_clone_map.png`、`clone_cna_heatmap.png`、`clone_cna_profiles.png`、`spot_clone_labels.tsv`、`summary.tsv`（含 `depth_confound_r_PC1_vs_logcov`，应 ≈0）。
- **`region_significance_maps.png`** ⭐（新）—— 论文 Fig 2c 式：默认 chr8/17/18 的**空间显著性图**（signed −log10 p，红=覆盖显著高/缺失低）；`region_zscores.tsv` = 每 spot 每区域 z 值。

`stage3_out_perspot/<sample>/`（**逐 spot CNV**）：
- **`per_spot_cnv_matrix.tsv.gz`** ⭐ —— 核心交付：每行一个 spot（`x_id,y_id,total_frags`），每列一个 1Mb bin，值=拷贝数（基线 2）。2500 spots × ~2646 bins。
- **`spatial_cn_<chr>.png`** —— 把某条染色体的逐 spot 拷贝数画回 50×50 网格（红=增益、蓝=缺失，中心=2）。默认 chr8/17/18。
- `bin_info.tsv` —— 矩阵列顺序对应的 bin → chr,start,end。
- `per_spot_cnv_burden.tsv` —— 每 spot 的 `|CN−2|>cna_thresh` 基因组占比（CNV 负荷）。
- `example_spot_profiles.png` —— 几个示例 spot 的全基因组拷贝数曲线。
- **原理**：单个 50µm spot 在 1Mb bin 上 ~1 fragment，太稀疏；故每个 spot 的 CNV 用**空间邻域平滑**（`--spatial_win`，默认 7×7 ≈ 论文 kNN≈49）从「自身 + 邻近 spot」聚合算出。**有效分辨率 = 邻域大小**：调小更接近单点但更噪，调大更平滑但更粗。

## 4a. normal-anchor 接口（`--normal_spots`）与一个重要教训
`stage3_bulk_cna.py` 和 `stage3_per_spot_cnv.py` 都支持 `--normal_spots <x_id,y_id 的 tsv>`：用这些"正常/参考" spot 当 pseudo-normal 基线（bulk 是把肿瘤聚合谱**除以**正常聚合谱，类比论文 ÷bulk-WGS；逐 spot 是把参考改成只在这些 spot 上取中位）。这是给**病理确认的正常区 / Stage-6 对照样本 / Stage-4 RNA 定的间质**预留的接口。

`pick_normal_spots.py` 可从逐 spot CNV 负荷**自动挑最平坦的一簇** spot 生成 `normal_spots.tsv`（数据驱动，不需标注）。
```bash
python3 pick_normal_spots.py --burden <perspot>/<sample>/per_spot_cnv_burden.tsv --out <out>/<sample> --quantile 0.30
```
> ⚠️ **实测教训（520）**：本切片**低纯度 + 肿瘤铺满**，不存在可数据驱动识别的内部正常 spot。用自动挑的"最平坦"簇做锚定，会把真信号一并除掉——log2 SD 从 0.205 塌到 0.035（17%），TP53/SMAD4/MYC 等已知 CNA 全归零。`stage3_bulk_cna.py` 内置 **signal-collapse 告警**会自动检出这种情形并提示"保留未锚定结果"。
> **结论**：对均匀低纯度切片，**正确的正常参考必须来自外部**（Stage-6 对照 / 真实正常区 / RNA 间质），不能用内部 spot 自锚。520 的可发表 bulk CNA 仍以**基因组中位归一的未锚定版**为准。

## 4e. 逐样本空间异质性（`stage3_spatial_heterogeneity.py`）
量化"是否存在真实空间 CNA 异质",用 Moran's I 但**先扣两个混杂**：
```bash
# 先生成 win=1（不平滑）矩阵作诚实基线，再跑异质性（同时传平滑与 win1 两个 root）：
python3 stage3_per_spot_cnv.py --stage2_dir <s2> --samples ... --ref_dir ./ref --outroot <perspot_win1> --spatial_win 1
python3 stage3_spatial_heterogeneity.py --perspot_root <perspot> --raw_perspot_root <perspot_win1> --samples ... --out <out>
```
输出每样本 `cnv_burden_map.png`/`spatial_subregion_map.png`(探索性)/`subregion_cnv_profiles.png` + 队列 `cohort_morans_I.png`。
> ⚠️ **两个必须扣掉的混杂**：① **平滑**——逐 spot CNV 用 N×N 邻域平滑生成,相邻 spot 共享数据,平滑矩阵的 Moran's I（~0.7–0.9）是构造产物,**不是生物学**；必须用 win=1 矩阵测真实自相关。② **覆盖**——组织密度本身空间自相关;真实 CNA 异质的 `cna_beyond_coverage` 须 = win1 负荷 Moran > 覆盖 Moran。**实测 4 样本全 False**（详见 `stage3_4samples/HETEROGENEITY_RESULTS.md`）：低纯度+50µm 下无稳健空间 CNA 异质,聚类图仅探索性。

## 4d. 多样本队列对比（`compare_samples_cna.py`）
跑完多个样本的 bulk CNA 后，做跨样本队列分析（找复发 vs 私有 CNA、样本相似度）：
```bash
python3 compare_samples_cna.py --bulk_root <stage3_out_bulk> --samples 520_520 525_525 A_A B_B --out <out>
```
输出：`cohort_cna_profiles.png`（各样本谱对齐）、`recurrence_track.png`（每 bin 多少样本增/缺）、`sample_correlation.png` + `.tsv`、`recurrence.tsv`、`recurrent_genes.tsv`、`inferred_sex.tsv`。
> ⚠️ 两个易踩的坑（已内置处理/告警）：① **样本相似度默认排除 chrX**——否则会按**性别**聚类而非肿瘤生物学（实测 520/B_B 女、525/A_A 男，含 chrX 时假性配对 0.94/0.82，排除后均匀 0.72–0.80）；② **着丝粒附近的复发性强缺失（chr6/7/10 着丝粒）多为技术伪影**，不是复发驱动事件，解读复发轨时要剔除着丝粒/低 mappability 区。

## 4b. 已融入的论文方法（slide-DNA-seq, Nature 2022；详见根目录 `slide_dna_seq_论文流程精髓.md`）
- **双重平滑**（`double_smooth_pcs`，论文 pc_scores_smo_both）：PCA 后在 **PC 空间 + xy 空间**各 kNN 平滑再联合，然后才 KMeans——对抗稀疏、不抹克隆边界。
- **置换 z-score 空间显著性**（`permutation_zscore_region`，论文 Fig 2c）：随机重分配 fragment→spot 建经验零分布 → z → 减 normal 均值 → 空间 p 值图。参数 `--zscore_chroms chr8,chr17,chr18 --n_perm 100`。
- **bulk-WGS 除法接口**：`fig1f_aggregate_cna.py --control <bulk_wgs_1Mb.tsv>`（待对照/正常 bulk 数据到位即可逐位复刻论文 Fig 1f 的 ÷control 步）。

## 4c. 论文 Fig 1f/1j 复刻与对比工具
`fig1f_aggregate_cna.py`：全 fragment 1Mb 聚合 → 拷贝数（基线2）沿染色体分布，输出 RAW vs GC+mappability 校正两版（论文 Fig 1f/1j 样式）。
```bash
python3 fig1f_aggregate_cna.py --frag <stage2>/<sample>/top50x50_fragments_with_cb.tsv.gz \
  --ref_dir ./ref --sample <id> --out <out>/<id>_aggregate_1Mb   # 有 bulk-WGS 再加 --control
```
> 我们的图 ≈ 论文 **Fig S6b**（做了 GC/map 校正、未除 bulk-WGS）；走到干净的 Fig 1f 就差"÷ bulk-WGS"。

## 5. 方法学（流程与校正逻辑）
1. fragment → (x,y)：用 `matched_spot_barcodes.tsv` 把每条 fragment 落到 50×50 网格。
2. 1Mb 分箱 → spot×bin 矩阵。
3. **bin 质控**：去 blacklist（`ref/hg38-blacklist.v2.bed.gz`）、低 mappability（`ref/k100.umap.bed.gz`，k100<0.5 弃）、极端 GC、低覆盖 bin。
4. **GC 校正**（必需）：按 GC 分位做 count 去偏——否则 GC-rich 小染色体会冒充"增益/缺失"。
5. **mappability 校正**：按 mappability 分位去偏——清掉小染色体残留偏差。
6. **bulk CNA**：全 spot 加总 → 校正 → 除以基因组中位（设 modal=二倍体）→ log2。
7. **空间负荷**：每 spot profile（GC/map 校正 + 空间平滑 + 每 bin 全 spot 中位为参考）与 bulk 肿瘤谱求相关 → 画回网格。
8. 克隆检测额外做：**回归掉 log 覆盖**消除深度混杂（这是 v1→v3 的关键教训：不校正会按测序深度而非 CNA 聚类）。

> 关键教训：稀疏空间 DNA 的 CNV 必须先做 **GC + mappability + 深度混杂** 校正，否则覆盖/GC 伪影会冒充克隆结构（v1 全蓝伪影即因此）。

## 6. 结果回传
输出都是小文件（PNG/PDF/TSV，每样本几 MB），跑完直接把 `stage3_out_bulk/` 和 `stage3_out_clones/` 打包传回即可。
