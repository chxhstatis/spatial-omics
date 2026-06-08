# Stage 5：WSI / H&E 配准 · 肿瘤标注 · 锚定 CNA（待数据）

**触发条件**：`samples.tsv` 某样本 `wsi=yes`（全片扫描 WSI / H&E 图像已到）。

> **为什么这步最关键**：本项目是低纯度 PDAC（间质多）。Stage-3 的 pseudo-normal 参考是"靠跨 spot 中位猜"的——只要肿瘤铺得开，参考就被肿瘤污染、CNA 幅度被压扁。**一旦有病理真实的肿瘤区/间质区标注，就能用"间质 spot"当真正的正常参考**，重算更干净、幅度更真实的 bulk CNA 和逐 spot CNV。这是把"温和、弥散"的当前结果提升为"有空间结构、可发表"的最大杠杆。

---

## 输入
- **病理图像**：WSI/H&E（`.svs/.ndpi/.tiff/.png`）。**三种采集情形分支处理**（见步骤 2）。
- **Stage-2**：`matched_spot_barcodes.tsv`（spot→x_id,y_id，50×50 网格坐标）。
- **Stage-3**：`per_spot_cnv_matrix.tsv.gz`（逐 spot 拷贝数，本流程新产物）、`spot_cna_burden.tsv`、`bulk_cna_profile.tsv`。
- **芯片几何**：DBiT 50×50，50 µm pitch → 捕获区 **2.5 mm × 2.5 mm**；fiducial / 边界 / 四角对位信息（若有）。

## 已确定的情形（用户确认 2026-06-05）
- **A. 图像与切片关系 = 同片、芯片在位明场图**（能看到 50×50 捕获区/通道）→ **配准走最简单稳健路径：检测捕获区四角/通道直接求仿射**。这是本 stage 最确定的一步。
- **B. 无病理医生圈注** → 走**自动 + 数据驱动**两条腿，且对自动肿瘤判读**保持克制**：
  - **现实评估**：芯片在位明场图的组织形态常被通道/芯片叠加干扰，**纯图像自动判肿瘤/间质不可靠**，不能作为唯一依据。
  - **B1 图像侧（弱）**：先做稳健能做的——组织 vs 背景掩膜、核密度/着色强度作为"细胞密度"代理（HoVer-Net 可选，但通道伪影需谨慎），**只当辅助证据**。
  - **B2 数据驱动锚定（已实现 + 已实测，结论：本切片不适用）**：`pick_normal_spots.py` 能自动挑最平坦簇生成 `normal_spots.tsv`，`stage3` 用 `--normal_spots` 接收。**但 520 实测证明此路在本切片不通**：低纯度+肿瘤铺满 → 无内部正常 spot → 自锚把真信号除掉（log2 SD 0.205→0.035，已知 CNA 归零；`stage3_bulk_cna.py` 已内置 signal-collapse 告警自动检出）。**因此本切片的正常参考必须来自外部（B3 或 Stage-6 对照）**，不能内部自锚。
  - **B3 金标准（待 Stage-4 / Stage-6）**：① 同片 RNA 定间质/免疫 spot（Stage-4）；② **对照样本 Panel-of-Normals（Stage-6）**——对均匀低纯度切片，这才是最可靠的正常锚。两者都通过同一个 `--normal_spots` 接口接入。

---

## 计划步骤（含具体工具、产物文件）

### 1. 读图 + 缩放
- `openslide` / `tifffile` 读 WSI；取一个下采样层（thumbnail），分辨率换算到与 50×50 网格可比（已知 pitch=50 µm，每 spot 对应的像素尺寸 = 50µm ÷ 图像 µm/px）。
- 产物：`wsi_thumb.png` + `wsi_meta.json`（mpp、下采样倍数、层尺寸）。

### 2. 配准（spot 网格 ↔ 图像）→ 仿射变换【情形已定：芯片在位，直接做】
- 目标：求 **`spot_to_image` 2×3 仿射矩阵**，把任意 (x_id,y_id) 映射到图像像素坐标。
- 做法：检测捕获区 —— ① 自动：`cv2` 边缘/Hough 检测 50×50 通道边界 → 拟合外接四边形取四角；② 兜底：在 thumbnail 上手动点选捕获区四角（最稳，一次即可）。四角 → `cv2.getAffineTransform`/`estimateAffinePartial2D`。
- 产物：`spot_to_image_affine.npy`、`spot_pixel_coords.tsv`（x_id,y_id → px,py）、**`chip_grid_overlay.png`**（原图上画出 50×50 网格做目视 QC）、每 spot 切出 50µm tile。
- **脚本**：`stage5_register.py`

### 3. 组织学代理 → 每个 spot 一个"密度/掩膜"标签【自动，弱证据】
- **不做"自信的肿瘤判读"**（无 H&E/无圈注）。只做稳健的：
  - 组织 vs 背景掩膜（Otsu/染色强度）→ 标记每 spot 是否落在组织上、是否折叠/坏死区。
  - 细胞密度代理：核密度（HoVer-Net 可选）或着色强度 → `spot_density.tsv`，仅作辅助。
- 产物：**`spot_tissue_label.tsv`**（x_id,y_id → on_tissue / fold / background + density）。**不输出肿瘤/间质硬标签**，避免过度解读。
- **脚本**：`stage5_annotate_to_spots.py`

### 4. 锚定 pseudo-normal → 回灌 Stage-3 重算（**核心收益**）
- **主路 B2（数据驱动，不依赖标注，立即可做）**：对逐 spot CNV 矩阵聚类，取**最平坦/近二倍体/CNV 负荷最低的一簇** spot → `normal_spots.tsv`。用步骤 3 的组织密度做合理性交叉验证（正常参考应落在低密度间质区）。
- **给 stage3 加 `--normal_spots <tsv>` 接口**：用这些参考 spot 的 per-bin 中位作参考（替代当前"全 spot 中位"），重算 bulk CNA + 逐 spot CNV。
  - `stage3_bulk_cna.py --normal_spots normal_spots.tsv ...`
  - `stage3_per_spot_cnv.py --normal_spots normal_spots.tsv ...`
- 预期：肿瘤区 CNA 幅度回升、参考区趋平 → 信噪比上升。
- **代码改动点（落地即做，不依赖图像）**：stage3 的 `ref=np.median(libn[present],0)` 改为"若提供 normal_spots 则只在这些 spot 上取中位"。
- **后续 B3**：Stage-4 RNA 到位后，用 RNA 定的间质/免疫 spot 替换 B2 的数据驱动参考（金标准）。

### 5. 验证：组织/密度 ↔ CNA 空间一致性
- 把组织掩膜 + 密度图 与 `spatial_cna_burden_map` / 逐 spot `cnv_burden` 叠加在配准后的原图上。
- 定量：高 CNV-负荷区是否落在组织高密度区（空间相关）；数据驱动正常簇是否落在低密度区（交叉验证 B2 合理性）。
- 产物：**`image_vs_cnaburden_overlay.png`** + `concordance_stats.tsv`。
- **脚本**：`stage5_anchor_and_validate.py`

### 6. 肿瘤内空间异质性（数据驱动）
- 对逐 spot CNV 矩阵做聚类（PCA+KMeans），扣掉数据驱动正常簇后，看剩余 spot 是否有空间结构化的 CNA 亚群。
- 产物：`intratumor_cnv_clusters.png`、`spot_cnv_cluster_labels.tsv`。
- > 注：低纯度下空间亚克隆能否稳健分出仍是开放问题；此步如不显著，如实报告，待 RNA/更高覆盖再议。

---

## 要新建/改的文件
```
stage5_wsi_registration/
  stage5_register.py            # 步骤1-2：读图+配准，出 affine + 网格叠加图
  stage5_annotate_to_spots.py   # 步骤3：GeoJSON/概率图 → spot_histology_label.tsv
  stage5_anchor_and_validate.py # 步骤4-5：normal_spots 提取 + 一致性验证
  run_stage5.sh                 # 串起来
stage3_cna_clone/
  stage3_bulk_cna.py            # +--normal_spots 接口
  stage3_per_spot_cnv.py        # +--normal_spots 接口
```

## 依赖（落数据时加进 env/）
`openslide-python`、`tifffile`、`opencv-python`、`shapely`、`scikit-image`、（手动 landmark 可选 `napari`）；自动标注可选 `HoVer-Net`/`CLAM`（外部）；圈注用 QuPath（外部）。

## 与上下游接口
- **上游**：Stage-3 的逐 spot CNV 矩阵 = 步骤 5/6 的输入。
- **下游**：`normal_spots.tsv` 回灌 Stage-3；`spot_histology_label.tsv` 同时是 **Stage-4（RNA）方差分解**里"肿瘤密度/微环境"标签的金标准校验。

## 产出（交付物）
1. `chip_grid_overlay.png` —— 配准 QC（原图叠 50×50 网格）。
2. `spot_tissue_label.tsv` + `spot_density.tsv` —— spot→组织掩膜/密度（弱证据，辅助）。
3. `normal_spots.tsv`（数据驱动）+ **锚定后重算的 bulk/逐spot CNA**（更干净、幅度更真）。
4. `image_vs_cnaburden_overlay.png` + 一致性统计。
5. 肿瘤内 CNA 异质性图（如显著）。

---

## 落地顺序（情形已定）
- **已做（不依赖图像）**：stage3 已加 `--normal_spots` 接口 + `pick_normal_spots.py` 自动选参考 + signal-collapse 告警。**520 实测结论：内部自锚不适用本切片**（见上），520 的 bulk CNA 仍以未锚定版为准。接口已就绪，等外部正常参考。
- **病理图一到**：跑 `stage5_register.py`（四角/通道→仿射，必要时手动点四角）→ `stage5_annotate_to_spots.py`（组织掩膜+密度）→ `stage5_anchor_and_validate.py`（叠加验证）。配准+组织密度图是确定能拿到的；肿瘤硬标签不强求。
- **Stage-6 对照 / Stage-4 RNA 一到**：用外部 PoN 或 RNA 间质标签经同一 `--normal_spots` 接口做真正的 normal-anchor —— 这才是本低纯度切片提升 CNA 的正确杠杆。

> 唯一可能要你做的人工动作：若四角自动检测不准，在 thumbnail 上点选捕获区 4 个角（一次，几秒）。
