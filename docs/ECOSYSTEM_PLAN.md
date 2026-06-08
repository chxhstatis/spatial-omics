# 长期生态方案：把空间 DNA-seq 流程做成"Seurat 级"的开放工具 + 平台推广引擎

**日期：** 2026-06-06
**对象：** 空间 DNA-seq（DBiT 50×50 → 拷贝数/克隆，对标 slide-DNA-seq, Nature 2022）
**现状基线：** `spatial_omics_pipeline/`（config 驱动、stage1–6、内置参考、spot×bin 稀疏矩阵数据模型、三个算法卖点）
**目标：** 让"买/用你这套测序技术的人"当天就能装上工具、跑通教程、出可发表的图——从而把"分析门槛"从平台采用的障碍变成你的护城河。

---

## 0. 战略定位（先想清楚这三件事，再动手）

| 问题 | 你的答案决定了 | 默认建议 |
|---|---|---|
| **谁是用户？** | 投入规模 | 早期：合作课题组 + 自己实验室；中期：买你芯片/试剂盒的湿实验用户 |
| **卖什么？** | 商业模式 | 软件**永远免费开源**（=引流），收入来自芯片/试剂盒/测序服务。软件是"杀软不收钱、卖剃须刀片"的剃须刀 |
| **护城河是什么？** | 能不能成 | **不是"又一个流程"**，而是①标准数据对象 ②三个别人没有的算法 ③"空间基因组分析"这个**空着的生态位** |

> 关键判断：scRNA 分析被 Seurat/Scanpy 占死，挤不进。但**"稀疏空间基因组（CNA/克隆）的标准分析包"目前是空白**——官方 slide-DNA-seq 代码是 MATLAB（要付费 license、装不上、没人能跑）。这是你真正的机会窗口。

---

## 1. 产品形态：三层结构（照搬 10X 的成功架构）

10X 不是只有 Seurat，它是**三条腿**。你要复制这个结构：

```
┌─ 第①层  上游引擎 = "Cell Ranger 等价物" ───────────────┐
│  FASTQ → spot×bin 矩阵（stage1–2）                       │
│  形态：一条 nextflow/snakemake 流水线 + Docker 镜像       │
│  交付：一行命令从原始数据到标准对象，HPC/云都能跑          │
├─ 第②层  下游分析包 = "Seurat 等价物" ★核心★ ───────────┤
│  标准数据对象 + 去偏 + 双重平滑 + 克隆 + 显著性 + 方差分解 │
│  形态：pip/conda 可装的 Python 包，import 进任何脚本       │
│  这是教程网站演示的主角，是论文的主角，是引用的来源        │
├─ 第③层  门户网站 = "satijalab.org 等价物" ─────────────┤
│  文档 + 可复制教程(vignette) + demo 数据 + API 手册       │
│  可选：无代码 web app = "Loupe Browser 等价物"            │
└──────────────────────────────────────────────────────────┘
```

**优先级：先做第②层（分析包），它是一切的中心。** 第①层你已经有脚本，后期容器化即可；第③层等第②层有 3 个能跑的教程后再上线。

---

## 2. 标准数据对象（这是整个工具的"宪法"，等同于 SeuratObject / AnnData）

Seurat 成功的隐形原因：**所有人都围绕同一个对象写代码**。你必须先定义一个标准对象，否则永远是散脚本。

**直接复用 `AnnData`（anndata 库）——不要自己造轮子。** 好处：免费继承 scanpy / squidpy 整个生态，用户零学习成本。

映射方式：

```
AnnData 结构              你的空间 DNA-seq 含义
─────────────────────────────────────────────────────────
.X                    →   spot × bin 的拷贝数矩阵（主矩阵）
.obs (行 = 2500 spots)→   x_id, y_id, total_frags, clone_label, cnv_burden
.var (列 = ~2646 bins)→   chr, start, end, gc, mappability, blacklist_flag
.obsm['spatial']      →   50×50 物理坐标（squidpy 直接能画空间图）
.layers['raw_counts'] →   原始 fragment 计数
.layers['gc_corrected']→  去偏后
.layers['smoothed']   →   双重平滑后
.obsp['spatial_knn']  →   xy 空间近邻图
.obsp['pc_knn']       →   PC 基因组相似空间近邻图
.uns['sample']        →   样本元数据、config 快照、参考版本（可溯源）
.uns['permutation']   →   置换显著性 z 分数图
```

**行动**：写一个 `from_pipeline()` 读你现有的 `*_counts.npz` + `matched_spot_barcodes.tsv` → 返回标准 AnnData。这一步把你"散脚本时代"和"产品时代"接上。

---

## 3. 算法差异化：把三个卖点做成"一行 API + 一张惊艳的图"

工具的灵魂不是包装，是**别人 demo 一跑就 wow** 的算法。你笔记里已经识别出这三招，现在要把它们变成干净的 API：

| API | 做什么 | 为什么是卖点 | 现状 |
|---|---|---|---|
| `pp.correct_bias(adata)` | GC + mappability + blacklist + 深度 层层除偏 | 稀疏空间 CNV 的前提，做错就全是伪影 | ✅ 已有(stage3)，需封装 |
| `tl.dual_smooth(adata)` | **xy 物理空间 + PC 基因组空间 双重平滑** | **对抗稀疏的灵魂**，目前只做了 xy，补 PC 即领先 | ⚠️ 只做了一半 |
| `tl.permutation_significance(adata)` | fragment→bead 置换经验零分布 → 空间 p 值图 | "某扩增在组织哪里显著"，比单纯聚类有说服力 | ❌ 待加 |
| `tl.call_clones(adata)` | PCA + kmeans + CH 准则客观选 k | de novo 克隆，无需参考 | ✅ 已有 |
| `tl.variance_decomposition(adata, rna)` | stepwiselm：表达 ~ 亚克隆 + 肿瘤密度 + 免疫密度 | **multi-modal 的最大贡献**，归因内因 vs 外因 | ❌ 待 RNA(stage4) |

> 命名按 scanpy 惯例分 `pp`(预处理) / `tl`(分析) / `pl`(画图) 三个模块——用户秒上手。

**这一层每个函数都要配一张"主图"**：双重平滑前后的克隆边界对比、置换显著性的空间 p 值图、方差分解的内因/外因堆叠条形图。教程和论文就靠这几张图。

---

## 4. 需要补充哪些数据（你明确问的）

分四类，按紧急程度排：

### 🔴 必需：可公开分享的 Demo 数据（没有它，教程网站根本上不了线）
病人数据不能公开 → 教程不能用真实数据。两条路（建议都做）：
1. **合成数据生成器** `datasets.simulate(grid=50, n_clones=3, coverage=0.4)`——程序化造一份带已知克隆结构的假数据，体积小、可随包分发、教程秒跑。**这是最快能让网站上线的方案。**
2. **接公开基准数据**：slide-DNA-seq 原文（Zhao et al., *Nature* 2022）的小鼠/人数据在 GEO 公开。下采样一份做"真实数据教程"，同时**证明你的工具能复现顶刊结果**（极强的信任背书）。

### 🟠 重要：对照/正常基线（stage6 缺的）
现在缺 normal/control 的 CNA 基线，导致拷贝数标定靠"自身众数×2"——有风险。补一份**配套 bulk-WGS 或正常组织**，做 `--control` 通道（你已留接口）。这是从"自用"到"可发表/可信"的关键。

### 🟡 增值：同片 RNA（解锁 stage4 = 方差分解 = 全文灵魂）
没有 RNA，方差分解（你最大的差异化）做不了。补**同一张片子的 spatial RNA**，stage4 才能落地。这决定你的工具是"又一个 CNV 工具"还是"真正的 multi-modal 平台"。

### 🟡 增值：WSI / H&E 图（stage5 配准）
肿瘤区标注 → 锚定验证空间负荷。有它，图的说服力上一个台阶。

### 🟢 基础：参考数据版本化
现在 `ref/` 内置了 hg38 的 GC/mappability/blacklist（很好）。需要：①补 rep-timing track（hg38 可下）；②做成**带版本号、带 checksum** 的可下载参考包（reproducibility 要求）；③预留 mm10/T2T 接口。

---

## 5. 怎么搭建（工程化：从脚本 → 可安装的包）

### 仓库结构（重构现有 `spatial_omics_pipeline/`）
```
spacednatools/                 # ← 新的 Python 包名（见下方命名建议）
├── pyproject.toml             # 现代打包（取代 setup.py）
├── src/spacednatools/
│   ├── io.py                  # from_pipeline(), read/write h5ad
│   ├── pp/                    # correct_bias, qc, normalize  （← stage3 拆进来）
│   ├── tl/                    # dual_smooth, call_clones, permutation_sig, variance_decomp
│   ├── pl/                    # 所有画图（spatial_cn, clone_map, significance, vardecomp）
│   ├── datasets/              # simulate(), 公开数据下载器
│   └── _ref/                  # 参考数据（或运行时下载）
├── workflow/                  # ← stage1-2 的 nextflow/snakemake + Dockerfile（第①层）
├── tests/                     # pytest，每个算法一个单测（保证不回归）
├── docs/                      # 网站源码（见第6节）
│   ├── tutorials/*.ipynb      # 可执行教程 = vignette
│   └── api/                   # 自动从 docstring 生成
├── .github/workflows/ci.yml   # 自动测试 + 自动建文档 + 自动发版
├── CITATION.cff               # 让人怎么引用
├── CHANGELOG.md
└── README.md
```

### 打包与分发渠道（从易到难，逐步上）
1. **PyPI**：`pip install spacednatools` —— 最低门槛，先做这个。
2. **conda / bioconda**：生信用户习惯 conda。进 bioconda channel 还自带"被领域认证"的信号。
3. **Docker / BioContainers**：上游流水线（stage1-2）容器化，推到 quay.io，HPC/云一键拉取。
4. **Zenodo**：每个 release 自动存档 → 拿到 **DOI**，论文可引、永久可复现。

### 质量基建（这是"维护得起"的前提）
- **CI（GitHub Actions）**：每次提交自动跑 pytest + 在合成数据上跑全流程 + 自动重建文档。
- **语义化版本**（v0.1.0 起步）+ CHANGELOG。
- **单元测试**：尤其去偏和置换检验，必须有"已知输入→已知输出"的测试，防止改坏。
- **类型注解 + docstring**：docstring 直接生成 API 文档，一份功夫两用。

---

## 6. 怎么"在线"（门户网站 = satijalab.org 等价物）

### 推荐技术栈（小团队友好、免费、好看）
- **文档引擎**：`mkdocs` + `mkdocs-material` 主题（Seurat 那种干净观感）+ `mkdocs-jupyter`（让 .ipynb 教程直接渲染成网页）。
  - 备选：`sphinx` + `myst-nb`（scanpy/squidpy 用这套，更强但更重）。**小团队选 mkdocs-material。**
- **托管**：**GitHub Pages 免费**（`mkdocs gh-deploy` 一行部署），绑自定义域名。
  - 备选：ReadTheDocs（自动从每次 git push 重建，带版本切换）。
- **教程 = 可执行 notebook**：每个 vignette 是一个能从头跑到尾的 `.ipynb`，CI 里自动执行确保"教程永远是能跑通的"——这正是 Seurat 网站的精髓（不是介绍，是能复制粘贴的真实分析）。

### 网站必备的页面（按 Seurat 的成功结构）
1. **首页**：一句话定位 + 一张最惊艳的图（双重平滑的克隆图或方差分解图）+ 安装命令。
2. **Installation**：`pip` / `conda` / `docker` 三种，30 秒装上。
3. **Tutorials（核心）**：
   - *Getting started*：合成数据 10 分钟跑通（QC→去偏→克隆→出图）。
   - *Spatial copy number & clones*：你的主交付物。
   - *Permutation significance*：空间显著性图。
   - *Multi-modal variance decomposition*：RNA 到位后的旗舰教程。
   - *Reproducing slide-DNA-seq (Nature 2022)*：用公开数据复现顶刊 → 信任背书。
4. **From FASTQ（上游）**：nextflow 流水线怎么跑。
5. **API reference**：自动生成。
6. **About / Cite / 联系**：论文、DOI、平台介绍（这里软性推广你的测序技术）。

### 可选第三层：无代码 Web App = "Loupe Browser 等价物"
湿实验用户不写代码。可以用 `streamlit` 或 `shiny for python` 做一个**上传 h5ad → 交互看空间克隆图/显著性图**的网页，部署到 Streamlit Community Cloud（免费）或你自己的服务器。**这一步等核心包稳定后再做，是放大采用率的利器，但不是 MVP。**

---

## 7. 论文与引用策略（让工具"可被引用"才有生命力）

- **方法学预印本**：写好包 + 教程后，发 **bioRxiv 预印本**（先占位、先可引），再投 *Genome Biology* / *Bioinformatics* / *Nature Methods*（看创新度）。
- **CITATION.cff + Zenodo DOI**：每个 release 可引。
- **论文的卖点**：不是"我们做了个流程"，而是"我们提出**双重平滑 + 置换显著性 + 方差分解**解决稀疏空间基因组的三大难题，并提供首个 pip 可装的开放工具"。把算法当主角。
- **同时推平台**：方法学论文 + 一篇用该平台+工具的生物学应用论文（如你的 PDAC 课题），双论文互引，平台和工具一起立住。

---

## 8. 分阶段路线图（务实，约 12–18 个月）

| 阶段 | 时间 | 交付 | 里程碑 |
|---|---|---|---|
| **P0 内核** | 1–2 月 | 定义 AnnData 对象 + `from_pipeline()` + 把 stage3 拆成 `pp`/`tl`/`pl` API | `import spacednatools` 能跑通自己的数据 |
| **P1 算法补全** | 2–3 月 | 补 PC 双重平滑 + 置换显著性；合成数据生成器；pytest | 三个卖点都有 API + 主图 |
| **P2 上线 MVP** | 1 月 | PyPI 发包 + mkdocs 网站 + 3 个教程上 GitHub Pages | **网址能访问、`pip install` 能装、教程能跑** |
| **P3 可信背书** | 2 月 | 复现 slide-DNA-seq 公开数据教程 + control 基线（补数据后） | "能复现 Nature 2022"写进网站 |
| **P4 multi-modal** | 2–3 月 | 补同片 RNA → stage4 方差分解落地 + 旗舰教程 | 从"CNV 工具"升级为"multi-modal 平台" |
| **P5 论文 + 上游** | 2–3 月 | nextflow 上游容器化 + bioRxiv 预印本 + Zenodo DOI + bioconda | 可引用、可一键全流程、被领域看见 |
| **P6 放大** | 持续 | 无代码 web app + WSI 配准 + 多基因组 + 社区 | 网络效应启动 |

---

## 9. 维护与社区（最被低估、也最决定生死）

- **这是长期投入，不是做完就完**：开源工具死于"作者毕业/换方向后无人维护"。要有**明确的维护人 + 接班机制**。
- **GitHub Issues 响应**：用户提问 48h 内回——早期口碑全靠这个。
- **每季度一个 release** + 更新日志，让人感觉"活着"。
- **示范用户**：先拉 2–3 个合作课题组用起来、出论文、互相引用——网络效应的火种。
- **诚实的局限声明**（你 README 里已经做得很好）：~0.4× 覆盖只能做 CNA/克隆、不能做 SNV；低纯度限制——**写进文档**，过度宣称会反噬信任。

---

## 10. 成本与团队（诚实评估）

| 资源 | 最低配置 | 理想配置 |
|---|---|---|
| 人 | 1 个会 Python 打包 + 生信的人，全职 6 个月做出 MVP | 1 开发 + 1 算法 + 1 文档/支持 |
| 钱 | 几乎零（GitHub/PyPI/Pages/Zenodo 全免费）；域名 ~¥100/年 | + web app 服务器、设计 |
| 数据 | 合成数据（免费）+ 公开 slide-DNA-seq | + 自产 RNA/control/WSI 配套数据 |

> 真正的成本不是钱，是**持续的人力**和**一个真能打动人的算法+数据类型**。如果只是自己实验室用，做到 P0–P1（内核+API）即可，不必上网站/发版。**只有当你确实要靠它推广平台、有外部用户时，P2 之后的投入才值得。**

---

## 一句话总结
> 把 `pp.correct_bias → tl.dual_smooth → tl.call_clones → tl.permutation_significance → tl.variance_decomposition` 这条链做成一个 **AnnData 原生、pip 可装、教程能跑通** 的包，用合成数据 + 复现 Nature 2022 建立信任，配一个 mkdocs 网站和一篇方法学论文——这就是空间基因组版的 Seurat，也是你测序平台的采用引擎。

---

## 附录 A：竞品生态位调研结论（2026-06-06）

完整调研要点（含来源 URL）：

1. **空白真实且干净**：现有所有 CNV/克隆工具（inferCNV、CopyKAT、Numbat、STARCH、SpatialInferCNV、infercnvpy）都是**从 RNA 表达反推 CNV**，没有一个为「测得的 DNA reads」设计。→ **这就是我们的差异化定位：直接分析 DNA。**
2. **官方 slide-DNA-seq 代码已废弃**：`buenrostrolab/slide_dna_seq_analysis` 是 MATLAB（87%），2021 后无更新、4 stars、不能 pip 装、无教程。这正是要替代的对象。
3. **inferCNV 已停止维护**（README 明示），CopyKAT 停滞（2022），Numbat 活跃但基于 RNA/allele。
4. **数据对象**：事实标准在收敛到 scverse 体系（`AnnData` → `SpatialData`）。**结论：站在 scverse 肩上，做成兼容插件**（我们的包已是 AnnData 原生）——而不是自造对象。
5. **对标剧本 = Stereo-seq + Stereopy**（华大平台方自养小团队做工具+readthedocs+Nat Commun 论文推自家芯片）。**MVP 起步 1–3 人即可。**
6. **公开基准数据**：slide-DNA-seq 原始数据在 **SRA `PRJNA768453`**（Zhao et al. Nature 2021），可下载、可做"复现顶刊"旗舰教程。
7. **命名**：PyPI 实测 `spatialcnv` / `karyospace` / `clonespace` / `cnvspace` 均空闲；避开 `spatula`、`spclone`（已被无关项目占用）。工作名用了 `spacecna`（待最终定名）。

**值不值得做（一句话判断）**：值得，但要当成「平台市场赌注的软件杠杆」而非独立软件创业——技术空白干净、1–3 人可做出事实标准 MVP；**真正的成败不取决于代码，而取决于这套 DBiT 空间 DNA 芯片能否卖出装机量**。最大风险是平台采用度，其次是被 scverse 官方顺手吃掉（对策：从第一天就做成 scverse 兼容插件）。

---

## 附录 B：P0 内核已落地（2026-06-06）

**已搭好可运行的包骨架**，位置：`spatial_pdac/spacecna/`（工作名 `spacecna`）。

- ✅ **标准数据对象**：AnnData 原生（`io.from_pipeline` 接你现有 stage1-2 产物 / `io.from_matrix` 接矩阵）。
- ✅ **scanpy 风格 API**：`pp`（bin_qc / correct_bias / normalize）、`tl`（dual_smooth / call_clones / copy_number / permutation_significance）、`pl`（4 种图）、`datasets.simulate`。
- ✅ **三个算法卖点全部移植自你的 stage3**（双重平滑、置换显著性已在你代码里实现，本次封装进 API）。修复了一个稀疏致命 bug：~0.4× 覆盖下伪正常参考的跨 spot 中位数恒为 0，必须先做空间致密化（已加 `normalize` 的 `spatial_k`）。
- ✅ **合成数据生成器**：无需病人数据即可跑通全流程、上教程、做单测。
- ✅ **验证**：合成数据克隆恢复 **ARI ≈ 0.96**，chr8 增益、无深度混杂；6 个单测全过；`examples/quickstart.py` 出 4 张图。
- ✅ **工程化**：`pyproject.toml`、pytest、mkdocs-material 网站（首页/安装/3 教程/自动 API）、GitHub Actions CI（多 Python 版本测试 + 自动跑教程 + 自动发 Pages）、README/CHANGELOG/LICENSE/CITATION/.gitignore。

**下一步（P1→P2）**：补 RNA 后做 `tl.variance_decomposition`（stage4 灵魂）；接 `PRJNA768453` 公开数据做复现教程；定名 + `git init` 推 GitHub + 发 PyPI + `mkdocs gh-deploy` 上线网址。
