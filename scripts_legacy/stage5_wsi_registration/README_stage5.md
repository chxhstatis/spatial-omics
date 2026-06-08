# Stage 5：H&E ↔ DBiT 网格配准（chip-in-place ROI 框）

把"芯片在位"的 H&E 概览图(带绘制的矩形 ROI 框,标出 50×50 捕获区)配准到 spot 网格,并读出每 spot 的组织密度。

## 运行
```bash
python3 stage5_register.py \
  --he <sample>.jpg \
  --coverage <stage3_out_perspot>/<sample>/per_spot_cnv_burden.tsv \
  --sample <id> --out <out>/<id>
```
- `--he`：H&E 概览 jpg（`.kfb` 全片是江丰专有格式,openslide 读不了;jpg 概览足够配准）。
- `--coverage`：任意含 `x_id,y_id,total_frags` 的表（用作朝向校验的 DNA 覆盖）。

## 方法
1. **检测 ROI 框**：按颜色(红/蓝/绿)取掩膜 → 最大轮廓 → `cv2.minAreaRect` 得 4 角（**支持旋转框**,如样本 A 的蓝色斜框）。
2. **网格→图像仿射/透视**：用 4 角把单位网格 [0,1]² 映到 ROI 四边形。
3. **朝向消歧（8 选 1）**：网格贴框有 4 旋转×2 翻转 = 8 种朝向,无 fiducial 无法先验确定。**尝试用数据定向**——正确朝向应让"H&E 组织包络"与"DNA 覆盖包络"最相关（先各自中值平滑以压制 DBiT 通道条纹）。取相关最高者；并报告 `orientation_resolved`（r>0.25 且与次优差>0.1 才算 RESOLVED）。

## ⚠️ 本数据的实测限制（务必知晓）
- **几何配准稳健**：3 个样本(520/525/A)的网格都正确覆盖 ROI 框,`chip_grid_overlay.png` 可直接用于可视化（朝向无关）。
- **spot 朝向无法从数据自动确定**：DBiT 每-spot 覆盖被**微流控通道效率条纹**主导,不反映组织；且 520/525 组织几乎铺满 ROI 框 → 组织密度近均匀。结果 8 朝向相关都很弱（520/525 r≈0.05–0.17,margin≈0.003 → AMBIGUOUS）。A 因组织未铺满框有部分信号（r≈0.50）但残留 180° 翻转歧义。
- **结论**：要定死"哪个角是 (x_id=1,y_id=1)、x 沿哪方向",需 **fiducial 或实验方的芯片坐标约定**（用户知道 stage2 的 X/Y barcode 与物理芯片的对应）。`spot_density.tsv` 依赖朝向,定向后才最终化；`chip_grid_overlay.png` 不依赖朝向,现可用。

## 产物
`spot_to_image_perspective.npy`(3×3)、`spot_pixel_coords.tsv`(x_id,y_id,px,py,tissue_density)、`spot_density.tsv`、`chip_grid_overlay.png`(H&E+网格+框,QC)、`density_vs_coverage.png`(校验:密度 vs 覆盖网格图+散点)、`register_summary.tsv`(框颜色/角/朝向/校验 r/all_orientation_r)。

## 下一步（定向后）
- `stage5_annotate_to_spots.py`：组织掩膜 + 密度（弱证据,辅助）。
- `stage5_anchor_and_validate.py`：叠加 CNA 负荷做一致性验证。
- 注意:B_B 无 H&E,本 stage 仅 520/525/A。
