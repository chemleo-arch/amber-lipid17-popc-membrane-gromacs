---
name: amber-lipid17-popc-membrane-gromacs
description: 免 CHARMM-GUI 构建纯 POPC 双层膜并转 GROMACS 的完整管线：AmberTools 自带 packmol-memgen（LIPID17 参数化）→ parmed API 转 top/gro → 脂质重原子 posre → 半各向异性 C-rescale 平衡。当用户需要搭建 POPC/脂质双层膜、准备膜-配体/膜-蛋白插入模拟的前置膜体系、或把 AMBER 膜拓扑转成 GROMACS 运行时使用。触发词：POPC、双层膜、Lipid17、packmol-memgen、膜平衡、膜搭建、membrane bilayer、GROMACS 膜。
version: 1.0.0
---

# 免 CHARMM-GUI 构建纯 POPC 双层并转 GROMACS（LIPID17）

## 适用场景
需要 GROMACS 脂质双层（纯 POPC 起步，可推广到混合脂）作为膜-配体/膜-蛋白插入模拟的宿主，但无法/不想用 CHARMM-GUI 时。本管线全本地运行，力场与 AMBER 家族（GAFF2 配体 + TIP3P 水）天然一致，适合"AMBER 参数化配体 + LIPID17 膜"的插入模拟。

## 前置条件
- 已安装 AmberTools（实测 AmberTools 2018 版自带 packmol-memgen，python2.7）
- 已安装 GROMACS（本技能 mdrun/grompp 命令以 gmx 或 gmx_mpi 为例）
- 无需网络、无需 CHARMM-GUI 账号、无需 pip 安装

## 工作流

### 1. 环境准备（关键：PYTHONPATH）
packmol-memgen **无法用 pip 安装**（实测 `pip install packmol-memgen` 失败），它内置在 AmberTools 的 python2.7 site-packages 里：

```bash
export AMBERHOME=/path/to/amber        # 按实际安装路径修改
export PATH=$AMBERHOME/bin:$PATH
export PYTHONPATH=$AMBERHOME/lib/python2.7/site-packages:$PYTHONPATH
packmol-memgen --help   # 能出帮助即环境就绪
```

注意：脚本 shebang 是 `#!/usr/bin/python`（python2.7）。若 `module load amber` 后 AMBERHOME 未设置，需手动 export。AmberTools 自带的脂质 PDB 库在 `$AMBERHOME/lib/python2.7/site-packages/packmol_memgen/data/pdbs/`（59 种脂质，含 POPC.pdb）。

### 2. packmol-memgen 构建双层（LIPID17 参数化）
```bash
packmol-memgen --lipids POPC --distxy_fix 70 --leaflet 20 \
    --parametrize --salt --saltcon 0.15 --output popc_bilayer
```
参数要点：
- `--lipids POPC`：脂质类型；混合脂用冒号分隔（如 `POPC:CHL1`）
- `--distxy_fix 70`：膜 x/y 尺寸（Å），**只收单值**（同时设定 x、y）；数字越大脂质越多
- `--leaflet 20`：leaflet 宽度（Å，默认 23）
- `--parametrize`：调用 tleap 参数化，脂质走 **LIPID17**（有蛋白时蛋白走 ff14SB）
- `--salt --saltcon 0.15`：加 0.15 M 盐（KCl），自动中和
- 水模型默认 TIP3P

产出（`--output popc_bilayer`）：`popc_bilayer.pdb`（打包后坐标）、`popc_bilayer_lipid.top`（AMBER prmtop，头为 `%VERSION VERSION_STAMP`）、`popc_bilayer_lipid.crd`（inpcrd）、`popc_bilayer_lipid.pdb`。
实测 70 Å 盒 ≈ 154 POPC（77/leaflet）+ 5736 WAT + 14 K+ + 14 Cl- ≈ 3.8 万原子。packmol 打包含放水，需数分钟，耐心等待 memgen 日志出现 "DONE!"。

### 3. parmed 转 GROMACS（用 Python API，不要 heredoc）
**不要**用 `parmed xxx << EOF` 交互 heredoc 方式——经 ssh 传输会被吞、输出为空（实测踩坑）。改用 Python API 脚本（见 scripts/amber2gmx.py）：

```python
import parmed as pmd
s = pmd.load_file('popc_bilayer_lipid.top', xyz='popc_bilayer_lipid.crd')
s.save('popc_gmx.top', overwrite=True)
s.save('popc_gmx.gro', overwrite=True)
```

运行前同样要 `export PYTHONPATH=$AMBERHOME/lib/python2.7/site-packages`（parmed 是 python2.7 包）。转换后建议核验：总电荷≈0、盒子尺寸、残基统计。
要点：LIPID17 中每条 POPC 拆成 **PA/OL/PC 三个残基**（palmitoyl/oleoyl/phosphocholine 头），parmed 会按键连关系合并为单个 moleculetype（实测命名 `system1`，138 原子/脂）——不要把三残基当成三个独立分子。

### 4. 插入脂质重原子位置约束
在拓扑第一个 `[ atoms ]`（脂质 moleculetype）之后插入 `#ifdef POSRES` 的 `[ position_restraints ]` 块，fc=1000。重原子判定：原子类型不以 h 开头（实测 POPC 每脂 52 重原子 = C42+N1+O8+P1）。直接用 scripts/add_posre.py（已实测）：
```bash
cp topol.top topol_noposre.bak.top   # 先备份
python add_posre.py
```
注意：拓扑里已有 `[ exclusions ]`，脚本只插入 posre 块、不动其他；所有脂质同构，一份 posre 块全局有效。

### 5. 平衡 mdp（半各向异性 + C-rescale）
四阶段（mdp 文件见 mdp/ 目录，sbatch 模板见 templates/eq_mem.sbatch）：
1. EM（steep，50k 步）
2. NVT 100 ps：`define = -DPOSRES`，V-rescale 300 K
3. NPT 100 ps：`define = -DPOSRES`，**`pcoupl = C-rescale` + `pcoupltype = semiisotropic`**（膜必须半各向异性：xy 平面与 z 方向独立控压）
4. NPT 无约束 5 ns（起步）：同上半各向异性，去掉 POSRES

膜平衡关键设置：`compressibility = 4.5e-5 4.5e-5`、`ref_p = 1.0 1.0`、`tau_p = 2.0`、PME rcoulomb=1.2、rvdw=1.2、DispCorr=EnerPres、constraints=h-bonds，且三阶段 mdp 都要 **`comm-mode = none`**（有位置约束时必须关闭质心移除，否则 grompp 报 COM 警告且可能引入伪影）。

### 6. 验证与警告判读
- 逐级 `grompp`（EM→NVT→NPT→NPT free）通过
- **"net charge -0.000314" 警告：浮点舍入，无害**（AMBER 电荷精度导致，体系实为中性）
- "Center of mass motion removal" 警告：加上 `comm-mode = none` 后消失
- 残留 1 NOTE + 1 WARNING（电荷）不致命，`-maxwarn 2` 可过

## 关键陷阱（非显而易见修正，均来自实测）
1. pip 装不到 packmol-memgen → 用 AmberTools 内置版（python2.7 + 手动 PYTHONPATH）
2. `--distxy_fix` 只收**单值**：写 `70`，写 `70 70` 会报 "unrecognized arguments: 70"
3. parmed heredoc 经 ssh 被吞 → 用 Python API 脚本（scp 上传后执行）
4. 膜平衡必须 semiisotropic + C-rescale + 脂质重原子 posre + comm-mode=none，缺一不可
5. LIPID17 脂质是多残基表示（PA/OL/PC），parmed 后合并为一个 moleculetype
6. 净电荷警告是 AMBER 电荷舍入，非真带电
7. 含中文注释的脚本在 python2 下跑需加 `# -*- coding: utf-8 -*-` 编码声明

## 验证清单
- [ ] `packmol-memgen --help` 正常输出（PYTHONPATH 生效）
- [ ] 产出 prmtop(.top 头为 %VERSION) + inpcrd(.crd)
- [ ] parmed 转换后 EM grompp 通过
- [ ] posre 块插入后 NVT grompp（-DPOSRES）通过，重原子数正确（POPC=52/脂）
- [ ] NPT 半各向异性跑通，盒子 x/y 与 z 独立变化
- [ ] grompp 无 COM 移除警告（comm-mode=none 已生效）

## 辅助文件索引
- `reference.md`：详细操作说明与参数解释
- `scripts/build_membrane.sh`：packmol-memgen 构建脚本模板（改 AMBERHOME 即用）
- `scripts/amber2gmx.py`：parmed API 转换脚本
- `scripts/add_posre.py`：脂质重原子 posre 插入脚本（实测版）
- `mdp/em.mdp`、`mdp/nvt_mem.mdp`、`mdp/npt_mem.mdp`、`mdp/npt_mem_free.mdp`：四阶段平衡 mdp（含 comm-mode=none）
- `templates/eq_mem.sbatch`：Slurm 平衡作业模板（分区/资源参数按集群修改）
