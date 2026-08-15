# reference：免 CHARMM-GUI 纯 POPC 双层构建与 GROMACS 转换（LIPID17）

本文件是 SKILL.md 的展开说明，记录完整操作细节、参数解释与实测踩坑记录。管线在 HPC 集群（远程 ssh + Slurm）上实测通过；本机有 AmberTools 时同样适用。

## 背景与动机

膜-配体插入模拟需要先有一个平衡好的脂质双层。常规做法是 CHARMM-GUI Membrane Builder（需注册、交互式网页），或本地 packmol/insane（还要自己解决力场参数化）。本管线的优势：

- 全本地/集群执行，无网络依赖
- packmol-memgen 一次性完成"打包 + 水/离子 + tleap 参数化"
- 力场统一走 AMBER 家族：配体 GAFF2、膜 LIPID17、水 TIP3P——三者参数天然一致，避免混合力场问题

## 实测环境

- AmberTools 2018（HPC 版），packmol-memgen 版本 0.9.9（2018 内置）
- GROMACS 2023.2（DCU 移植版，`gmx_mpi` 二进制）
- 系统 python2.7 驱动 packmol-memgen 与 parmed

## 详细步骤

### Step 0：工具定位（为什么找不到 packmol-memgen）

`pip install packmol-memgen` 会失败——该工具不发布到 PyPI（在 conda-forge 或随 AmberTools 分发）。AmberTools 安装目录下已有：

```
$AMBERHOME/bin/packmol-memgen      # 主脚本（shebang: #!/usr/bin/python，即 python2.7）
$AMBERHOME/bin/packmol              # packmol 二进制
$AMBERHOME/lib/python2.7/site-packages/packmol_memgen/
    data/pdbs/POPC.pdb              # 脂质模板库，59 种脂质
    amber.py / example/example.sh   # 官方示例
```

不设 PYTHONPATH 直接跑会报 `No module named packmol_memgen`。

### Step 1：构建命令与参数

```bash
export AMBERHOME=/public/software/apps/amber/2018/hpcx-2.4.1-gcc-7.3.1   # 示例路径
export PATH=$AMBERHOME/bin:$PATH
export PYTHONPATH=$AMBERHOME/lib/python2.7/site-packages:$PYTHONPATH

packmol-memgen --lipids POPC --distxy_fix 70 --leaflet 20 \
    --parametrize --salt --saltcon 0.15 --output popc_bilayer
```

| 参数 | 取值 | 说明 |
|---|---|---|
| `--lipids` | `POPC` | 冒号分隔可混合，如 `POPC:CHL1` |
| `--distxy_fix` | `70` | 膜 x/y 尺寸（Å），**单值**同时设定两方向 |
| `--leaflet` | `20` | leaflet 宽度（Å），默认 23 |
| `--parametrize` | — | tleap 参数化：脂质 LIPID17，蛋白 ff14SB |
| `--salt` / `--saltcon` | `0.15` | 0.15 M KCl，自动中和 |
| `--output` | `popc_bilayer` | 输出名，packmol PDB 与参数化产物共用前缀 |

产出文件：

```
popc_bilayer.pdb             # packmol 打包结果（水/离子已加入，参数化前）
popc_bilayer_lipid.top       # AMBER prmtop（头两行 %VERSION VERSION_STAMP）
popc_bilayer_lipid.crd       # AMBER inpcrd（坐标）
popc_bilayer_lipid.pdb       # 参数化后坐标
```

实测 70 Å 盒：154 POPC（77/leaflet）+ 5736 WAT + 14 K+ + 14 Cl-，37,872 原子。打包（放水）较慢，memgen 日志出现 "DONE!" 才算完。

### Step 2：parmed 转换（避免 heredoc 陷阱）

交互 heredoc 在 ssh 里会出问题（实测 `parmed top crd << EOF ... EOF` 输出为空/被吞）。正确方式：写 Python API 脚本，scp 上传后执行。

```python
import parmed as pmd
s = pmd.load_file('popc_bilayer_lipid.top', xyz='popc_bilayer_lipid.crd')
print('ATOMS', len(s.atoms), 'RESIDUES', len(s.residues))
from collections import Counter
print(dict(Counter(r.name for r in s.residues).most_common(8)))
print('TOTAL_CHARGE', sum(a.charge for a in s.atoms))
s.save('popc_gmx.top', overwrite=True)
s.save('popc_gmx.gro', overwrite=True)
```

注意：parmed 同样是 python2.7 包，运行前保持 `PYTHONPATH` 设置。实测版本 3.0.0。

转换后核验项：
- 总电荷 ≈ 0（实测 -3.1e-4 量级，舍入）
- 残基统计：OL/PC/PA 各 154（=154 条 POPC），WAT 5736，K+/Cl- 各 14
- 盒 7.7×7.8×8.0 nm
- GROMACS 拓扑里 `[ moleculetype ]`：`system1`（脂质，138 原子/分子）+ K+ + Cl- + WAT

### Step 3：脂质重原子 posre

目标：约束脂质重原子（轻氢不约束）防止膜结构在平衡早期散架。做法：在拓扑第一个 `[ atoms ]`（属于脂质 moleculetype）块后插入：

```
#ifdef POSRES
[ position_restraints ]
; ai  funct  fcx  fcy  fcz
<重原子行：%6d  1  1000 1000 1000>
#endif
```

用 add_posre.py 自动完成：解析第一个 `[ atoms ]`，原子类型不以 h 开头者判为重原子。实测 POPC = 52 重原子/脂（C42 + N1 + O8 + P1）。插入前先备份原拓扑。

### Step 4：平衡四阶段

| 阶段 | mdp | 时长 | POSRES | 控压 |
|---|---|---|---|---|
| EM | em.mdp | 50k 步 steep | 无 | 无 |
| NVT | nvt_mem.mdp | 100 ps | `-DPOSRES` | 无 |
| NPT | npt_mem.mdp | 100 ps | `-DPOSRES` | Berendsen semiisotropic（tau_p=5）+ refcoord-scaling=com |
| NPT free | npt_mem_free.mdp | 5 ns（起步） | 无 | Parrinello-Rahman semiisotropic（tau_p=5） |

膜平衡特有要点：

1. **半各向异性控压**：`pcoupltype = semiisotropic`，`compressibility = 4.5e-5 4.5e-5`，`ref_p = 1.0 1.0`，`tau_p = 5.0`。各向同性会错误约束膜面积；semiisotropic 让 xy 面积与 z 厚度独立响应。**受限 NPT 用 `Berendsen`（弱耦合、不震荡）并设 `refcoord-scaling = com`**（否则压耦 + 绝对位置约束打架），**无约束生产段用 `Parrinello-Rahman`**。切勿在初始膜 + posres 下用 C-rescale：本会话实测盒子单步缩放 mu 冲到 1.45，拉坏坐标触发 DCU VMFault 崩作业。
2. **`comm-mode = none`**：三阶段 mdp 都要。有位置约束时质心移除（默认）会与约束冲突，grompp 会警告；关闭后警告消失，也避免伪影。
3. 力场一致性：AMBER 参数（gen_pairs=yes、fudgeLJ=0.5、fudgeQQ=0.8333）由 parmed 转换自动带入，无需手改。
4. NVT→NPT 用 `-t` 续接（.cpt），NPT→NPT free 同样。

### Step 5：提交与验证

- `grompp` 每阶段先登录节点试跑，确认只出 1 NOTE + 1 WARNING（电荷舍入）再提交
- mdrun 命令示例（GROMACS 2023+ 语法）：`gmx_mpi mdrun -deffnm nvt -ntomp 8 -nb gpu`
- 平衡完成判据：面积/脂质（APL）收敛、盒 z 稳定

## 陷阱速查表

| # | 现象 | 原因 | 修正 |
|---|---|---|---|
| 1 | `No module named packmol_memgen` | pip 无此包；未设 PYTHONPATH | `export PYTHONPATH=$AMBERHOME/lib/python2.7/site-packages` |
| 2 | `unrecognized arguments: 70` | `--distxy_fix` 只收单值 | 写 `--distxy_fix 70` |
| 3 | parmed 交互无输出 | heredoc 经 ssh 被吞 | 用 Python API 脚本 + scp 上传 |
| 4 | grompp 报 COM 移除警告 | 位置约束 + 默认 comm 模式 | 加 `comm-mode = none` |
| 5 | grompp 报 net charge -0.0003 | AMBER 电荷浮点舍入 | 无害，`-maxwarn 2` |
| 6 | 中文注释脚本 python2 报错 | 编码声明缺失 | 文件头加 `# -*- coding: utf-8 -*-` |
| 7 | Lipid17 残基数吓人（PA/OL/PC 各 154） | 脂质多残基表示 | parmed 已合并为单 moleculetype，勿手动拆分 |
| 8 | 受限 NPT 崩 DCU VMFault（`pressure scaling more than 1%`，mu>1.4） | C-rescale 在初始膜 + posres 下震荡拉坏坐标 | 受限段改 Berendsen（tau_p=5）+ `refcoord-scaling = com`，无约束段 Parrinello-Rahman |
| 9 | genion 报水组不存在或 `Water is not continuous` | 水组应选 `Water`（非 WAT）；水被聚集体切成两段不连续 | 选 `Water`；保证水连续（重排到末尾）或体系中性时跳过 genion |

## 后续衔接（膜-配体插入）

膜平衡完成后，可：
1. 用 packmol/`gmx insert-molecules` 把配体聚集体放到膜上方
2. 合并拓扑：GAFF2 配体 itp + LIPID17 膜 top（两者原子类型零冲突：GAFF2 全小写、LIPID17 大写/混合；`[ defaults ]` 均为 AMBER 参数，可直接拼接）
3. 再次 posre + 半各向异性平衡，再做无偏插入或伞采样 PMF
4. 加离子注意：genion 水组选 `Water`（非 `WAT`），且水必须连续（聚集体别插在膜水与后加水之间）；体系电中性时直接跳过 genion 最稳妥

## 原始来源

本技能整理自 2026-08-14 卤素二肽-膜插入模拟项目（chatId mssitrzx55bt1ak6）的实测流程。命令与脚本均为集群上跑通版本。
