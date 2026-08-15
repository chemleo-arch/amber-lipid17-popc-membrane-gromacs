# amber-lipid17-popc-membrane-gromacs

免 CHARMM-GUI 构建纯 POPC 双层膜并转 GROMACS 的 Agent Skill（LIPID17 力场）。

## 用途

当需要搭建 POPC/脂质双层膜、准备膜-配体/膜-蛋白插入模拟的前置膜体系，或把 AMBER 膜拓扑转成 GROMACS 运行时，使用本技能。全管线本地执行，无需网络、无需 CHARMM-GUI 账号，力场与 AMBER 家族（GAFF2 配体 + TIP3P 水）天然一致。

## 管线概览

1. AmberTools 自带 packmol-memgen（python2.7 + PYTHONPATH）构建双层
2. `--parametrize` 走 tleap + LIPID17，产出 AMBER prmtop/inpcrd
3. parmed Python API 转 GROMACS top/gro
4. 插入脂质重原子位置约束（`#ifdef POSRES`）
5. 半各向异性 C-rescale 平衡（EM → NVT → NPT → NPT free）

## 目录结构

```
.
├── SKILL.md                     # 技能主文件
├── reference.md                 # 详细操作说明与踩坑记录
├── scripts/
│   ├── build_membrane.sh        # packmol-memgen 构建脚本
│   ├── amber2gmx.py             # parmed API 转换脚本
│   └── add_posre.py             # 脂质重原子 posre 插入脚本
├── mdp/
│   ├── em.mdp
│   ├── nvt_mem.mdp
│   ├── npt_mem.mdp
│   └── npt_mem_free.mdp
└── templates/
    └── eq_mem.sbatch            # Slurm 平衡作业模板
```

## 安装为技能

将本仓库目录复制到 `~/.qwenworkcn/skills/amber-lipid17-popc-membrane-gromacs/` 即可。

## 关键陷阱（实测修正）

- pip 装不到 packmol-memgen，需用 AmberTools 内置版并设 `PYTHONPATH`
- `--distxy_fix` 只收单值（`70` 而非 `70 70`）
- parmed heredoc 经 ssh 会被吞，改用 Python API
- 膜平衡必须半各向异性 + C-rescale + 脂质重原子 posre + `comm-mode=none`
