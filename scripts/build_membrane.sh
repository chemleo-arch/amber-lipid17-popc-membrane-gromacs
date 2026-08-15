#!/bin/bash
# ============================================================
# 免 CHARMM-GUI 构建纯 POPC 双层（LIPID17 参数化）
# 用法：修改下方 AMBERHOME 为实际路径后，直接 bash build_membrane.sh
# 产出：popc_bilayer_lipid.top (AMBER prmtop) / .crd (inpcrd) / .pdb
# 实测：AmberTools 2018 内置 packmol-memgen 0.9.9（python2.7）
# ============================================================
set -e

# TODO: 改成你的 AmberTools 安装路径（例：/public/software/apps/amber/2018/hpcx-2.4.1-gcc-7.3.1）
export AMBERHOME=/path/to/amber
export PATH=$AMBERHOME/bin:$PATH
export PYTHONPATH=$AMBERHOME/lib/python2.7/site-packages:$PYTHONPATH

cd "$(dirname "$0")" || exit 1

echo "== 环境自检 =="
packmol-memgen --help >/dev/null 2>&1 && echo "packmol-memgen OK" || { echo "packmol-memgen 不可用，检查 AMBERHOME/PYTHONPATH"; exit 1; }

echo "== 构建 POPC 双层 ~7nm (distxy_fix 单值!) =="
# --distxy_fix 只收单值（70 = x/y 各 70 Å）；--parametrize 走 tleap+LIPID17
packmol-memgen --lipids POPC --distxy_fix 70 --leaflet 20 \
    --parametrize --salt --saltcon 0.15 \
    --output popc_bilayer > memgen.log 2>&1

echo "memgen exit=$?"
echo "== 产物 =="
ls -lh popc_bilayer*.top popc_bilayer*.crd popc_bilayer*.pdb 2>/dev/null || true
echo "DONE: 下一步用 scripts/amber2gmx.py 转 GROMACS"
