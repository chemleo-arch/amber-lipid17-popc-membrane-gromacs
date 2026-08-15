#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 在 GROMACS 拓扑的脂质 moleculetype 后插入重原子位置约束（#ifdef POSRES）
# 用法: 在 topol.top 同目录执行 `python add_posre.py`
# 重原子判定: 原子类型不以 h 开头（实测 POPC = 52 重原子/脂, fc=1000）
# 注意: 运行前先备份原拓扑（cp topol.top topol_noposre.bak.top）
import re

lines = open('topol.top').read().split('\n')

atoms_start = None
for i, l in enumerate(lines):
    if l.strip().startswith('[ atoms ]'):
        atoms_start = i
        break
assert atoms_start is not None, '未找到 [ atoms ] 段'

j = atoms_start + 1
heavy = []
last_atom_line = atoms_start
while j < len(lines):
    s = lines[j].strip()
    if s.startswith('['):
        break
    if s and not s.startswith(';'):
        parts = s.split()
        if len(parts) >= 6 and parts[0].isdigit():
            atomnr = int(parts[0])
            atype = parts[1]
            if not atype.lower().startswith('h'):
                heavy.append(atomnr)
            last_atom_line = j
    j += 1
print('heavy atoms in lipid:', len(heavy))

block = ['#ifdef POSRES', '[ position_restraints ]', '; ai  funct  fcx  fcy  fcz']
for a in heavy:
    block.append('%6d  1  1000 1000 1000' % a)
block.append('#endif')

newlines = lines[:last_atom_line + 1] + block + lines[last_atom_line + 1:]
open('topol.top', 'w').write('\n'.join(newlines))
print('inserted posre block, total lines', len(newlines))
