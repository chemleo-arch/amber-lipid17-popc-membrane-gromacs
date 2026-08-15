#!/usr/bin/env python
# -*- coding: utf-8 -*-
# AMBER prmtop/inpcrd -> GROMACS top/gro（parmed Python API）
# 用法: python amber2gmx.py <prmtop> <inpcrd> <out.top> <out.gro>
# 运行前: export PYTHONPATH=$AMBERHOME/lib/python2.7/site-packages
# 注意: 不要用 parmed 交互 heredoc（经 ssh 会被吞），一律用本脚本
import sys
import parmed as pmd
from collections import Counter

top, crd, out_top, out_gro = sys.argv[1:5]

s = pmd.load_file(top, xyz=crd)
print('ATOMS', len(s.atoms), 'RESIDUES', len(s.residues))
c = Counter(r.name for r in s.residues)
print('RES COUNTS', dict(c.most_common(8)))
qtot = sum(a.charge for a in s.atoms)
print('TOTAL_CHARGE %.3f' % qtot)
try:
    print('BOX', s.box)
except Exception as e:
    print('box err', e)

s.save(out_top, overwrite=True)
s.save(out_gro, overwrite=True)
print('SAVED', out_top, out_gro)
