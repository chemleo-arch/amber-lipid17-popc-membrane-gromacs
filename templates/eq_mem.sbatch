#!/bin/bash
# ============================================================
# POPC 膜平衡作业模板：EM → NVT(约束) → NPT(约束,半各向异性) → NPT(无约束,5ns)
# 用法：sbatch eq_mem.sbatch（在膜工作目录下提交）
# 注意：分区/资源参数为示例（Hygon DCU 集群），按目标集群修改：
#       -p 分区、--gres=dcu:1（GPU 集群）或删掉、--cpus-per-task、module load 行
# 前置：conf.gro + topol.top（含 posre 块）+ 四个 mdp 已在当前目录
# ============================================================
#SBATCH -J MEMEQ
#SBATCH -p wzhdtest
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=dcu:1
#SBATCH -t 72:00:00
#SBATCH -o memeq_%j.out
#SBATCH -e memeq_%j.err

source /etc/profile
module purge
module load apps/Gromacs-DCU2/2023.2-mpi-hpcx2.4.1-dtk23.10   # 按集群 GROMACS 模块修改
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}
export GMX_MAXBACKUPFILES=-1
cd "$SLURM_SUBMIT_DIR" || exit 1

G="gmx_mpi mdrun -v -ntomp $OMP_NUM_THREADS -nb gpu"   # 无 GPU 时用 "gmx mdrun -v"
set -e

echo "== EM =="
gmx_mpi grompp -f em.mdp -c conf.gro -p topol.top -o em.tpr -maxwarn 2
$G -deffnm em

echo "== NVT posres =="
gmx_mpi grompp -f nvt_mem.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr -maxwarn 2
$G -deffnm nvt

echo "== NPT posres semiiso =="
gmx_mpi grompp -f npt_mem.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr -maxwarn 2
$G -deffnm npt

echo "== NPT free semiiso 5ns =="
gmx_mpi grompp -f npt_mem_free.mdp -c npt.gro -t npt.cpt -p topol.top -o mdfree.tpr -maxwarn 2
$G -deffnm mdfree

echo MEM_EQ_DONE
