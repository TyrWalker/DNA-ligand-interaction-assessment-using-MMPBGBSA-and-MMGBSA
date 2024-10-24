# module load cuda/11.2
# module load amber/22
import os
import argparse
import subprocess
import multiprocessing as mp
import pandas as pd
import numpy as np
import random
try:
    import mdtraj
except:
    from prody import *


def count_residue(pdb_file):
    try:
        t = mdtraj.load_pdb(pdb_file)
        return t.n_residues
    except:
        t = parsePDB(pdb_file)
        return t.getResnums().max()


def write_file(output_file, outline):
    buffer = open(output_file, 'w')
    buffer.write(outline)
    buffer.close()
    if os.path.exists(output_file):
        print(f'{output_file}  Have been written!')
    else:
        print(f'{output_file}: Failed, not created!')


def dir_check(dir: str, create_dir: bool = True):
    if os.path.exists(dir):
        print(f'{dir}   Exists')
    else:
        if create_dir:
            os.system(f'mkdir -p {dir}')
            print(f'{dir}  Create Already!')
        else:
            print(f'{dir}  Not Exists Yet')


def tleap_dl(ff, pdbid, output_dir):
    receptor_prmtop = f"{os.path.join(output_dir, pdbid)}_receptor_leap_{ff.split('.')[-1]}.prmtop"
    receptor_inpcrd = f"{os.path.join(output_dir, pdbid)}_receptor_leap_{ff.split('.')[-1]}.inpcrd"
    receptor_pdb = f"{os.path.join(output_dir, pdbid)}_receptor_leap_{ff.split('.')[-1]}.pdb"
    receptor1_prmtop = f"{os.path.join(output_dir, pdbid)}_receptor1_leap_{ff.split('.')[-1]}.prmtop"
    receptor1_inpcrd = f"{os.path.join(output_dir, pdbid)}_receptor1_leap_{ff.split('.')[-1]}.inpcrd"
    receptor1_pdb = f"{os.path.join(output_dir, pdbid)}_receptor1_leap_{ff.split('.')[-1]}.pdb"
    ligand1_prmtop = f"{os.path.join(output_dir, pdbid)}_ligand1_leap_{ff.split('.')[-1]}.prmtop"
    ligand1_inpcrd = f"{os.path.join(output_dir, pdbid)}_ligand1_leap_{ff.split('.')[-1]}.inpcrd"
    ligand1_pdb = f"{os.path.join(output_dir, pdbid)}_ligand1_leap_{ff.split('.')[-1]}.pdb"
    ligand2_prmtop = f"{os.path.join(output_dir, pdbid)}_ligand2_leap_{ff.split('.')[-1]}.prmtop"
    ligand2_inpcrd = f"{os.path.join(output_dir, pdbid)}_ligand2_leap_{ff.split('.')[-1]}.inpcrd"
    ligand2_pdb = f"{os.path.join(output_dir, pdbid)}_ligand2_leap_{ff.split('.')[-1]}.pdb"
    complex_prmtop = f"{os.path.join(output_dir, pdbid)}_complex_leap_{ff.split('.')[-1]}.prmtop"
    complex_inpcrd = f"{os.path.join(output_dir, pdbid)}_complex_leap_{ff.split('.')[-1]}.inpcrd"
    complex_pdb = f"{os.path.join(output_dir, pdbid)}_complex_leap_{ff.split('.')[-1]}.pdb"
    complex_solv_prmtop = f"{os.path.join(output_dir, pdbid)}_complex_solv_leap_{ff.split('.')[-1]}.prmtop"
    complex_solv_inpcrd = f"{os.path.join(output_dir, pdbid)}_complex_solv_leap_{ff.split('.')[-1]}.inpcrd"
    complex_solv_pdb = f"{os.path.join(output_dir, pdbid)}_complex_solv_leap_{ff.split('.')[-1]}.pdb"
    complex12_prmtop = f"{os.path.join(output_dir, pdbid)}_complex12_leap_{ff.split('.')[-1]}.prmtop"
    complex12_inpcrd = f"{os.path.join(output_dir, pdbid)}_complex12_leap_{ff.split('.')[-1]}.inpcrd"
    complex12_pdb = f"{os.path.join(output_dir, pdbid)}_complex12_leap_{ff.split('.')[-1]}.pdb"
    complex12_solv_prmtop = f"{os.path.join(output_dir, pdbid)}_complex12_solv_leap_{ff.split('.')[-1]}.prmtop"
    complex12_solv_inpcrd = f"{os.path.join(output_dir, pdbid)}_complex12_solv_leap_{ff.split('.')[-1]}.inpcrd"
    complex12_solv_pdb = f"{os.path.join(output_dir, pdbid)}_complex12_solv_leap_{ff.split('.')[-1]}.pdb"
    return (
        receptor_prmtop, receptor_inpcrd, receptor_pdb, receptor1_prmtop, receptor1_inpcrd, receptor1_pdb,
        ligand1_prmtop, ligand1_inpcrd, ligand1_pdb,
        ligand2_prmtop, ligand2_inpcrd, ligand2_pdb, complex_prmtop,
        complex_inpcrd, complex_pdb, complex_solv_prmtop, complex_solv_inpcrd, complex_solv_pdb, complex12_prmtop,
        complex12_inpcrd, complex12_pdb, complex12_solv_prmtop, complex12_solv_inpcrd, complex12_solv_pdb)


def mmgbsa(flag, interior, output_dir, complex_solv_prmtop, complex_prmtop, receptor_prmtop, ligand_prmtop, inpcrd,
           complex_pdb, receptor_pdb, ligand_pdb, solvation_models):
    if flag == 'pbsa':
        output = f'''Input file for running PB
&general
  interval=1, verbose=1,
  receptor_mask = ":1-{count_residue(receptor_pdb)}"
  ligand_mask = ":{count_residue(complex_pdb)}"
/
&pb
  indi={interior}, radiopt=0,
/
'''
        dir_check(output_dir)
        mmgbsain_file = os.path.join(output_dir, f'mmpbsa_{interior}.in')
        mmgbsalog_file = os.path.join(output_dir, f'mmpbsa_{interior}.log')
        write_file(mmgbsain_file, output)

        if solvation_models in ['no', 'implicit1', 'implicit2', 'implicit5', 'implicit7', 'implicit8']:
            command = f'mpirun -np 48 MMPBSA.py -O -i {mmgbsain_file} -cp {complex_prmtop} -rp {receptor_prmtop} -lp {ligand_prmtop} -y {inpcrd} -o {os.path.join(output_dir, "FINAL_RESULTS_MMPBSA.dat")} -prefix {os.path.join(output_dir, "_MMPBSA_")}  > {mmgbsalog_file} && '
            command += f'cd %s && MMPBSA.py --clean > /dev/null' % output_dir  ## 清除临时文件

        elif solvation_models == 'explicit':
            command = f'mpirun -np 48 MMPBSA.py -O -i {mmgbsain_file} -sp {complex_solv_prmtop} -cp {complex_prmtop} -rp {receptor_prmtop} -lp {ligand_prmtop} -y {inpcrd} -o {os.path.join(output_dir, "FINAL_RESULTS_MMPBSA.dat")} -prefix {os.path.join(output_dir, "_MMPBSA_")} > {mmgbsalog_file} && '
            command += f'cd %s && MMPBSA.py --clean > /dev/null' % output_dir  ## 清除临时文件

        if not os.path.exists(f"{os.path.join(output_dir, 'FINAL_RESULTS_MMPBSA.dat')}"):
            proc = subprocess.Popen(
                command,
                shell=True,
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE
            )
            p = proc.communicate()[0]
        else:
            command = f'cd %s && MMPBSA.py --clean > /dev/null' % output_dir  ## 清除之前存在的临时文件
            os.system(command)
            print('mmgbsa have been done before!')
            return

        if os.path.exists(f"{os.path.join(output_dir, 'FINAL_RESULTS_MMPBSA.dat')}"):
            print('successful mmgbsa')
        else:
            print('Failed')
            print(command)
    elif flag in ['gbsa1', 'gbsa2', 'gbsa5', 'gbsa7', 'gbsa8']:
        output = f'''Input file for running GB
&general
  interval=1, verbose=1,
  receptor_mask = ":1-{count_residue(receptor_pdb)}"
  ligand_mask = ":{count_residue(complex_pdb)}"
/
&gb
  igb={flag[-1]},
/
&decomp
  idecomp=1, dec_verbose=0,
/
'''
        dir_check(output_dir)
        mmgbsain_file = os.path.join(output_dir, f'mmgbsa{flag[-1]}_{interior}.in')
        # mmgbsalog_file = os.path.join(output_dir, f'mmpbsa{flag[-1]}_{interior}.log')
        mmgbsalog_file = os.path.join(output_dir, f'mmgbsa{flag[-1]}_{interior}.log')
        write_file(mmgbsain_file, output)

        output = f'''File generated by MMPBSA.py 
&cntrl
  ntb=0, nsnb=99999, cut=999.0, imin=5,
  ncyc=0, igb={flag[-1]}, intdiel={interior}, gbsa=2, surften=0.0072, extdiel=80.0, idecomp=1,
  dec_verbose=0,
/
Residues considered as REC
RRES 1 {count_residue(receptor_pdb)}
END
Residues considered as LIG
LRES {count_residue(receptor_pdb) + 1} {count_residue(complex_pdb)}
END 
Residues to print
RES 1 {count_residue(complex_pdb)}
END 
END 
'''
        gb_decomp_com_file = os.path.join(output_dir, '_MMPBSA_gb_decomp_com.mdin')
        write_file(gb_decomp_com_file, output)

        output = f'''File generated by MMPBSA.py 
&cntrl
  ntb=0, nsnb=99999, cut=999.0, imin=5,
  ncyc=0, igb={flag[-1]}, intdiel={interior}, gbsa=2, surften=0.0072, extdiel=80.0, idecomp=1,
  dec_verbose=0,
/
Residues considered as REC
RRES 1 {count_residue(receptor_pdb)}
END
Residues to print
RES 1 {count_residue(receptor_pdb)}
END 
END 
'''
        gb_decomp_rec_file = os.path.join(output_dir, '_MMPBSA_gb_decomp_rec.mdin')
        write_file(gb_decomp_rec_file, output)

        output = f'''File generated by MMPBSA.py 
&cntrl
  ntb=0, nsnb=99999, cut=999.0, imin=5,
  ncyc=0, igb={flag[-1]}, intdiel={interior}, gbsa=2, surften=0.0072, extdiel=80.0, idecomp=1,
  dec_verbose=0,
/
Residues considered as LIG
RRES 1 {count_residue(ligand_pdb)}
END
Residues to print
RES 1 {count_residue(ligand_pdb)}
END 
END 
'''
        gb_decomp_lig_file = os.path.join(output_dir, '_MMPBSA_gb_decomp_lig.mdin')
        write_file(gb_decomp_lig_file, output)

        if solvation_models in ['no', 'implicit1', 'implicit2', 'implicit5', 'implicit7', 'implicit8']:
            command = f'mpirun -np 48 MMPBSA.py -O -i {mmgbsain_file} -cp {complex_prmtop} -rp {receptor_prmtop} -lp {ligand_prmtop} -y {inpcrd} -o {os.path.join(output_dir, "FINAL_RESULTS_MMPBSA.dat")} -prefix {os.path.join(output_dir, "_MMPBSA_")} -do {os.path.join(output_dir, "FINAL_DECOMP_MMPBSA.dat")} -use-mdins > {mmgbsalog_file} && '
            command += f'cd %s && MMPBSA.py --clean > /dev/null' % output_dir  ## 清除临时文件
        elif solvation_models == 'explicit':
            command = f'mpirun -np 48 MMPBSA.py -O -i {mmgbsain_file} -sp {complex_solv_prmtop} -cp {complex_prmtop} -rp {receptor_prmtop} -lp {ligand_prmtop} -y {inpcrd} -o {os.path.join(output_dir, "FINAL_RESULTS_MMPBSA.dat")} -prefix {os.path.join(output_dir, "_MMPBSA_")} -do {os.path.join(output_dir, "FINAL_DECOMP_MMPBSA.dat")} -use-mdins > {mmgbsalog_file} && '
            command += f'cd %s && MMPBSA.py --clean > /dev/null' % output_dir  ## 清除临时文件
        if not os.path.exists(f"{os.path.join(output_dir, 'FINAL_RESULTS_MMPBSA.dat')}"):
            proc = subprocess.Popen(
                command,
                shell=True,
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE
            )
            p = proc.communicate()[0]
        else:
            command = f'cd %s && MMPBSA.py --clean > /dev/null' % output_dir  ## 清除之前存在的临时文件
            os.system(command)
            print('mmgbsa have been done before!')
            return

        if os.path.exists(f"{os.path.join(output_dir, 'FINAL_RESULTS_MMPBSA.dat')}"):
            print('successful mmgbsa')
        else:
            print('Failed')
            print(command)
    else:
        print('mmgbsa Wrong input')
        return


def main(num_parts, part_id):
    if os.path.exists('/home/dejun/workspace/NLmmgbsa/scirpts/remained_mmgbsa_task.npy'):
        remained_dict = np.load('/home/dejun/workspace/NLmmgbsa/scirpts/remained_mmgbsa_task.npy',  allow_pickle=True).item()
    else:
        p1s, p2s, p3s, p4s, p5s, p6s, p7s, p8s, p9s, p10s, p11s, p12s = [], [], [], [], [], [], [], [], [], [], [], []
        pdbids = ['3S4P', '4LVW', '4LVY', '2BE0']
        types = ['rna' for _ in pdbids]
        for pdbid, type in zip(pdbids[:], types[:]):
            if type == 'rna':
                force_fields = {'OL3': 'leaprc.RNA.OL3',
                                'LJbb': 'leaprc.RNA.LJbb',
                                'YIL': 'leaprc.RNA.YIL',
                                'ff99bsc0': 'oldff/leaprc.ff99bsc0',
                                'ROC': 'leaprc.RNA.ROC',
                                'Shaw': 'leaprc.RNA.Shaw'}
            else:
                # DNA fource filed
                force_fields = {'ff99bsc0': 'oldff/leaprc.ff99bsc0',
                                'bsc1': 'leaprc.DNA.bsc1',
                                'OL15': 'leaprc.DNA.OL15'}
            ff_keys = force_fields.keys()
            # 确保pdbid在当前机器上！
            if os.path.exists('/home/dejun/workspace/NLmmgbsa/results/%s/%s' % (type, pdbid)):
                for force_field in ff_keys:
                    if type == 'rna':
                        solvation_models_ls = ['implicit8', 'explicit']
                    else:
                        solvation_models_ls = ['implicit7', 'explicit']
                    for solvation_models in solvation_models_ls:
                        output_dir = '/home/dejun/workspace/NLmmgbsa/results/%s/%s' % (type, pdbid)
                        # tleap
                        ff_dir = os.path.join(output_dir, force_field)
                        receptor_prmtop, receptor_inpcrd, receptor_pdb, receptor1_prmtop, receptor1_inpcrd, receptor1_pdb, ligand1_prmtop, ligand1_inpcrd, ligand1_pdb, ligand2_prmtop, ligand2_inpcrd, ligand2_pdb, complex_prmtop, complex_inpcrd, complex_pdb, complex_solv_prmtop, complex_solv_inpcrd, complex_solv_pdb, complex12_prmtop, complex12_inpcrd, complex12_pdb, complex12_solv_prmtop, complex12_solv_inpcrd, complex12_solv_pdb = tleap_dl(force_fields[force_field], pdbid, ff_dir)

                        for md_or_mini in ['md', 'no']:
                            if md_or_mini == 'md':
                                final_inpcrd = '/home/dejun/workspace/NLmmgbsa/results/%s/%s/%s/md/%s/equil2_last5ns.mdcrd' % (
                                type, pdbid, force_field, solvation_models)
                            else:
                                if solvation_models == 'explicit':  # 显示溶剂模型分三步能量最小化
                                    final_inpcrd = '/home/dejun/workspace/NLmmgbsa/results/%s/%s/%s/mini/%s/minimization3.inpcrd' % (
                                    type, pdbid, force_field, solvation_models)  # 用能量最小化的inpcrd文件进行mmgbsa打分
                                else:  # 隐式溶剂模型一步能量最小化
                                    final_inpcrd = '/home/dejun/workspace/NLmmgbsa/results/%s/%s/%s/mini/%s/minimization.inpcrd' % (
                                        type, pdbid, force_field, solvation_models)  # 用能量最小化的inpcrd文件进行mmgbsa打分

                            for free_energy_calculation in ['gbsa1', 'gbsa2', 'gbsa5', 'gbsa7', 'gbsa8', 'pbsa']:
                                for interior_dielectric_constants in [1, 2, 4, 8, 12, 16, 20]:
                                    # mmpbsa/mmgbsa
                                    solvation_dir = os.path.join(os.path.join(ff_dir, md_or_mini), solvation_models)
                                    mmgbsa_working_dir = os.path.join(solvation_dir, free_energy_calculation + '_' + str(
                                        interior_dielectric_constants))
                                    if not os.path.exists(f"{os.path.join(mmgbsa_working_dir, 'FINAL_RESULTS_MMPBSA.dat')}"):
                                        p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12 = free_energy_calculation, interior_dielectric_constants, mmgbsa_working_dir, complex12_solv_prmtop, complex12_prmtop, receptor1_prmtop, ligand2_prmtop, final_inpcrd, complex12_pdb, receptor1_pdb, ligand2_pdb, solvation_models
                                        p1s.append(p1)
                                        p2s.append(p2)
                                        p3s.append(p3)
                                        p4s.append(p4)
                                        p5s.append(p5)
                                        p6s.append(p6)
                                        p7s.append(p7)
                                        p8s.append(p8)
                                        p9s.append(p9)
                                        p10s.append(p10)
                                        p11s.append(p11)
                                        p12s.append(p12)


        # 对列表进行shuffle, 确保每个机器分配到的任务比较均衡
        p1s, p2s, p3s, p4s, p5s, p6s, p7s, p8s, p9s, p10s, p11s, p12s = np.array(p1s), np.array(p2s), np.array(p3s), np.array(p4s), np.array(p5s), np.array(p6s), np.array(p7s), np.array(p8s), np.array(p9s), np.array(p10s), np.array(p11s), np.array(p12s)
        ls = list(range(len(p1s)))
        random.shuffle(ls)
        p1s, p2s, p3s, p4s, p5s, p6s, p7s, p8s, p9s, p10s, p11s, p12s = list(p1s[ls]), list(p2s[ls]), list(p3s[ls]), list(p4s[ls]), list(p5s[ls]), list(p6s[ls]), list(p7s[ls]), list(p8s[ls]), list(p9s[ls]), list(p10s[ls]), list(p11s[ls]), list(p12s[ls])

        remained_dict = {'p1': p1s, 'p2': p2s, 'p3': p3s, 'p4': p4s, 'p5': p5s, 'p6': p6s,
                         'p7': p7s, 'p8': p8s, 'p9': p9s, 'p10': p10s, 'p11': p11s, 'p12': p12s}
        np.save('/home/dejun/workspace/NLmmgbsa/scirpts/remained_mmgbsa_task.npy', remained_dict)

    p1s, p2s, p3s, p4s, p5s, p6s, p7s, p8s, p9s, p10s, p11s, p12s = remained_dict['p1'], remained_dict['p2'], remained_dict['p3'], remained_dict['p4'], remained_dict['p5'], remained_dict['p6'], remained_dict['p7'], remained_dict['p8'], remained_dict['p9'], remained_dict['p10'], remained_dict['p11'], remained_dict['p12']
    print('the total number of remained mmgbsa  claculation task is: %s' % len(p1s))
    num_per_part = len(p1s) // num_parts
    if part_id != num_parts:
        p1s_sel, p2s_sel, p3s_sel, p4s_sel, p5s_sel, p6s_sel, p7s_sel, p8s_sel, p9s_sel, p10s_sel, p11s_sel, p12s_sel = p1s[num_per_part*(part_id-1):num_per_part*part_id], p2s[num_per_part*(part_id-1):num_per_part*part_id], p3s[num_per_part*(part_id-1):num_per_part*part_id], \
        p4s[num_per_part*(part_id-1):num_per_part*part_id], p5s[num_per_part*(part_id-1):num_per_part*part_id], p6s[num_per_part*(part_id-1):num_per_part*part_id], \
        p7s[num_per_part*(part_id-1):num_per_part*part_id], p8s[num_per_part*(part_id-1):num_per_part*part_id], p9s[num_per_part*(part_id-1):num_per_part*part_id], \
        p10s[num_per_part*(part_id-1):num_per_part*part_id], p11s[num_per_part*(part_id-1):num_per_part*part_id], p12s[num_per_part*(part_id-1):num_per_part*part_id]
    else:  # part_id == num_parts
        p1s_sel, p2s_sel, p3s_sel, p4s_sel, p5s_sel, p6s_sel, p7s_sel, p8s_sel, p9s_sel, p10s_sel, p11s_sel, p12s_sel = p1s[num_per_part*(part_id-1):], p2s[num_per_part*(part_id-1):], p3s[num_per_part*(part_id-1):], p4s[num_per_part*(part_id-1):], p5s[num_per_part*(part_id-1):], p6s[num_per_part*(part_id-1):], \
        p7s[num_per_part*(part_id-1):], p8s[num_per_part*(part_id-1):], p9s[num_per_part*(part_id-1):], \
        p10s[num_per_part*(part_id-1):], p11s[num_per_part*(part_id-1):], p12s[num_per_part*(part_id-1):]
    print('the selected number of remained mmgbsa  claculation task is: %s' % len(p1s_sel))

    for p1_sel, p2_sel, p3_sel, p4_sel, p5_sel, p6_sel, p7_sel, p8_sel, p9_sel, p10_sel, p11_sel, p12_sel in zip(p1s_sel, p2s_sel, p3s_sel, p4s_sel, p5s_sel, p6s_sel, p7s_sel, p8s_sel, p9s_sel, p10s_sel, p11s_sel, p12s_sel):
        mmgbsa(p1_sel, p2_sel, p3_sel, p4_sel, p5_sel, p6_sel, p7_sel, p8_sel, p9_sel, p10_sel, p11_sel, p12_sel)


if __name__ == '__main__':
    parser = argparse.ArgumentParser('mmgbsa  claculation')
    parser.add_argument("num_parts", type=int, default=4, help='how many parts')
    parser.add_argument("part_id", type=int, default=1, help="part_id")
    args = parser.parse_args()
    num_parts = args.num_parts
    part_id = args.part_id
    main(num_parts, part_id)
