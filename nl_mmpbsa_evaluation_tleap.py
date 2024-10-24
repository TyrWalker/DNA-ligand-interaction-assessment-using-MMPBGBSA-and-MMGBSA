# module load amber/22
import os
import subprocess


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


def tleap(ff, pdbid, output_dir, type):
    # 创建新目录
    dir_check(output_dir)
    receptor_prmtop = f"{os.path.join(output_dir, pdbid)}_receptor_leap_{ff.split('.')[-1]}.prmtop"
    receptor_inpcrd = f"{os.path.join(output_dir, pdbid)}_receptor_leap_{ff.split('.')[-1]}.inpcrd"
    receptor_pdb = f"{os.path.join(output_dir, pdbid)}_receptor_leap_{ff.split('.')[-1]}.pdb"
    ligand_prmtop = f"{os.path.join(output_dir, pdbid)}_ligand_leap_{ff.split('.')[-1]}.prmtop"
    ligand_inpcrd = f"{os.path.join(output_dir, pdbid)}_ligand_leap_{ff.split('.')[-1]}.inpcrd"
    ligand_pdb = f"{os.path.join(output_dir, pdbid)}_ligand_leap_{ff.split('.')[-1]}.pdb"
    complex_prmtop = f"{os.path.join(output_dir, pdbid)}_complex_leap_{ff.split('.')[-1]}.prmtop"
    complex_inpcrd = f"{os.path.join(output_dir, pdbid)}_complex_leap_{ff.split('.')[-1]}.inpcrd"
    complex_pdb = f"{os.path.join(output_dir, pdbid)}_complex_leap_{ff.split('.')[-1]}.pdb"
    complex_solv_prmtop = f"{os.path.join(output_dir, pdbid)}_complex_solv_leap_{ff.split('.')[-1]}.prmtop"
    complex_solv_inpcrd = f"{os.path.join(output_dir, pdbid)}_complex_solv_leap_{ff.split('.')[-1]}.inpcrd"
    complex_solv_pdb = f"{os.path.join(output_dir, pdbid)}_complex_solv_leap_{ff.split('.')[-1]}.pdb"

    if type == 'rna':
        ligand_frcmod = '/home/dejun/workspace/NLmmgbsa/results/%s/%s/OL3/%s_ligand.frcmod' % (type, pdbid, pdbid)
        ligand_mol2 = '/home/dejun/workspace/NLmmgbsa/results/%s/%s/OL3/%s_ligand.mol2' % (type, pdbid, pdbid)
        receptor1_pdb_file = '/home/dejun/workspace/NLmmgbsa/results/%s/%s/OL3/%s_receptor1.pdb' % (type, pdbid, pdbid)
    else:
        ligand_frcmod = '/home/dejun/workspace/NLmmgbsa/results/%s/%s/ff99bsc0/%s_ligand.frcmod' % (type, pdbid, pdbid)
        ligand_mol2 = '/home/dejun/workspace/NLmmgbsa/results/%s/%s/ff99bsc0/%s_ligand.mol2' % (type, pdbid, pdbid)
        receptor1_pdb_file = '/home/dejun/workspace/NLmmgbsa/results/%s/%s/ff99bsc0/%s_receptor1.pdb' % (
        type, pdbid, pdbid)

    if os.path.exists(ligand_frcmod) and os.path.exists(receptor1_pdb_file):
        pass

    output = f'''source {ff}
source leaprc.water.tip3p
source leaprc.gaff2
loadamberparams {ligand_frcmod}
set default pbradii mbondi2
receptor=loadpdb {receptor1_pdb_file}
saveamberparm receptor {receptor_prmtop} {receptor_inpcrd}
savepdb receptor {receptor_pdb}
ligand=loadmol2 {ligand_mol2}
saveamberparm ligand {ligand_prmtop} {ligand_inpcrd}
savepdb ligand {ligand_pdb}
complex=combine {{receptor ligand}}
saveamberparm complex {complex_prmtop} {complex_inpcrd}
savepdb complex {complex_pdb}
solvateBox complex TIP3PBOX 8.0
addIons complex Cl- 0
addIons complex Na+ 0
saveamberparm complex {complex_solv_prmtop} {complex_solv_inpcrd}
savepdb complex {complex_solv_pdb}
quit
'''
    leapin_file = os.path.join(output_dir, f'leap_{ff.split(".")[-1]}.in')
    leaplog_file = os.path.join(output_dir, f'leap_{ff.split(".")[-1]}.log')
    write_file(leapin_file, output)

    command = f'cd %s && tleap -f {leapin_file} > {leaplog_file}' % output_dir
    if not os.path.exists(complex_solv_pdb):
        proc = subprocess.Popen(
            command,
            shell=True,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )
        p = proc.communicate()[0]
    else:
        print('tleap have been prepared before!')
        return (
            receptor_prmtop, receptor_inpcrd, receptor_pdb, ligand_prmtop, ligand_inpcrd, ligand_pdb, complex_prmtop,
            complex_inpcrd, complex_pdb, complex_solv_prmtop, complex_solv_inpcrd, complex_solv_pdb)

    if os.path.exists(
            f"{os.path.join(output_dir, pdbid)}_complex_solv_leap_{ff.split('.')[-1]}.pdb"):
        print('successful tleap')
        return (
            receptor_prmtop, receptor_inpcrd, receptor_pdb, ligand_prmtop, ligand_inpcrd, ligand_pdb, complex_prmtop,
            complex_inpcrd, complex_pdb, complex_solv_prmtop, complex_solv_inpcrd, complex_solv_pdb)
    else:
        print('Failed')
        print(command)


def main():
    pdbids = ['1ARJ', '1FYP', '1Q8N', '1UTS', '1UUD', '2AU4', '2F4S', '2F4T', '2F4U', '2KGP', '2KTZ', '2KU0', '2KX8',
              '2L94', '2LWK', '2MXS', '2N0J', '4LVZ', '4LW0', '6HAG', '1QD3',  '2O3W', '2O3X', '3NPN',  '1BYJ']
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
            force_fields = {'ff99bsc0': 'oldff/leaprc.ff99bsc0',
                            'bsc1': 'leaprc.DNA.bsc1',
                            'OL15': 'leaprc.DNA.OL15'}
        ff_keys = force_fields.keys()
        for ff in ff_keys:
            output_dir = '/home/dejun/workspace/NLmmgbsa/results/%s/%s' % (type, pdbid)
            ff_dir = os.path.join(output_dir, ff)
            tleap(force_fields[ff], pdbid, ff_dir, type)


if __name__ == '__main__':
    main()
