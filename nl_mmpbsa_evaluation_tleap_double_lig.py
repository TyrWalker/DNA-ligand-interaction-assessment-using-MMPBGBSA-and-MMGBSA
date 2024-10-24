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

    if type == 'rna':
        ligand1_frcmod = '/home/dejun/workspace/NLmmgbsa/results/%s/%s/OL3/%s_ligand1.frcmod' % (type, pdbid, pdbid)
        ligand2_frcmod = '/home/dejun/workspace/NLmmgbsa/results/%s/%s/OL3/%s_ligand2.frcmod' % (type, pdbid, pdbid)
        ligand1_mol2 = '/home/dejun/workspace/NLmmgbsa/results/%s/%s/OL3/%s_ligand1.mol2' % (type, pdbid, pdbid)
        ligand2_mol2 = '/home/dejun/workspace/NLmmgbsa/results/%s/%s/OL3/%s_ligand2.mol2' % (type, pdbid, pdbid)
        receptor1_pdb_file = '/home/dejun/workspace/NLmmgbsa/results/%s/%s/OL3/%s_receptor1.pdb' % (type, pdbid, pdbid)
    else:
        ligand1_frcmod = '/home/dejun/workspace/NLmmgbsa/results/%s/%s/ff99bsc0/%s_ligand1.frcmod' % (type, pdbid, pdbid)
        ligand2_frcmod = '/home/dejun/workspace/NLmmgbsa/results/%s/%s/ff99bsc0/%s_ligand2.frcmod' % (type, pdbid, pdbid)
        ligand1_mol2 = '/home/dejun/workspace/NLmmgbsa/results/%s/%s/ff99bsc0/%s_ligand1.mol2' % (type, pdbid, pdbid)
        ligand2_mol2 = '/home/dejun/workspace/NLmmgbsa/results/%s/%s/ff99bsc0/%s_ligand2.mol2' % (type, pdbid, pdbid)
        receptor1_pdb_file = '/home/dejun/workspace/NLmmgbsa/results/%s/%s/ff99bsc0/%s_receptor1.pdb' % (type, pdbid, pdbid)

    if os.path.exists(ligand1_frcmod) and os.path.exists(ligand2_frcmod) and os.path.exists(receptor1_pdb_file):
        pass

### leap1
    output = f'''# receptor1, ligand2, complex12
source {ff}
source leaprc.water.tip3p
source leaprc.gaff2
loadamberparams {ligand1_frcmod}
loadamberparams {ligand2_frcmod}
set default pbradii mbondi2
receptor = loadpdb {receptor1_pdb_file}
ligand1 = loadmol2 {ligand1_mol2}
ligand2 = loadmol2 {ligand2_mol2}
receptor1 = combine {{receptor ligand1}}
saveamberparm receptor1 {receptor1_prmtop} {receptor1_inpcrd}
savepdb receptor1 {receptor1_pdb}
saveamberparm ligand1 {ligand1_prmtop} {ligand1_inpcrd}
savepdb ligand1 {ligand1_pdb}
saveamberparm ligand2 {ligand2_prmtop} {ligand2_inpcrd}
savepdb ligand2 {ligand2_pdb}
complex12 = combine {{receptor ligand1 ligand2}}
saveamberparm complex12 {complex12_prmtop} {complex12_inpcrd}
savepdb complex12 {complex12_pdb}
solvateBox complex12 TIP3PBOX 8.0
addIons complex12 Cl- 0
addIons complex12 Na+ 0
saveamberparm complex12 {complex12_solv_prmtop} {complex12_solv_inpcrd}
savepdb complex12 {complex12_solv_pdb}
quit
'''
    leapin_file = os.path.join(output_dir, f'leap_{ff.split(".")[-1]}_1.in')
    leaplog_file = os.path.join(output_dir, f'leap_{ff.split(".")[-1]}_1.log')
    write_file(leapin_file, output)

    command = f'cd %s && tleap -f {leapin_file} > {leaplog_file}' % output_dir
    if not os.path.exists(complex12_solv_pdb):
        proc = subprocess.Popen(
            command,
            shell=True,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )
        p = proc.communicate()[0]

    else:
        print('tleap1 have been prepared before!')

### leap2
    output = f'''# receptor, ligand2, complex
source {ff}
source leaprc.water.tip3p
source leaprc.gaff2
loadamberparams {ligand1_frcmod}
loadamberparams {ligand2_frcmod}
set default pbradii mbondi2
receptor = loadpdb {receptor1_pdb_file}
ligand1 = loadmol2 {ligand1_mol2}
ligand2 = loadmol2 {ligand2_mol2}
saveamberparm receptor {receptor_prmtop} {receptor_inpcrd}
savepdb receptor {receptor_pdb}
complex = combine {{receptor ligand2}}
saveamberparm complex {complex_prmtop} {complex_inpcrd}
savepdb complex {complex_pdb}
solvateBox complex TIP3PBOX 8.0
addIons complex Cl- 0
addIons complex Na+ 0
saveamberparm complex {complex_solv_prmtop} {complex_solv_inpcrd}
savepdb complex {complex_solv_pdb}
quit
'''
    leapin_file = os.path.join(output_dir, f'leap_{ff.split(".")[-1]}_2.in')
    leaplog_file = os.path.join(output_dir, f'leap_{ff.split(".")[-1]}_2.log')
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
        print('tleap2 have been prepared before!')

    if os.path.exists(
            f"{os.path.join(output_dir, pdbid)}_complex_solv_leap_{ff.split('.')[-1]}.pdb") and os.path.exists(
            f"{os.path.join(output_dir, pdbid)}_complex12_solv_leap_{ff.split('.')[-1]}.pdb"):
        print('successful tleap')
        return (
            receptor_prmtop, receptor_inpcrd, receptor_pdb, receptor1_prmtop, receptor1_inpcrd, receptor1_pdb, ligand1_prmtop, ligand1_inpcrd, ligand1_pdb,
            ligand2_prmtop, ligand2_inpcrd, ligand2_pdb, complex_prmtop,
            complex_inpcrd, complex_pdb, complex_solv_prmtop, complex_solv_inpcrd, complex_solv_pdb, complex12_prmtop,
            complex12_inpcrd, complex12_pdb, complex12_solv_prmtop, complex12_solv_inpcrd, complex12_solv_pdb,)
    else:
        print('Failed')
        print(command)


def main():
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
