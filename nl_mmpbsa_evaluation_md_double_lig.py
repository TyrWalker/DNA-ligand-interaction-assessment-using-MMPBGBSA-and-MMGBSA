# module load cuda/11.2
# module load amber/22
import os
import mdtraj
import subprocess


def count_residue(pdb_file):
    t = mdtraj.load_pdb(pdb_file)
    return t.n_residues


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


def minimize_or_md(flag, solvation_models, output_dir, complex_prmtop, complex_incrd, complex_pdb):
    solvation_dir = os.path.join(os.path.join(os.path.join(output_dir, flag), solvation_models))
    if flag == 'no':
        dir_check(solvation_dir)
        print('Pass minimize_or_md')
        return complex_incrd

    elif flag == 'mini':
        dir_check(os.path.join(os.path.join(output_dir, flag)))
        if solvation_models == 'explicit':
            dir_check(solvation_dir)
            output = f'''Minimization input file in explicit solvent
 &cntrl
    imin=1,
    ntx=1,
    irest=0,
    maxcyc=5000,
    ncyc=2000,
    cut=10.0,
    ntpr=5000,
    ntwx=0,
    ntr=1,
    restraint_wt=1000,
    restraintmask=':1-{count_residue(complex_pdb)-2}',
/
'''
            minimize1_mdin = os.path.join(solvation_dir, 'minimization1.mdin')
            minimize1_inpcrd = os.path.join(solvation_dir, 'minimization1.inpcrd')
            minimize1_mdout = os.path.join(solvation_dir, 'minimization1.mdout')
            minimize1_mdinfo = os.path.join(solvation_dir, 'minimization1.mdinfo')
            minimize1_log = os.path.join(solvation_dir, 'minimization1.log')
            write_file(minimize1_mdin, output)

            output = f'''Minimization input file in explicit solvent
 &cntrl
    imin=1,
    ntx=1,
    irest=0,
    maxcyc=5000,
    ncyc=2000,
    cut=10.0,
    ntpr=5000,
    ntwx=0,
    ntr=1,
    restraint_wt=500.0,
    restraintmask=':1-{count_residue(complex_pdb)-2}',
/
'''
            minimize2_mdin = os.path.join(solvation_dir, 'minimization2.mdin')
            minimize2_inpcrd = os.path.join(solvation_dir, 'minimization2.inpcrd')
            minimize2_mdout = os.path.join(solvation_dir, 'minimization2.mdout')
            minimize2_mdinfo = os.path.join(solvation_dir, 'minimization2.mdinfo')
            minimize2_log = os.path.join(solvation_dir, 'minimization2.log')
            write_file(minimize2_mdin, output)

            output = f'''Minimization input file in explicit solvent
 &cntrl
    imin=1,
    ntx=1,
    irest=0,
    maxcyc=5000,
    ncyc=2000,
    cut=10.0,
    ntpr=5000,
    ntwx=0,
/
'''
            minimize3_mdin = os.path.join(solvation_dir, 'minimization3.mdin')
            minimize3_inpcrd = os.path.join(solvation_dir, 'minimization3.inpcrd')
            minimize3_mdout = os.path.join(solvation_dir, 'minimization3.mdout')
            minimize3_mdinfo = os.path.join(solvation_dir, 'minimization3.mdinfo')
            minimize3_log = os.path.join(solvation_dir, 'minimization3.log')
            write_file(minimize3_mdin, output)

            command = f'pmemd.cuda -O -i {minimize1_mdin} -p {complex_prmtop} -c {complex_incrd} -o {minimize1_mdout} -r {minimize1_inpcrd} -inf {minimize1_mdinfo} -ref {complex_incrd} > {minimize1_log}'

            if not os.path.exists(minimize1_inpcrd):
                proc = subprocess.Popen(
                    command,
                    shell=True,
                    stdin=subprocess.PIPE,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE
                )
                p = proc.communicate()[0]
            else:
                print('mini1 have been done before!')

            command = f'pmemd.cuda -O -i {minimize2_mdin} -p {complex_prmtop} -c {minimize1_inpcrd} -o {minimize2_mdout} -r {minimize2_inpcrd} -inf {minimize2_mdinfo} -ref {complex_incrd} > {minimize2_log}'

            if not os.path.exists(minimize2_inpcrd):
                proc = subprocess.Popen(
                    command,
                    shell=True,
                    stdin=subprocess.PIPE,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE
                )
                p = proc.communicate()[0]
            else:
                print('mini2 have been done before!')

            command = f'pmemd.cuda -O -i {minimize3_mdin} -p {complex_prmtop} -c {minimize2_inpcrd} -o {minimize3_mdout} -r {minimize3_inpcrd} -inf {minimize3_mdinfo} -ref {complex_incrd} > {minimize3_log}'

            if not os.path.exists(minimize3_inpcrd):
                proc = subprocess.Popen(
                    command,
                    shell=True,
                    stdin=subprocess.PIPE,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE
                )
                p = proc.communicate()[0]
            else:
                print('mini3 have been done before!')
                return minimize3_inpcrd

            if os.path.exists(minimize3_inpcrd):
                print('successful mini')
                return minimize3_inpcrd
            else:
                print('Failed')
                print(command)
                print(p)

        elif solvation_models in ['implicit1', 'implicit2', 'implicit5', 'implicit7', 'implicit8']:
            dir_check(solvation_dir)
            output = f'''Minimization input file in implicit solvent
 &cntrl
    imin=1,
    ntx=1,
    irest=0,
    maxcyc=5000,
    ncyc=2000,
    cut=1000,
    ntpr=5000,
    ntwx=0,
    igb={solvation_models[-1]},
/
'''
            minimize_mdin = os.path.join(solvation_dir, 'minimization.mdin')
            minimize_inpcrd = os.path.join(solvation_dir, 'minimization.inpcrd')
            minimize_mdout = os.path.join(solvation_dir, 'minimization.mdout')
            minimize_mdinfo = os.path.join(solvation_dir, 'minimization.mdinfo')
            minimize_log = os.path.join(solvation_dir, 'minimization.log')
            write_file(minimize_mdin, output)

            command = f'pmemd.cuda -O -i {minimize_mdin} -p {complex_prmtop} -c {complex_incrd} -o {minimize_mdout} -r {minimize_inpcrd} -inf {minimize_mdinfo} -ref {complex_incrd} > {minimize_log}'

            if not os.path.exists(minimize_inpcrd):
                proc = subprocess.Popen(
                    command,
                    shell=True,
                    stdin=subprocess.PIPE,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE
                )
                p = proc.communicate()[0]
            else:
                print('mini have been done before!')
                return minimize_inpcrd

            if os.path.exists(minimize_inpcrd):
                print('successful mini')
                return minimize_inpcrd
            else:
                print('Failed')
                print(command)
                print(p)

        else:
            print('solvation_models Wrong input')

    elif flag == 'md':
        dir_check(os.path.join(os.path.join(output_dir, flag)))
        if solvation_models == 'explicit':
            dir_check(solvation_dir)
            output = f'''Heat.in: 500ps of heat equilibration
 &cntrl
    imin = 0, irest = 0, ntx = 1,
    nstlim = 250000, dt = 0.002,
    ntc = 2, ntf = 2,
    cut = 10.0, ntb = 1,
    ntpr = 500, ntwx = 500,
    ntt = 3, gamma_ln = 2.0,
    tempi = 0.0, temp0 = 300.0,
    ntr = 1, restraintmask = ':1-{count_residue(complex_pdb)}&!@H=', restraint_wt = 500.0,
    nmropt = 1,
    iwrap = 1, ioutfm = 1,
/
 &wt 
    TYPE='TEMP0', istep1=0, istep2=250000,
    value1=0.1, value2=300.0 
/
 &wt TYPE='END' /
'''
            heat_mdin = os.path.join(solvation_dir, 'heat.mdin')
            heat_rst = os.path.join(solvation_dir, 'heat.rst')
            heat_mdout = os.path.join(solvation_dir, 'heat.mdout')
            heat_mdinfo = os.path.join(solvation_dir, 'heat.mdinfo')
            heat_log = os.path.join(solvation_dir, 'heat.log')
            write_file(heat_mdin, output)

            output = f'''Equil1.in: 1ns of constant pressure equilibration at 300K, BER
 &cntrl
    imin = 0, irest = 1, ntx = 5,
    nstlim = 500000, dt = 0.002,
    ntc = 2, ntf = 2, barostat = 1, 
    cut = 10.0, ntb = 2, ntp = 1, taup = 2.0,
    ntpr = 1000, ntwx = 1000, ntwr = 1000,
    ntt = 3, gamma_ln = 2.0,
    temp0 = 300.0,
    iwrap = 1, ioutfm = 1,
/
'''
            equil1_mdin = os.path.join(solvation_dir, 'equil1.mdin')
            equil1_rst = os.path.join(solvation_dir, 'equil1.rst')
            equil1_mdout = os.path.join(solvation_dir, 'equil1.mdout')
            equil1_mdinfo = os.path.join(solvation_dir, 'equil1.mdinfo')
            equil1_log = os.path.join(solvation_dir, 'equil1.log')
            write_file(equil1_mdin, output)

            output = f'''Equil2.in: 6ns of constant pressure equilibration at 300K, MC 
 &cntrl
    imin = 0, irest = 1, ntx = 5,
    nstlim = 3000000, dt = 0.002,
    ntc = 2, ntf = 2, barostat = 2, 
    cut = 10.0, ntb = 2, ntp = 1, taup = 2.0,
    ntpr = 1000, ntwx = 1000,
    ntt = 3, gamma_ln = 2.0,
    temp0 = 300.0,
    iwrap = 1, ioutfm = 1,
/
'''
            equil2_mdin = os.path.join(solvation_dir, 'equil2.mdin')
            equil2_rst = os.path.join(solvation_dir, 'equil2.rst')
            equil2_mdout = os.path.join(solvation_dir, 'equil2.mdout')
            equil2_mdinfo = os.path.join(solvation_dir, 'equil2.mdinfo')
            equil2_mdcrd = os.path.join(solvation_dir, 'equil2.mdcrd')
            equil2_log = os.path.join(solvation_dir, 'equil2.log')
            write_file(equil2_mdin, output)

            output = f'''trajin {equil2_mdcrd} 501 3000 10
trajout {equil2_mdcrd.split('.')[0]}_last5ns.mdcrd crd nobox
'''
            trajtake_in = os.path.join(solvation_dir, 'trajtake.in')
            trajtake_log = os.path.join(solvation_dir, 'trajtake.log')
            write_file(trajtake_in, output)

            if os.path.exists(
                    os.path.join(os.path.join(os.path.join(os.path.join(output_dir, "mini"), solvation_models)),
                                 "minimization3.inpcrd")):
                command = f'pmemd.cuda -O -i {heat_mdin} -p {complex_prmtop} -c {os.path.join(os.path.join(os.path.join(os.path.join(output_dir, "mini"), solvation_models)), "minimization3.inpcrd")} -o {heat_mdout} -r {heat_rst} -inf {heat_mdinfo} -ref {os.path.join(os.path.join(os.path.join(os.path.join(output_dir, "mini"), solvation_models)), "minimization3.inpcrd")} > {heat_log}'

                if not os.path.exists(heat_rst):
                    proc = subprocess.Popen(
                        command,
                        shell=True,
                        stdin=subprocess.PIPE,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE
                    )
                    p = proc.communicate()[0]
                else:
                    print('heat have been done before!')

                command = f'pmemd.cuda -O -i {equil1_mdin} -p {complex_prmtop} -c {heat_rst} -o {equil1_mdout} -r {equil1_rst} -inf {equil1_mdinfo} -ref {heat_rst} > {equil1_log}'

                if not os.path.exists(equil1_rst):
                    proc = subprocess.Popen(
                        command,
                        shell=True,
                        stdin=subprocess.PIPE,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE
                    )
                    p = proc.communicate()[0]
                else:
                    print('equil1 have been done before!')

                command = f'pmemd.cuda -O -i {equil2_mdin} -p {complex_prmtop} -c {equil1_rst} -o {equil2_mdout} -r {equil2_rst} -inf {equil2_mdinfo} -ref {equil1_rst} -x {equil2_mdcrd} > {equil2_log}'

                if not os.path.exists(equil2_rst):
                    proc = subprocess.Popen(
                        command,
                        shell=True,
                        stdin=subprocess.PIPE,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE
                    )
                    p = proc.communicate()[0]
                else:
                    print('equil2 have been done before!')
                command = f'cpptraj -p {complex_prmtop} < {trajtake_in} > {trajtake_log}'
                if not os.path.exists(f"{equil2_rst.split('.')[0]}_last5ns.mdcrd"):
                    proc = subprocess.Popen(
                        command,
                        shell=True,
                        stdin=subprocess.PIPE,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE
                    )
                    p = proc.communicate()[0]
                else:
                    print('trajintake have been done before!')
                    return f"{equil2_rst.split('.')[0]}_last5ns.mdcrd"

                if os.path.exists(f"{equil2_rst.split('.')[0]}_last5ns.mdcrd"):
                    print('successful md')
                    return f"{equil2_rst.split('.')[0]}_last5ns.mdcrd"
                else:
                    print('Failed')
                    print(command)
                    print(p)

            else:
                print('md Need mini first!')

        elif solvation_models in ['implicit1', 'implicit2', 'implicit5', 'implicit7', 'implicit8']:
            dir_check(solvation_dir)
            output = f'''Heat.in: 500ps of heat equilibration
 &cntrl
    imin = 0, irest = 0, ntx = 1,
    nstlim = 250000, dt = 0.002,
    ntc = 2, ntf = 2,
    cut = 1000, igb = {solvation_models[-1]},
    ntpr = 500, ntwx = 500,
    ntt = 3, gamma_ln = 2.0,
    tempi = 0.0, temp0 = 300.0,
    ntr = 1, restraintmask = ':1-{count_residue(complex_pdb)}&!@H=', restraint_wt = 500.0,
    nmropt = 1,
    ioutfm = 1,
/
 &wt 
    TYPE='TEMP0', istep1=0, istep2=250000,
    value1=0.1, value2=300.0 
/
 &wt TYPE='END' /
'''
            heat_mdin = os.path.join(solvation_dir, 'heat.mdin')
            heat_rst = os.path.join(solvation_dir, 'heat.rst')
            heat_mdout = os.path.join(solvation_dir, 'heat.mdout')
            heat_mdinfo = os.path.join(solvation_dir, 'heat.mdinfo')
            heat_log = os.path.join(solvation_dir, 'heat.log')
            write_file(heat_mdin, output)

            output = f'''Equil1.in: 1ns of constant pressure equilibration at 300K, BER
 &cntrl
    imin = 0, irest = 1, ntx = 5,
    nstlim = 500000, dt = 0.002,
    ntc = 2, ntf = 2, igb = {solvation_models[-1]},
    cut = 1000, taup = 2.0, barostat = 1, 
    ntpr = 1000, ntwx = 1000, ntwr = 1000,
    ntt = 3, gamma_ln = 2.0,
    temp0 = 300.0,
    ioutfm = 1,
/
'''
            equil1_mdin = os.path.join(solvation_dir, 'equil1.mdin')
            equil1_rst = os.path.join(solvation_dir, 'equil1.rst')
            equil1_mdout = os.path.join(solvation_dir, 'equil1.mdout')
            equil1_mdinfo = os.path.join(solvation_dir, 'equil1.mdinfo')
            equil1_log = os.path.join(solvation_dir, 'equil1.log')
            write_file(equil1_mdin, output)

            output = f'''Equil2.in: 6ns of constant pressure equilibration at 300K, MC
 &cntrl
    imin = 0, irest = 1, ntx = 5,
    nstlim = 3000000, dt = 0.002,
    ntc = 2, ntf = 2, igb = {solvation_models[-1]},
    cut = 1000, taup = 2.0, barostat = 2,
    ntpr = 1000, ntwx = 1000,
    ntt = 3, gamma_ln = 2.0,
    temp0 = 300.0,
    ioutfm = 1,
/
'''
            equil2_mdin = os.path.join(solvation_dir, 'equil2.mdin')
            equil2_rst = os.path.join(solvation_dir, 'equil2.rst')
            equil2_mdout = os.path.join(solvation_dir, 'equil2.mdout')
            equil2_mdinfo = os.path.join(solvation_dir, 'equil2.mdinfo')
            equil2_mdcrd = os.path.join(solvation_dir, 'equil2.mdcrd')
            equil2_log = os.path.join(solvation_dir, 'equil2.log')
            write_file(equil2_mdin, output)

            output = f'''trajin {equil2_mdcrd} 501 3000 10
trajout {equil2_mdcrd.split('.')[0]}_last5ns.mdcrd crd nobox
'''
            trajtake_in = os.path.join(solvation_dir, 'trajtake.in')
            trajtake_log = os.path.join(solvation_dir, 'trajtake.log')
            write_file(trajtake_in, output)

            if os.path.exists(
                    os.path.join(os.path.join(os.path.join(os.path.join(output_dir, "mini"), solvation_models)),
                                 "minimization.inpcrd")):
                command = f'pmemd.cuda -O -i {heat_mdin} -p {complex_prmtop} -c {os.path.join(os.path.join(os.path.join(os.path.join(output_dir, "mini"), solvation_models)), "minimization.inpcrd")} -o {heat_mdout} -r {heat_rst} -inf {heat_mdinfo} -ref {os.path.join(os.path.join(os.path.join(os.path.join(output_dir, "mini"), solvation_models)), "minimization.inpcrd")} > {heat_log}'

                if not os.path.exists(heat_rst):
                    proc = subprocess.Popen(
                        command,
                        shell=True,
                        stdin=subprocess.PIPE,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE
                    )
                    p = proc.communicate()[0]
                else:
                    print('heat have been done before!')

                command = f'pmemd.cuda -O -i {equil1_mdin} -p {complex_prmtop} -c {heat_rst} -o {equil1_mdout} -r {equil1_rst} -inf {equil1_mdinfo} -ref {heat_rst} > {equil1_log}'

                if not os.path.exists(equil1_rst):
                    proc = subprocess.Popen(
                        command,
                        shell=True,
                        stdin=subprocess.PIPE,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE
                    )
                    p = proc.communicate()[0]
                else:
                    print('equil1 have been done before!')

                command = f'pmemd.cuda -O -i {equil2_mdin} -p {complex_prmtop} -c {equil1_rst} -o {equil2_mdout} -r {equil2_rst} -inf {equil2_mdinfo} -ref {equil1_rst} -x {equil2_mdcrd} > {equil2_log}'

                if not os.path.exists(equil2_rst):
                    proc = subprocess.Popen(
                        command,
                        shell=True,
                        stdin=subprocess.PIPE,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE
                    )
                    p = proc.communicate()[0]
                else:
                    print('equil2 have been done before!')

                command = f'cpptraj -p {complex_prmtop} < {trajtake_in} > {trajtake_log}'
                if not os.path.exists(f"{equil2_rst.split('.')[0]}_last5ns.mdcrd"):
                    proc = subprocess.Popen(
                        command,
                        shell=True,
                        stdin=subprocess.PIPE,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE
                    )
                    p = proc.communicate()[0]
                else:
                    print('trajintake have been done before!')
                    return f"{equil2_rst.split('.')[0]}_last5ns.mdcrd"

                if os.path.exists(f"{equil2_rst.split('.')[0]}_last5ns.mdcrd"):
                    print('successful md')
                    return f"{equil2_rst.split('.')[0]}_last5ns.mdcrd"
                else:
                    print('Failed')
                    print(command)
                    print(p)

            else:
                print('md Need mini first!')

        else:
            print('solvation_models Wrong input')

    else:
        print('minimize_or_md Wrong input')
        return


def main():
    pdbids = ['3S4P', '4LVW', '4LVY', '2BE0']
    types = ['rna' for _ in pdbids]
    md_or_mini = 'mini'
    # md_or_mini = 'md'  # 进行短时间动力学模拟
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
        for force_field in ff_keys:
            if type == 'rna':
                solvation_models_ls = ['implicit8', 'explicit']
            else:
                solvation_models_ls = ['implicit7', 'explicit']
            for solvation_models in solvation_models_ls:
                output_dir = '/home/dejun/workspace/NLmmgbsa/results/%s/%s' % (type, pdbid)
                ff_dir = os.path.join(output_dir, force_field)
                receptor_prmtop, receptor_inpcrd, receptor_pdb, receptor1_prmtop, receptor1_inpcrd, receptor1_pdb, ligand1_prmtop, ligand1_inpcrd, ligand1_pdb, ligand2_prmtop, ligand2_inpcrd, ligand2_pdb, complex_prmtop, complex_inpcrd, complex_pdb, complex_solv_prmtop, complex_solv_inpcrd, complex_solv_pdb, complex12_prmtop, complex12_inpcrd, complex12_pdb, complex12_solv_prmtop, complex12_solv_inpcrd, complex12_solv_pdb = tleap_dl(force_fields[force_field], pdbid, ff_dir)
                # mini/md/no
                if solvation_models in ['no', 'implicit1', 'implicit2', 'implicit5', 'implicit7', 'implicit8']:
                    final_inpcrd = minimize_or_md(md_or_mini, solvation_models, ff_dir, complex12_prmtop, complex12_inpcrd, complex12_pdb)
                else:
                    final_inpcrd = minimize_or_md(md_or_mini, solvation_models, ff_dir, complex12_solv_prmtop, complex12_solv_inpcrd, complex12_pdb)


if __name__ == '__main__':
    main()
