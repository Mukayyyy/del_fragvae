from .conversion import mol_to_smiles

def adock(receptor_input,
        smiles,
        ligand_name,
        center_x=7.750,
        center_y=-14.556,
        center_z=6.747,
        size_x=20,
        size_y=20,
        size_z=20,
        seed=None,
        cpu=1,
        lig_dir = './new_BAS/ligand_files/',
        out_dir = './new_BAS/output/',
        log_dir = './new_BAS/log/',
        conf_dir = './new_BAS/config/'):

    #Imports
    import os
    import subprocess
    import psutil
    import re
    from vina import Vina

    timeout_duration = 600
    score = 0

    #mkdir
    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(conf_dir, exist_ok=True)
    os.makedirs(log_dir, exist_ok=True)
    os.makedirs(lig_dir, exist_ok=True)

    ligand = lig_dir + ligand_name + '.pdbqt'
    output = out_dir + ligand_name + '_out.pdbqt'
    config = conf_dir + ligand_name + '.conf'
    log = log_dir + ligand_name + '_log.txt'

    v = Vina(sf_name='vina', cpu=cpu)
    v.set_receptor(receptor_input)

    #Convert smiles
    if not os.path.isfile(ligand):
        with subprocess.Popen('obabel -:"' + smiles + '" -O ' + ligand + ' -h --gen3d' + ' > /dev/null 2>&1', \
                shell=True, start_new_session=True) as proc:
            try:
                proc.wait(timeout=timeout_duration)
            except subprocess.TimeoutExpired:
                # os.killpg(os.getpgid(proc.pid), signal.SIGTERM)
                p = psutil.Process(proc.pid)
                p.kill()
    else:
        print(f'Ligand file: {ligand!r} already exists.')

    if not os.path.isfile(output):
        v.set_ligand_from_file(ligand)
        v.compute_vina_maps(center=[center_x, center_y, center_z], box_size=[size_x, size_y, size_z])
        v.dock(exhaustiveness=32, n_poses=20)
        v.write_poses(output, n_poses=5, overwrite=True)
    else:
        print(f'Output file: {output!r} already exists.')

    try:
        score = float("inf")
        with open(output, 'r') as f:
            for line in f.readlines():
                if "REMARK VINA RESULT" in line:
                    new_score = re.findall(r'([-+]?[0-9]*\.?[0-9]+)', line)[0]
                    break
                    #average_score = max(score, float(new_score))
        score = new_score
    except FileNotFoundError:
        print('Could not open output file')
        score = 0
    
    # energy = v.score()
    # score = energy[0]
    # print('score before minimization:', score)
    # energy_minimized = v.optimize()
    # score = energy_minimized[0]
    # v.write_pose(ligand, overwrite=True)
    # print('score after minimization:', score)

    return (score)

def calculateDockingScore(mol):
    protein_surface = './DATA/protein_files/6rqu.pdbqt'
    smi = mol_to_smiles(mol)
    ligand_name = smi.replace('(', '{').replace(')', '}')
    return adock(protein_surface, smi, ligand_name)
