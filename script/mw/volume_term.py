from script.basic.basic_func import *
import os
from rdkit import Chem
from rdkit.Chem import AllChem

def vol_calc(outputdir, pqr, clusterdir):
    convert_pdb = outputdir+'/'+str(os.path.basename(pqr))+'.pdb'
    t_file(outputdir+'/'+str(os.path.basename(pqr))+'.pdb')
    poc = clusterdir+pqr
    tmp_num = 0
    for i in open(poc,'r').readlines():
        tmp_num += 1
        with open(convert_pdb,'a')as out:
            print(i[0:6]+'{:5}'.format(tmp_num)+'   C '+i[16:54]+'  1.00  0.00           C', file=out)
    mol = Chem.MolFromPDBFile(convert_pdb,proximityBonding=False)
    vol = AllChem.ComputeMolVolume(mol, gridSpacing=1.0)

    return vol



