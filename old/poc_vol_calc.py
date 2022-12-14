import numpy as np
from script.basic.basic_func import lat_gen
from script.basic.basic_func import t_file
import glob
import os
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem


def poc_vol_calc(outdirname):
    cls_list = glob.glob(outdirname+'/*.pqr')
    for cls in cls_list:
        path = os.path.splitext(cls)[0]
        name = (path.split('/')[-1])
        #print(name)
        lat_gen(cls, 'ATOM  ', outdirname+'/lat_'+name+'.pdb')

    estimate_pocvol = {}
    lat_list = glob.glob(outdirname+'/lat*')
    for lat in lat_list:
        with open(lat,'r')as file:
            estimate_pocvol[lat]=sum(1 for line in file)
    return estimate_pocvol

def poc_cbn_conv(outdirname):
    cls_list = glob.glob(outdirname+'/*.pqr')
    for cls in cls_list:
        path = os.path.splitext(cls)[0]
        name = (path.split('/')[-1])
        #print(name)
        num = 0
        t_file(outdirname+'/'+name+'.pdb')
        with open(cls,'r')as inp:
            for line in inp:
                num += 1
                with open(outdirname+'/'+name+'.pdb','a')as pdb:
                    print(line[0:6]+'{:5}'.format(num)+'   C '+line[16:54]+'  1.00  0.00           C', file=pdb)
    estimate_pocvol = {}
    cpoc_list = glob.glob(outdirname+'/*.pdb')
    for cpoc in cpoc_list:
        mol = Chem.MolFromPDBFile(cpoc,proximityBonding=False)
        estimate_pocvol[cpoc]=AllChem.ComputeMolVolume(mol, gridSpacing=1.0)
    return estimate_pocvol



