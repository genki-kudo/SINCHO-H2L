from script_sincho.poc_vol_calc import poc_cbn_conv
from script_sincho.clus_mw_estimate import poc_mw_calc
from script_sincho.clus_mw_estimate import lig_poc_dist_list
from script_sincho.clus_mw_estimate import estimate_idealmw
from script_sincho.clus_mw_estimate import extend_startatom
from script_sincho.points_rest import r_candidate
from script_sincho.basic.basic_func import *
import sys
import json
import os
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem

args = sys.argv
data = args[1]

#absolute path of script file 
absscriptpath = os.path.abspath(__file__)
#stat-correlation-set
mw_from_vol = json.loads(open('/home/user01/GIT-repository/02-EXTEND-POINT/comp_stat_dataset/vol_mw_fitting.json').read())
mw_from_dist = json.loads(open('/home/user01/GIT-repository/02-EXTEND-POINT/comp_stat_dataset/dist_passmw_fitting.json').read())

txt = []
with open(data, 'r')as inp:
    for line in inp:
        lst = line.split()
        if lst[2] !='0':

            pocdir = lst[0]+'/p2c_output/cluster/cluster'+lst[2]+'.pqr'
            lig = lst[0]+'/'+lst[0]+'_hitmod.pdb'

            t_file(lst[0]+'/p2c_output/cluster/cluster'+lst[2]+'.pdb')
            num = 0
            with open(pocdir,'r')as inp2:
                for line in inp2:
                    num += 1
                    with open(lst[0]+'/p2c_output/cluster/cluster'+lst[2]+'.pdb','a')as pdb:
                        print(line[0:6]+'{:5}'.format(num)+'   C '+line[16:54]+'  1.00  0.00           C', file=pdb)
            mol = Chem.MolFromPDBFile(lst[0]+'/p2c_output/cluster/cluster'+lst[2]+'.pdb',proximityBonding=False)
            vol = AllChem.ComputeMolVolume(mol, gridSpacing=1.0)
            est_poc_mw = mw_from_vol[0]*vol+mw_from_vol[1]
            subprocess.call(['rm', lst[0]+'/p2c_output/cluster/cluster'+lst[2]+'.pdb'])

            with open(lig,'r')as ll:
                for atom in ll:
                    if lst[1]+' ' in atom:
                        xyz = np.array([float(atom[31:38]), float(atom[39:46]), float(atom[47:54])])
            dist_min = 1000000
            with open(pocdir,'r')as poc:
                for a in poc:
                    vec_a = vec_xyz(a)
                    dist_min = dist_cf(vec_a, xyz, dist_min)
            est_dist_mw = mw_from_dist[0]*dist_min+mw_from_dist[1]
            mw_pred = est_poc_mw+est_dist_mw
            txt.append(mw_pred)
            with open('txt.txt','a')as oo:
                print(mw_pred, est_poc_mw, est_dist_mw, file=oo)
        else:
            txt.append('0')
            with open('txt.txt','a')as oo:
                print('0', '0', '0',file=oo)




        
   