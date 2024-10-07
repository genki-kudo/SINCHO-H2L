#!/usr/bin/env python

from script_sincho.basic.basic_func import *
from script_sincho.atom_select.min_theta_dist import *

def theta_distance_calc(ligand, pqr):
    candidate_atom_and_hydro = {}
    pdb_num = {}
    for i in open(ligand,'r').readlines():
        if i[0:6]=='HETATM' or i[0:6]=='ATOM  ':
            bondh_coor = coordinate_bonding_hydrogen(ligand, i[12:16].replace(' ',''))
            #print("####################################################################")
            #print(i[12:16].replace(' ','') , bondh_coor)
            if len(bondh_coor)!=0:
                candidate_atom_and_hydro[str(int(i[6:11]))+'_'+i[12:16].replace(' ','')] = bondh_coor
                pdb_num[str(int(i[6:11]))] = i[12:16].replace(' ','')
    #print(candidate_atom_and_hydro)

    atom_theta_dist = min_theta_dist(candidate_atom_and_hydro, ligand, pqr)
    atom_theta_under90 = []
    for i in atom_theta_dist:
        if i[1]<90:
            atom_theta_under90.append(i)
    atom_dist_sort = sorted(atom_theta_under90, key=lambda x: x[2])
    #print(pqr, atom_dist_sort)

    return atom_dist_sort, pdb_num
    

    



def atom_reactivity(ligand, poc_atom_lst):
    return
    





        
