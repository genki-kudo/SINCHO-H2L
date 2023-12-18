import numpy as np
from script_sincho.basic.basic_func import *
from Bio.PDB import *
from script_sincho.extend_score.constant import *
from script_sincho.extend_score.protein_check import *
from script_sincho.extend_score.pqr_and_involveatom import *
from script_sincho.extend_score.term1_hydro_norm import *
from script_sincho.extend_score.term2_maxdist import *
from script_sincho.extend_score.term3_asa import *

def di_calc(ligand, cal_res_lst, clusterdir):

    ####################################################################################################
    #print(cal_res_lst)
    for clus in cal_res_lst:
        pqr = clus[0]
        cx,cy,cz = 0,0,0
        cn = 0
        for i in open(clusterdir+pqr,'r').readlines():
            if i[0:6]=='ATOM  ' or i[0:6]=='HETATM':
                cn += 1
                cx += float(i[30:38])
                cy += float(i[38:46])
                cz += float(i[46:54])
        cent = np.array([cx/cn, cy/cn, cz/cn])
        cent_novec = [cx/cn, cy/cn, cz/cn]

        ##suiso tuki atom to suiso no zisho##
        candidate_atom_and_hydro = {}
        for i in open(ligand,'r').readlines():
            if i[0:6]=='HETATM' or i[0:6]=='ATOM  ':
                #print(i[12:16])
                bondh_coor = coordinate_bonding_hydrogen(ligand, i[12:16].replace(' ',''))
                if len(bondh_coor)!=0:
                    candidate_atom_and_hydro[i[12:16].replace(' ','')] = bondh_coor
        #print(candidate_atom_and_hydro)
        for a, h in candidate_atom_and_hydro.items():
            dis_min = 50000
            heavyatom_coor = coordinate_of_atom(ligand, a)
            dtt = float(np.linalg.norm(cent - heavyatom_coor))
            
            if dtt < dis_min:
                for bh in h:
                    if theta_calc(bh, heavyatom_coor, cent_novec)<90:
                        dis_min =dtt
                        break
                ###kiwotsukete###
                #dis_min =dtt
                ###kiwotsukete###
        clus.append(dis_min)
    ####################################################################################################
    return cal_res_lst