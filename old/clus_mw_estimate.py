import glob
import os
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem
from script.basic.basic_func import *
from Bio.PDB import *
from numpy import linalg as LA

def poc_mw_calc(estimate_poc_vol,res1,):
    estimate_poc_mw = {}
    #print(res1)
    for i in estimate_poc_vol:
        path = os.path.splitext(i)[0]
        name = (path.split('/')[-1])
        #name = name[4:]
        vol = estimate_poc_vol[i]
        mw = res1[0]*vol+res1[1]
        estimate_poc_mw[name] = mw
    return estimate_poc_mw

def lig_poc_dist_list(outdirname, ligand, atom_hydro_vec, cos_value):
    cls_list = glob.glob(outdirname+'/cluster*')
    dist_dictionary = {}
    dist_dict_re = {}
    dist_cos = {}
    for file in cls_list:
        path = os.path.splitext(file)[0]
        name = (path.split('/')[-1])
        with open(ligand, 'r')as lig:
            for lline in lig:
                dis_min = 50000
                aname = lline[13:16]
                vec_lig = vec_xyz(lline)
                if vec_lig != 'None':
                    with open(file, 'r')as cls:
                        for line in cls:
                            vec_as = vec_xyz(line)
                            if vec_as != 'None':
                                distance = float(np.linalg.norm(vec_as - vec_lig))
                                if distance <= dis_min:
                                    dis_min = distance
                                    atom_poc_vec = (vec_as - vec_lig)
                    inner = np.inner(atom_poc_vec, atom_hydro_vec[aname])
                    norm = LA.norm(atom_poc_vec) * LA.norm(atom_hydro_vec[aname])
                    cos = inner / norm
                    dist_dictionary[name+'_'+aname] = dis_min
                    if name+'_'+aname in dist_cos:
                        if dist_cos[name+'_'+aname] <= cos:
                            dist_cos[name+'_'+aname] = cos
                    else:
                        dist_cos[name+'_'+aname] = cos
    #print(dist_dictionary)
    #print(dist_cos)

    for i in cls_list:
        path = os.path.splitext(i)[0]
        name = (path.split('/')[-1])
        #print(name)
        anum = 0
        temp_dict = {}
        for j in dist_dictionary:
            if name in j:
                anum += 1
                temp_dict[j] = dist_dictionary[j]
        temp_dict = sorted(temp_dict.items(), key=lambda x:x[1])
        a2num = 0
        for k in temp_dict:
            a2num +=1
            #dist_dict_re[k[0]] = k[1]
            if a2num <=(anum/2)+0.5 and dist_cos[k[0]] >= cos_value:
                dist_dict_re[k[0]] = k[1]
    return dist_dict_re

def estimate_idealmw(estimate_pocmw, ligand, dist_dictionary, mw_from_dist):
    extend = {}
    for i in estimate_pocmw:
        with open(ligand,'r')as lig:
            for line in lig:
                if line[0:6] =='HETATM' and (i+'_'+line[13:16] in dist_dictionary.keys()):
                    #ideal_mw = estimate_pocmw[i]+(27*float(dist_dictionary[i+'_'+line[13:16]])/1.53)-27
                    #ideal_mw = estimate_pocmw[i]+(9*float(dist_dictionary[i+'_'+line[13:16]]))-13.77
                    ideal_mw = -35+estimate_pocmw[i]+(mw_from_dist[0]*float(dist_dictionary[i+'_'+line[13:16]])+mw_from_dist[1])
                    extend[i+'_'+line[13:16]] = ideal_mw
    return extend

def extend_startatom(extend, ideal_mw):
    rem_key = 'tmp'
    point = []
    pocket = []
    mw = []
    for rank in range(20):
        extend.pop(rem_key, None)
        mw_all = list(extend.values())
        best = getNearestValue(mw_all, ideal_mw)
        key = [k for k, v in extend.items() if v == best]
        rem_key = key[0]
        point.append(key[0].split('_')[1])
        pocket.append(key[0].split('_')[0])
        mw.append(best)
        print(rank+1, best, key[0].split('_')[0], key[0].split('_')[1])
    return point, pocket, mw





