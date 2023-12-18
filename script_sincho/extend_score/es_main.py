from script_sincho.extend_score.constant import *
from script_sincho.extend_score.protein_check import *
from script_sincho.extend_score.pqr_and_involveatom import *
from script_sincho.extend_score.term1_hydro_norm import *
from script_sincho.extend_score.term2_maxdist import *
from script_sincho.extend_score.term3_asa import *
from script_sincho.extend_score.druggability_score import *
from script_sincho.extend_score.distance import *
from script_sincho.atom_select.theta_distance_main import *
from script_sincho.extend_score.sascore.sa_main import *
import logging



def es_calc(all_pqr,fpocket_output, ligand, protein, weight, clusterdir):
    cal_res_lst = ds_calc(all_pqr, fpocket_output, protein)
    print("##################################druggability")
    """
    cal_res_lst = di_calc(ligand, cal_res_lst, clusterdir)
    print(cal_res_lst)
    with open('test.txt','a')as out:
        for i in cal_res_lst:
            print(i, file=out)
    results = []
    for cal in cal_res_lst:
        result=[]
        result.append(cal[0])
        #atom_dist_sort = theta_distance_calc(ligand, pqr)
        es = (weight*float(cal[5])) + ((1-weight)*(np.log10(1/float(cal[4]))))
        result.append(es)
        results.append(result)
    es_res = sorted(results, key=lambda x: x[1])
    #print(es_res)
    """
    ##########does not use 2023/05/16##########
    results = []
    for cal in cal_res_lst:
        for j in range(3):
            result=[]
            result.append(cal[0])
            atom_dist_sort = theta_distance_calc(ligand, clusterdir+cal[0])
            print(atom_dist_sort)
            cx,cy,cz = 0,0,0
            cn = 0
            for i in open(clusterdir+cal[0],'r').readlines():
                if i[0:6]=='ATOM  ' or i[0:6]=='HETATM':
                    cn += 1
                    cx += float(i[30:38])
                    cy += float(i[38:46])
                    cz += float(i[46:54])
            cent = np.array([cx/cn, cy/cn, cz/cn])
            heavyatom_coor = coordinate_of_atom(ligand, atom_dist_sort[j][0])
            dtt = float(np.linalg.norm(cent - heavyatom_coor))
            for i in open(ligand,'r'):
                if (i[0:6]=='HETATM'or i[0:6]=='ATOM  ')and i[12:16].replace(' ','')==atom_dist_sort[j][0]:
                    num = int(i[6:11])
            sa_ave = sa_main(num, ligand)
            es = (0.2*float(dtt)) + (0.4*(np.log10(1/float(cal[4]))))+(0.4*float(sa_ave))
            result.append(es)
            result.append(num)
            results.append(result)
    es_res = sorted(results, key=lambda x: x[1])
    for x in range(10):
        print(es_res[x])
    #print(es_res)
    ####################
    return es_res

def pocket_sel(es_res, topnum, logger):
    if len(es_res)>=topnum:
        es_top = es_res[:topnum]
    else:
        logger.warning('number of selected pockets:'+str(topnum)+' is over all pockets around ligand')
        es_top = es_res
    sel_pocs = []
    for i in es_top:
        sel_pocs.append(i[0])
    
    return sel_pocs

def es_calc_each(dist, ds, sa):
    es = (0.1*float(dist)) + (0.3*(np.log10(1/ds)))+(0.6*float(sa))
    return es


