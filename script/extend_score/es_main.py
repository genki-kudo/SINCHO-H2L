from script.extend_score.constant import *
from script.extend_score.protein_check import *
from script.extend_score.pqr_and_involveatom import *
from script.extend_score.term1_hydro_norm import *
from script.extend_score.term2_maxdist import *
from script.extend_score.term3_asa import *
from script.extend_score.druggability_score import *
from script.extend_score.distance import *
import logging



def es_calc(all_pqr,fpocket_output, ligand, protein, weight, clusterdir):
    cal_res_lst = ds_calc(all_pqr, fpocket_output, protein)
    cal_res_lst = di_calc(ligand, cal_res_lst, clusterdir)
    #print(cal_res_lst)
    results = []
    for cal in cal_res_lst:
        result=[]
        result.append(cal[0])
        es = (weight*float(cal[5])) + ((1-weight)*(np.log10(1/float(cal[4]))))
        result.append(es)
        results.append(result)
    es_res = sorted(results, key=lambda x: x[1])
    print(es_res)

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

