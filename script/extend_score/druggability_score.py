import os
import numpy as np
import math
from script.basic.basic_func import *
from Bio.PDB import *
from subprocess import run

from script.extend_score.constant import *
from script.extend_score.protein_check import *
from script.extend_score.pqr_and_involveatom import *
from script.extend_score.term1_hydro_norm import *
from script.extend_score.term2_maxdist import *
from script.extend_score.term3_asa import *



bash=lambda x:run(x,shell=True)

def ds_calc(all_pqr, fpocket_output, protein):

    pro_lst = procheck(protein)

    cal_res_lst = []
    for pqr in all_pqr:
        cal_res_lst_tmp = []
        #pqr, term1, term2, term3
        cal_res_lst_tmp.append(os.path.basename(pqr))

        pqr_inf, involve, set_atm, surround = involve_set(pqr, fpocket_output, pro_lst)

        cal_res_lst_tmp.append(term1_tmp(pqr_inf))
        cal_res_lst_tmp.append(term2(pqr_inf))
        cal_res_lst_tmp.append(term3(pqr_inf, involve, set_atm, surround))
        cal_res_lst.append(cal_res_lst_tmp)
    cal_res_lst = term1(cal_res_lst)

    #summing up X1, X2, X3
    a0 = -9.5698768
    a1 = 7.479844
    a2 = 0.3696134
    a3 = -0.04671833
    #print(cal_res_lst)
    for f in cal_res_lst:
        d_score = 1.0/(1.0+math.exp(-(a0+(a1*f[1])+(a2*f[2])+(a3*f[3]))))
        f.append(d_score)

    return cal_res_lst
