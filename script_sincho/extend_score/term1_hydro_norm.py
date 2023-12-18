from script_sincho.basic.basic_func import *

def term1_tmp(pqr_inf):
    ################ calculate X1 mean_hyd_norm#######################
    lst = []
    pol_num = 0
    for i in pqr_inf:
        if i[12:16]=="APOL" or i[12:16]=="   C":
            l = []
            pol_num += 1
            l.append(vec_xyz(i))
            l.append(float(i[66:71]))
            lst.append(l)
    hyd = 0
    for i in lst:
        ok_num = 0
        for j in lst:
            if float(np.linalg.norm(i[0] - j[0])) <= float(i[1]+j[1]):
                ok_num += 1
        hyd += ok_num-1
    if hyd != 0 and pol_num != 0:
        mean_pol_dens = float(hyd/pol_num)
    else:
        mean_pol_dens = 0
    ################ calculate X1 mean_hyd_norm#######################

    return mean_pol_dens

def term1(cal_res_lst):
    ################ X1 normalization#######################
    norm_hyd = sorted(cal_res_lst, reverse=True, key=lambda x: x[1])
    max_hyd = norm_hyd[0][1]
    min_hyd = norm_hyd[-1][1]
    for each in cal_res_lst:
        each[1]= (each[1]-min_hyd)/(max_hyd - min_hyd)
    ################ X1 normalization#######################
    
    return cal_res_lst
