from script_sincho.basic.basic_func import *

def term2(pqr_inf):
    ################ calculate X2 max_dist_of_alpha-spheres #######################
    max_dis = 0
    for i in pqr_inf:
        for j in pqr_inf:
            max_dis = dist_cf_max(vec_xyz(i), vec_xyz(j), max_dis)
    ################ calculate X2 max_dist_of_alpha-spheres #######################

    return max_dis