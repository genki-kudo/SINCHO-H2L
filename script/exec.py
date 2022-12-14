from script.input.check_input import *
from script.extend_score.druggability_score import *
from script.extend_score.es_main import *
from script.atom_select.theta_distance_main import *
from script.mw.mw_main import *
import logging
import glob

def sincho_exec():

    input_list, logger = setting()
    checks(input_list, logger)

    if input_list[6]=='True':
        print('haha')
        ##### (A) P2C zikko shinasai #####
    
    # pocket check #
    all_pqr = [f for f in glob.glob(input_list[0]+'/*.pqr')]
    logger.info('Candidate pockets were selected>> '+input_list[0]+'/*.pqr')
    if len(all_pqr)==0:
        logger.error('No candidate pocket in '+input_list[0]+'! Prepare pocket pqr file!')
        exit()

    ##### (B) pocket selection #####
    logger.info('(B) extend score is calculating ...')
    # 1. extend score calc #
    es_res = es_calc(all_pqr, input_list[2], input_list[3], input_list[4], float(input_list[8]), input_list[0])
    # 2. pocket selection #
    sel_pocs = pocket_sel(es_res, int(input_list[9]),logger)

    ##### (C) start atom selection #####
    logger.info('(C) start atom is selecting ...')
    # 1. theta and distance check #
    poc_atom_lst = []
    for poc in sel_pocs:
        lst_tmp = []
        atom_dist_sort = theta_distance_calc(input_list[3], input_list[0]+poc)
        #print(atom_dist_sort)
        if len(atom_dist_sort) !=0:
            lst_tmp.append(poc)
            lst_tmp.append(atom_dist_sort[0][0])
            poc_atom_lst.append(lst_tmp)
    # 3. atom reactivity check #
    ##### ATODE TSUIKA SURUYO! #####
    atom_reactivity(input_list[3], poc_atom_lst)
    ##### ATODE TSUIKA SURUYO! #####
    # 4. atom selection and next input generation #
    ##### ATODE NAOSUYO! #####
    selected_poc_atom = poc_atom_lst[:10]
    print(selected_poc_atom)
    ##### ATODE NAOSUYO! #####

    ##### (D) molecular weight estimation #####
    logger.info('(D) molecular weight is estimating ...')
    # summing up molecular weight and adding sd-term  #
    ##20221209_correct term madadayo!
    for poc_atom in selected_poc_atom:
        mw_calc(input_list[3], poc_atom, input_list[0], input_list[1],logger)

    ##### (Z) output format #####
    logger.info('SINCHO FINISHED')




