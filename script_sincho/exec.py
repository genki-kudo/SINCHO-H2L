from script_sincho.input.check_input import *
from script_sincho.extend_score.druggability_score import *
from script_sincho.extend_score.es_main import *
from script_sincho.extend_score.sascore.sa_main import *
from script_sincho.atom_select.theta_distance_main import *
from script_sincho.mw.mw_main import *
from script_sincho.visualize.pse_gen import *
from script_sincho.logp.logp_main import *
import logging
import glob

def sincho_exec():

    input_list, logger = setting()
    checks(input_list, logger)
    print('input_list',input_list)

    #if input_list[6]=='True':
    logger.info("(A) P2C execution is skipped")
    # pocket check #
    all_pqr = [f for f in glob.glob(input_list[0]+'/*.pqr')]
    logger.info('Candidate pockets were selected>> '+input_list[0]+'/*.pqr')
    if len(all_pqr)==0:
        logger.error('No candidate pocket in '+input_list[0]+'! Prepare pocket pqr file!')
        exit()
    """
    else:
        #(A)P2C execution script#
        logger.info('(A) P2C execution')
    """

    ##### (B) pocket selection #####
    logger.info('(B) extend score is calculating ...')

    #########################################################################
    #230511_1 pocket-atom pairs select
    #230511_2 extend score calculation for each pair

    ligand = input_list[3]

    sa_ba = sa_base(ligand)
    #druggability calculation
    cal_res_lst = ds_calc(all_pqr, input_list[2], input_list[4])
    #print(cal_res_lst)

    es_list = []
    check_es = []
    check_es_terms = []

    for pqr in all_pqr:
        ds_res = [ s for s in cal_res_lst if pqr.split('/')[-1] in s]
        if ds_res:
            ds = ds_res[0][4]
            print(ds)
        else:
            print('ERROR! druggability score is not found')
        atom_dist_sort, pdb_num = theta_distance_calc(ligand, pqr)
        #print(atom_dist_sort)
        #[原子番号_原子ラベル, θ, distance]
        #print(pdb_num)
        #
        use_atoms = 2
        #
        if len(atom_dist_sort)!=0 and len(atom_dist_sort)>=use_atoms:
            for g in range(use_atoms):
                print(pqr,atom_dist_sort[g][0])
                sa_ave = sa_main(atom_dist_sort[g][0].split('_')[0],ligand)
                #print(pqr, atom_dist_sort[g][0], atom_dist_sort[g][2], ds ,sa_ave, sa_ba)
                print(pqr, atom_dist_sort[g][0], es_calc_each(atom_dist_sort[g][2], ds ,sa_ave-sa_ba))
                es_list.append([pqr, atom_dist_sort[g][0], es_calc_each(atom_dist_sort[g][2], ds ,sa_ave-sa_ba)])
                check_es.append([pqr.split('/')[-1], atom_dist_sort[g][0], round(atom_dist_sort[g][2],3), round(np.log10(1/ds),3), round(sa_ave-sa_ba,3), round(es_calc_each(atom_dist_sort[g][2], ds ,sa_ave-sa_ba),5)])
                check_es_terms.append([pqr.split('/')[-1], atom_dist_sort[g][0], round(0.1*float(atom_dist_sort[g][2]),3), round(0.3*float(np.log10(1/ds)),3), round(0.6*(float(sa_ave-sa_ba)),3), round(es_calc_each(atom_dist_sort[g][2], ds ,sa_ave-sa_ba),5)])

        elif len(atom_dist_sort)!=0 and len(atom_dist_sort)<use_atoms:
            for g in range(len(atom_dist_sort)):
                print(pqr,atom_dist_sort[g][0])
                sa_ave = sa_main(atom_dist_sort[g][0].split('_')[0],ligand)
                #print(pqr, atom_dist_sort[g][0], atom_dist_sort[g][2], ds ,sa_ave, sa_ba)
                print(pqr, atom_dist_sort[g][0], es_calc_each(atom_dist_sort[g][2], ds ,sa_ave-sa_ba))
                es_list.append([pqr, atom_dist_sort[g][0], es_calc_each(atom_dist_sort[g][2], ds ,sa_ave-sa_ba)])
                check_es.append([pqr.split('/')[-1], atom_dist_sort[g][0], round(atom_dist_sort[g][2],3), round(np.log10(1/ds),3), round(sa_ave-sa_ba,3), round(es_calc_each(atom_dist_sort[g][2], ds ,sa_ave-sa_ba),5)])
                check_es_terms.append([pqr.split('/')[-1], atom_dist_sort[g][0], round(0.1*float(atom_dist_sort[g][2]),3), round(0.3*float(np.log10(1/ds)),3), round(0.6*(float(sa_ave-sa_ba)),3), round(es_calc_each(atom_dist_sort[g][2], ds ,sa_ave-sa_ba),5)])

        else:
            print(pqr, 'None')
    
    #230511_3 sort es_minimum
    es_sorted = sorted(es_list, key=lambda x: x[2])
    check_sorted = sorted(check_es, key=lambda x: x[5])
    check_terms_sorted = sorted(check_es_terms, key=lambda x: x[5])
    with open("./sincho-output/check.txt",'a') as chs:
        print("cluster atom distance log(1/ds) delta_sa es",file=chs)
        for ch in check_sorted:
            print(" ".join(str(_) for _ in ch),file=chs)

    with open("./sincho-output/check_terms.txt",'a') as chs:
        print("cluster atom 0.1*distance 0.3*log(1/ds) 0.6*delta_sa es",file=chs)
        for ch in check_terms_sorted:
            print(" ".join(str(_) for _ in ch),file=chs)

    #print(es_sorted)

    #230511_4 mw and logP estimation
    logger.info('(D) molecular weight and logP are estimated ...')
    #logger.info('Correct Term of MW-Estimation ->'+str(input_list[10]))
    logger.info('##### RESULTS #####')
    # summing up molecular weight and adding sd-term  #
    ##20221209_correct term madadayo!
    if len(es_sorted)>=int(input_list[9]):
        for information in es_sorted[:int(input_list[9])]:
            idealmw = mw_calc(input_list[3], [information[0].split('/')[-1], information[1]], input_list[0], input_list[1],logger)
            ideallogp = logp_calc(input_list[4], [information[0].split('/')[-1], information[1]], input_list[0], input_list[1],logger)

            logger.info(information[0].split('/')[-1]+"_"+information[1]+" estimate-mw: "+ str(idealmw)+" estimate-logp: "+ str(ideallogp))

    elif len(es_sorted)==0:
        logger.info('ERROR! This complex cannot grow the compound!')
    else:
        for information in es_sorted:
            mw_calc(input_list[3], [information[0].split('/')[-1], information[1]], input_list[0], input_list[1],logger)



    

    logger.info('##### RESULTS #####')
    #########################################################################

    """
    # 1. extend score calc #
    es_res = es_calc(all_pqr, input_list[2], input_list[3], input_list[4], float(input_list[8]), input_list[0])
    
    # 2. pocket selection #
    sel_pocs = pocket_sel(es_res, int(input_list[9]),logger)
    #print(sel_pocs)

    ##### (C) start atom selection #####
    logger.info('(C) start atom is selecting ...')
    # 1. theta and distance check #
    poc_atom_lst = []
    pdb_num_list = []
    for poc in sel_pocs:
        lst_tmp = []
        atom_dist_sort, pdb_num = theta_distance_calc(input_list[3], input_list[0]+poc)
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
    #print(selected_poc_atom)
    ##### ATODE NAOSUYO! #####

    ##### (D) molecular weight estimation #####
    logger.info('(D) molecular weight is estimating ...')
    logger.info('Correct Term of MW-Estimation ->'+str(input_list[10]))
    logger.info('##### RESULTS #####')
    # summing up molecular weight and adding sd-term  #
    ##20221209_correct term madadayo!
    for poc_atom in selected_poc_atom:
        mw_calc(input_list[3], poc_atom, input_list[0], input_list[1],input_list[10],logger)
    logger.info('##### RESULTS #####')
    """
    ##### (Z) output format #####
    logger.info('SINCHO FINISHED')

    logger.info('PSE visualization...')
    visualize(input_list[4],input_list[3],check_terms_sorted, input_list[0], input_list[1], input_list[9])
    delete_file("subst.smi")




