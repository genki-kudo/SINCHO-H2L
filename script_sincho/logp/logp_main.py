#!/usr/bin/env python

from script_sincho.basic.basic_func import *
from script_sincho.comp_stat_dataset.stat_correlation import *
from script_sincho.logp.around_as import *
import logging
import os

#logP予測関数：
def logp_calc(protein, poc_and_atom, clusterdir, outputdir, ligand, logger):
    pocket = clusterdir+poc_and_atom[0]

    rf_tree = logp_rftree_load()

    distances = [4.5, 5.0, 5.5, 6.0, 6.5, 7.0]

    df = pd.DataFrame(columns=["plogP_4.5","plogP_5.0","plogP_5.5","plogP_6.0","plogP_6.5","plogP_7.0",
                               "psasa_4.5","psasa_5.0","psasa_5.5","psasa_6.0","psasa_6.5","psasa_7.0",
                               "ppsa_4.5","ppsa_5.0","ppsa_5.5","ppsa_6.0","ppsa_6.5","ppsa_7.0"])
    
    for dist in distances:
        poc_logP, pocket_sasa, pocket_psa = around_atoms_search(pocket, protein, dist, outputdir, ligand)
        #データフレームの該当場所に追加
        df.loc[0, f"plogP_{str(dist)}"] = poc_logP
        df.loc[0, f"psasa_{str(dist)}"] = pocket_sasa
        df.loc[0, f"ppsa_{str(dist)}"] = pocket_psa
    x = df.iloc[:,0:]

    #df.to_csv("hahaha.csv", index=True)
    #print(rf_tree.predict(x))


    return rf_tree.predict(x)[0]
    


    

