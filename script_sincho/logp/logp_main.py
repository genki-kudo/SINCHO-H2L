#!/usr/bin/env python

from script_sincho.basic.basic_func import *
from script_sincho.comp_stat_dataset.stat_correlation import *
from script_sincho.logp.around_as import *
from script_sincho.logp.polar_ratio import *
import logging
import os

def logp_calc(protein, poc_and_atom, clusterdir, outputdir, logger):
    pocket = clusterdir+poc_and_atom[0]

    coefficients = pocenv_complogp()
    a, b = coefficients[0], coefficients[1]
    if os.path.isfile(outputdir+"/pocket_environment/pocenv"+pocket.split("/")[-1].split(".")[0][7:]+".pdb"):
        around_atoms=outputdir+"/pocket_environment/pocenv"+pocket.split("/")[-1].split(".")[0][7:]+".pdb"
        estimate_pocenv = polar_ratio(around_atoms)
    else:
        around_atoms = around_as(pocket, protein, 6, outputdir)
        estimate_pocenv = polar_ratio(around_atoms)
    
    ideallogp = (a*float(estimate_pocenv))+b
    #logger.info(poc_and_atom[0]+"_"+poc_and_atom[1]+" estimate-mw: "+ str(idealmw))
    print(estimate_pocenv)
    return ideallogp

